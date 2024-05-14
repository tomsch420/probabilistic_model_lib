#pragma once

#include <vector>
#include <memory>
#include "probabilistic_model.h"
#include <cmath>

//FORWARD DECLARATIONS
class ProbabilisticCircuit;

typedef std::shared_ptr<ProbabilisticCircuit> ProbabilisticCircuitPtr_t;

class ProbabilisticCircuit : public ProbabilisticModel {
public:
    std::vector<ProbabilisticCircuitPtr_t> sub_circuits;

    /**
     * @return The string representing this component.
     */
    virtual std::string representation() const = 0;

    /**
     * Get the indices of the variables that are in both this circuit and the other circuit.
     *
     * The indices refer to the order of the variables in this circuit.
     *
     * @param other the other circuit.
     * @return the indices
     */
    std::vector<int> indices_of_intersection_with_other(const ProbabilisticCircuitPtr_t &other) const {
        std::vector<int> result;

        auto own_variables = get_variables();
        auto other_variables = other->get_variables();
        auto other_variables_iterator = other_variables->begin();
        int index = 0;
        for (auto const &variable : *own_variables) {
            if (variable == *other_variables_iterator){
                result.push_back(index);
                other_variables_iterator++;
            }
            index++;
        }
        return result;
    }

};

/**
 * Class for smooth sum units.
 */
class SmoothSumUnit : public ProbabilisticCircuit {
public:

    std::vector<double> weights;

    std::string representation() const override {
        return "+";
    }

    double likelihood(const FullEvidencePtr_t &event) const override {
        double sum = 0;
        auto current_weight = weights.begin();
        for (auto &sub_circuit: sub_circuits) {
            sum += *current_weight * sub_circuit->likelihood(event);
            current_weight++;
        }
        return sum;
    }

    double log_likelihood(const FullEvidencePtr_t &event) const override{
        double sum = 0;
        auto current_weight = weights.begin();
        for (auto &sub_circuit: sub_circuits) {
            sum += *current_weight * exp(sub_circuit->log_likelihood(event));
            current_weight++;
        }
        return log(sum);
    }

    void add_subcircuit(double weight, const ProbabilisticCircuitPtr_t &sub_circuit) {
        weights.push_back(weight);
        sub_circuits.push_back(sub_circuit);
    }

    template<typename... Args>
    static ProbabilisticCircuitPtr_t make_shared(Args &&... args) {
        return std::make_shared<SmoothSumUnit>(std::forward<Args>(args)...);
    };


    AbstractVariableSetPtr_t get_variables() const override {
        if (sub_circuits.empty()) {
            return make_shared_variable_set();
        }
        return sub_circuits[0]->get_variables();
    }

};

class DecomposableProductUnit : public ProbabilisticCircuit {
public:

    std::string representation() const override {
        return "*";
    }

    double log_likelihood(const FullEvidencePtr_t &event) const override {
        double product = 0;
        auto own_variables = get_variables();
        for (auto &sub_circuit: sub_circuits) {
            // get the indices of the subcircuit in this circuit
            auto indices = indices_of_intersection_with_other(sub_circuit);

            // only process the relevant part of this event
            auto sub_event = std::make_shared<FullEvidence>();
            for (auto index: indices) {
                sub_event->push_back(event->at(index));
            }
            product += sub_circuit->log_likelihood(sub_event);
        }
        return product;
    }

    void add_subcircuit(const ProbabilisticCircuitPtr_t &sub_circuit) {
        sub_circuits.push_back(sub_circuit);
    }

    AbstractVariableSetPtr_t get_variables() const override {
        auto result = make_shared_variable_set();
        for (auto &sub_circuit: sub_circuits) {
            auto sub_circuit_variables = sub_circuit->get_variables();
            result->insert(sub_circuit_variables->begin(), sub_circuit_variables->end());
        }
        return result;
    }

};
