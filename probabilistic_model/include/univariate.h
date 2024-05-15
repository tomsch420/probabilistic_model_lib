#pragma once

#include "probabilistic_model.h"
#include "probabilistic_circuit.h"
#include "variable.h"
#include "interval.h"
#include <cmath>
#include <utility>
#include <map>

//FORWARD DECLARATIONS
class UniformDistribution;
class DiracDeltaDistribution;

typedef std::shared_ptr<Interval<double>> ContinuousSupportPtr_t;
typedef std::shared_ptr<UniformDistribution> UniformDistributionPtr_t;
typedef std::shared_ptr<DiracDeltaDistribution> DiracDeltaDistributionPtr_t;


/**
 * Abstract Class for univariate distributions.
 */
class UnivariateDistribution : public ProbabilisticCircuit {
public:

    AbstractVariablePtr_t variable;

    /**
     * @return The support of the distribution.
     */
    virtual AbstractCompositeSetPtr_t get_support() const = 0;

    virtual std::string distribution_representation() const  = 0;

    AbstractVariableSetPtr_t get_variables() const override {
        auto result = make_shared_variable_set();
        result->insert(variable);
        return result;
    }

    std::string representation() const override {
        return *variable->name + " ~ " + distribution_representation();
    }

};


/**
 * Abstract Class for univariate, discrete distributions.
 */
class DiscreteDistribution : public UnivariateDistribution {
public:

    std::map<int, double> probabilities;

    explicit DiscreteDistribution(const AbstractVariablePtr_t &variable) {
        this->variable = variable;
        probabilities = std::map<int, double>();
    }

    DiscreteDistribution(const AbstractVariablePtr_t &variable, std::map<int, double> probabilities) {
        this->variable = variable;
        this->probabilities = std::move(probabilities);
    }

    double log_likelihood(const FullEvidencePtr_t &event) const override {
        return log(pmf(event->at(0)));
    }

    double pmf(int value) const {
        auto key = probabilities.find(value);
        if (key == probabilities.end()) {
            return 0;
        }
        return key->second;
    }

};

/**
 * Class for Symbolic Distributions.
 */
class SymbolicDistribution : public DiscreteDistribution {
public:

    SymbolicDistribution(const SymbolicPtr_t &variable, std::map<int, double> probabilities) : DiscreteDistribution(
            variable, std::move(probabilities)) {}

    AbstractCompositeSetPtr_t get_support() const override {
        auto all_elements = std::static_pointer_cast<Set>(variable->get_domain())->all_elements;
        auto result = variable->get_domain()->make_new_empty();
        for (auto &[value, probability]: probabilities) {
            if (probability == 0) {
                continue;
            }
            auto element = make_shared_set_element(value, all_elements);
            result = result->union_with(element);
        }
        return result;
    }

    std::string distribution_representation() const override{
        std::string result = "Nominal(";
        for (auto &[value, probability]: probabilities) {
            result += std::to_string(value) + ": " + std::to_string(probability) + ", ";
        }
        result += ")";
        return result;
    }

};


class IntegerDistribution : public DiscreteDistribution {
public:

    IntegerDistribution(const IntegerPtr_t &variable, std::map<int, double> probabilities) : DiscreteDistribution(
            variable, std::move(probabilities)) {}

    AbstractCompositeSetPtr_t get_support() const override {
        auto result = variable->get_domain()->make_new_empty();
        for (auto &[value, probability]: probabilities) {
            if (probability == 0) {
                continue;
            }
            auto element = singleton(value);
            result = result->union_with(element);
        }
        return result;
    }

    std::string distribution_representation() const override{
        std::string result = "Ordinal(";
        for (auto &[value, probability]: probabilities) {
            result += std::to_string(value) + ": " + std::to_string(probability) + ", ";
        }
        result += ")";
        return result;
    }
};


/**
 * Abstract Class for univariate, continuous distributions.
 */
class ContinuousDistribution : public UnivariateDistribution {
public:

    /**
     * The support of the distribution.
     */
    ContinuousSupportPtr_t support = reals();

    double log_likelihood(const FullEvidencePtr_t &event) const override {
        return log_pdf(event->at(0));
    }

    virtual double log_pdf(double value) const = 0;

    AbstractCompositeSetPtr_t get_support() const override {
        return support;
    }

};

/**
 * Class for the Dirac Delta distribution
 */
class DiracDeltaDistribution : public ContinuousDistribution {
public:
    double location;
    double density_cap;

    DiracDeltaDistribution(const AbstractVariablePtr_t &variable, double location,
                           double density_cap = std::numeric_limits<double>::infinity()) {
        this->variable = variable;
        this->location = location;
        this->density_cap = density_cap;
    }


    double pdf(double value) const {
        return value == location ? density_cap : 0;
    }

    double log_pdf(double value) const override {
        return log(pdf(value));
    }

    AbstractCompositeSetPtr_t get_support() const override {
        return singleton(location);
    }

    std::string distribution_representation() const override{
        return "δ(" + std::to_string(location) + ", " + std::to_string(density_cap) + ")";
    }

    template<typename... Args>
    static DiracDeltaDistributionPtr_t make_shared(Args &&... args) {
        return std::make_shared<DiracDeltaDistribution>(std::forward<Args>(args)...);
    };
};


/**
 * Class for the Uniform Distribution.
 */
class UniformDistribution : public ContinuousDistribution {
public:
    UniformDistribution(const ContinuousPtr_t &variable, const ContinuousSupportPtr_t &support) {
        this->variable = variable;
        this->support = support;
    }

    double pdf_value() const {
        return 1 / (support->upper() - this->support->lower());
    }

    double log_pdf(double value) const override {
        if (support->contains(value)) {
            return log(pdf_value());
        }
        return -std::numeric_limits<double>::infinity();
    }

    std::string distribution_representation() const override
    {
        return "U(" + *support->to_string() + ")";
    }

    template<typename... Args>
    static UniformDistributionPtr_t make_shared(Args &&... args) {
        return std::make_shared<UniformDistribution>(std::forward<Args>(args)...);
    };
};