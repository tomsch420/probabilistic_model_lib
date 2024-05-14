#pragma once

#include <set>
#include "variable.h"
#include "sigma_algebra.h"
#include <cmath>
#include <utility>


// TYPEDEFS
typedef std::shared_ptr<std::set<AbstractVariablePtr_t, PointerLess<AbstractVariablePtr_t>>> AbstractVariableSetPtr_t;
typedef std::vector<double> FullEvidence;
typedef std::shared_ptr<FullEvidence> FullEvidencePtr_t;


template<typename... Args>
std::shared_ptr<std::set<AbstractVariablePtr_t, PointerLess<AbstractVariablePtr_t >>>
make_shared_variable_set(Args &&... args) {
    return std::make_shared<std::set<AbstractVariablePtr_t, PointerLess<AbstractVariablePtr_t>>>(
            std::forward<Args>(args)...);
}

class ProbabilisticModel {
public:

    /**
     * @return The variables of the model.
     */
    virtual AbstractVariableSetPtr_t get_variables() const = 0;

    /**
     * The likelihood of an event.
     *
     * This method has by default calls the log_likelihood method. Hence one of the two has to be overloaded.
     *
     * @param event the event.
     * @return the likelihood of the event.
     */
    virtual double likelihood(const FullEvidencePtr_t &event) const {
        return exp(log_likelihood((event)));
    };

    /**
     * The log-likelihood of an event.
     *
     * This method has by default calls the likelihood method. Hence one of the two has to be overloaded.
     *
     * @param event the event.
     * @return the likelihood of the event.
     */
    virtual double log_likelihood(const FullEvidencePtr_t &event) const {
        return log(likelihood(event));
    };

};