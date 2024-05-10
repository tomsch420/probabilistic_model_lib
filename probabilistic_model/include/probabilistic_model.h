#pragma once

#include <set>
#include "variable.h"
#include "sigma_algebra.h"

// TYPEDEFS
typedef std::shared_ptr<std::set<AbstractVariablePtr_t>> AbstractVariableSetPtr_t;
typedef std::vector<float> FullEvidence;
typedef std::shared_ptr<FullEvidence> FullEvidencePtr_t;


template<typename... Args>
std::shared_ptr<std::set<AbstractVariablePtr_t>> make_shared_variable_set(Args &&... args) {
    return std::make_shared<std::set<AbstractVariablePtr_t>>(std::forward<Args>(args)...);
}

class ProbabilisticModel {
public:
    AbstractVariableSetPtr_t variables = nullptr;

    virtual double log_likelihood(FullEvidencePtr_t event) const = 0;

};