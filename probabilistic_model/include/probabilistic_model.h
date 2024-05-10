#pragma once

#include <set>
#include "variable.h"

using AbstractVariableSet_Ptr_t = std::shared_ptr<std::set<AbstractVariablePtr_t>>;

class ProbabilisticModel {
public:
    ProbabilisticModel() = default;
    AbstractVariableSet_Ptr_t variables = nullptr;

};