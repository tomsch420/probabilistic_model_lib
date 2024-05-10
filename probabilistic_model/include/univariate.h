#pragma once
#include "probabilistic_model.h"
#include "variable.h"
#include "interval.h"
#include <cmath>

class UnivariateDistribution : public ProbabilisticModel {
};

class ContinuousDistribution : public UnivariateDistribution {
public:

    /**
     * The support of the distribution.
     */
    IntervalPtr_t support;

    double log_likelihood(FullEvidencePtr_t event) const override {
        return log_pdf(event->at(0));
    }

    virtual double log_pdf(double value) const = 0;

};

class UniformDistribution : public ContinuousDistribution {

    UniformDistribution(ContinuousPtr_t variable, IntervalPtr_t support) {
        this->variables = make_shared_variable_set({std::static_pointer_cast<AbstractVariable>(variable));
        this->support = support;
    }

    double pdf_value() const {
      return  1 / (support->upper() - this->support->lower());
    }

    double log_pdf(double value) const override {
        if (support->contains((float) value)) {
            return log(pdf_value());
        }
        return -INFINITY;
    }
};