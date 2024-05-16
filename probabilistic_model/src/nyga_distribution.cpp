//
// Created by tom_sch on 15.05.24.
//
#include <include/nyga_distribution.h>

NygaDistribution::NygaDistribution(const ContinuousPtr_t &variable, size_t min_samples_per_quantile,
                                   double min_likelihood_improvement) {

    this->variable = variable;
    this->min_samples_per_quantile = min_samples_per_quantile;
    this->min_likelihood_improvement = min_likelihood_improvement;
}

NygaDistributionPtr_t  NygaDistribution::fit(const DataVectorPtr_t &data_p) {

    auto result = NygaDistribution::make_shared(variable, min_samples_per_quantile, min_likelihood_improvement);

    // sort the data
    std::sort(data_p->begin(), data_p->end());

    // create a map that maps from unique values to the frequency
    auto frequencies = std::map<double, int>();

    // initialize the previous element to a value that is not in the data
    double previous_element = data_p->at(0) - 1;

    // for every data point, update the frequencies
    for (auto &data_point: *data_p) {
        if (data_point == previous_element) {
            frequencies[data_point]++;
        } else {
            frequencies[data_point] = 1;
            previous_element = data_point;
        }
    }

    if (frequencies.size() == 1) {
        auto distribution = DiracDeltaDistribution::make_shared(variable, data_p->at(0));
        result->add_subcircuit(1., distribution);
        return result;
    }

    // create the data and weights vector
    auto weights_p = new WeightsVector();
    auto sorted_unique_data = new DataVector();
    sorted_unique_data->reserve(frequencies.size());
    weights_p->reserve(frequencies.size());

    // write data and weights to new vectors
    for (auto &[value, weight]: frequencies) {
        weights_p->emplace_back(log(frequencies[value]));
        sorted_unique_data->emplace_back(value);
    }

    auto initial_induction_step = InductionStep::make_shared(data_p, weights_p, 0, frequencies.size(), result);
    result = fit_with_initial_induction_step(initial_induction_step);

    // clean up
    delete weights_p;
    delete sorted_unique_data;

    return result;
}

NygaDistributionPtr_t  NygaDistribution::fit_with_initial_induction_step(const InductionStepPtr_t &initial_induction_step) {
    auto induction_steps = std::queue<InductionStepPtr_t>();
    induction_steps.push(initial_induction_step);

    while (!induction_steps.empty()) {
        auto induction_step = induction_steps.front();
        induction_steps.pop();
        auto result = induction_step->induce();
        if (result.has_value()) {
            auto [left, right] = result.value();
            induction_steps.push(left);
            induction_steps.push(right);
        }

    }

    return initial_induction_step->nyga_distribution_p;

}

double InductionStep::left_connecting_point_from_index(size_t index) const {
    if (index > 0) {
        return (data_p->at(index - 1) + data_p->at(index)) / 2;
    }
    return data_p->at(index);
}

double InductionStep::left_connecting_point() const {
    return left_connecting_point_from_index(begin_index);
}

double InductionStep::right_connecting_point_from_index(size_t index) const {
    if (index < data_p->size()) {
        return (data_p->at(index - 1) + data_p->at(index)) / 2;
    }
    return data_p->at(index - 1);
}

double InductionStep::right_connecting_point() const {
    return right_connecting_point_from_index(end_index);
}

UniformDistributionPtr_t
InductionStep::create_uniform_distribution_from_indices(size_t begin_index_, size_t end_index_) const {

    IntervalPtr_t<double> interval;

    if (end_index_ == data_p->size()) {
        interval = closed(left_connecting_point_from_index(begin_index_),
                          right_connecting_point_from_index(end_index_));
    } else {
        interval = closed_open(left_connecting_point_from_index(begin_index_),
                               right_connecting_point_from_index(end_index_));
    }
    auto variable = std::static_pointer_cast<Continuous>(nyga_distribution_p->variable);
    return UniformDistribution::make_shared(variable, interval);
}

UniformDistributionPtr_t InductionStep::create_uniform_distribution() const {
    return create_uniform_distribution_from_indices(begin_index, end_index);
}


std::tuple<double, int> InductionStep::compute_best_split() const {
    double maximum_log_likelihood = -std::numeric_limits<double>::infinity();
    int best_split_index = -1;

    auto right_connecting_point_ = right_connecting_point();
    auto left_connecting_point_ = left_connecting_point();

    for (int split_index = (int) (begin_index + nyga_distribution_p->min_samples_per_quantile);
         split_index < end_index - nyga_distribution_p->min_samples_per_quantile + 1; split_index++) {


        // Calculate the log likelihood of the left side
        auto log_likelihood_left = log_likelihood_of_split(split_index, left_connecting_point_);

        // Calculate the log likelihood of the right side
        auto log_likelihood_right = log_likelihood_of_split(split_index, right_connecting_point_);

        // Calculate the average log-likelihood
        auto average_likelihood = (log_likelihood_left + log_likelihood_right);

        // update the maximum likelihood and the best split index
        if (average_likelihood > maximum_log_likelihood) {
            maximum_log_likelihood = average_likelihood;
            best_split_index = split_index;
        }

    }
    return std::make_tuple(maximum_log_likelihood, best_split_index);
}

InductionStepPtr_t InductionStep::construct_left_induction_step(size_t split_index) const {
    return InductionStep::make_shared(data_p, log_weights_p, begin_index, split_index, nyga_distribution_p);
}

InductionStepPtr_t InductionStep::construct_right_induction_step(size_t split_index) const {
    return InductionStep::make_shared(data_p, log_weights_p, split_index, end_index, nyga_distribution_p);
}

std::optional<std::pair<InductionStepPtr_t, InductionStepPtr_t>> InductionStep::induce() {
    double summed_weights = sum_weights();
    double log_pdf = -log(right_connecting_point() - left_connecting_point());
    double log_likelihood_without_split = log_pdf + summed_weights;

    auto [best_log_likelihood, best_split_index] = compute_best_split();
    if (best_log_likelihood > log_likelihood_without_split +  log(1. + nyga_distribution_p->min_likelihood_improvement)) {
        return std::make_pair(construct_left_induction_step(best_split_index),
                              construct_right_induction_step(best_split_index));
    }

    // create uniform distribution and mount it into the nyga distribution
    auto distribution = create_uniform_distribution();
    nyga_distribution_p->add_subcircuit(summed_weights, distribution);
    return std::nullopt;
}

double InductionStep::sum_weights() const {
    return sum_weights_from_indices(begin_index, end_index);
}

double InductionStep::sum_weights_from_indices(size_t begin_index_, size_t end_index_) const {
    double result = 0;
    for (size_t i = begin_index_; i < end_index_; i++) {
        result += (*log_weights_p)[i];
    }
    return result;
}

double InductionStep::log_likelihood_of_split(size_t split_index, double connecting_point) const {
    // Calculate the split value
    auto split_value = (data_p->at(split_index - 1) + data_p->at(split_index)) / 2;

    // Calculate the density and the side of the split
    double density = split_value - connecting_point;
    bool is_left = density < 0;

    // Sum the logarithmic weights that are inside this split
    double log_weight_sum_of_split = is_left ? sum_weights_from_indices(begin_index, split_index)
                                         : sum_weights_from_indices(split_index, end_index);


    // Calculate the log density assuming a uniform distribution and accumulated it with the number of samples
    auto log_density = -log(fabs(density));

    // return weighted log likelihood
    return log_density + log_weight_sum_of_split;
}

const double InductionStep::number_of_samples() const {
    return (double) end_index - begin_index;
}
