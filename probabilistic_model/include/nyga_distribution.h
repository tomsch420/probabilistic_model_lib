#pragma once

#include <vector>
#include <memory>
#include <numeric>
#include "univariate.h"
#include "probabilistic_circuit.h"
#include "random_events/include/variable.h"
#include <optional>
#include <map>
#include <queue>

//FORWARD DECLARATIONS
class NygaDistribution;

class InductionStep;

// TYPEDEFS
typedef std::vector<double> WeightsVector;
typedef std::shared_ptr<WeightsVector> WeightsVectorPtr_t;
typedef std::vector<double> DataVector;
typedef std::shared_ptr<DataVector> DataVectorPtr_t;

typedef std::shared_ptr<NygaDistribution> NygaDistributionPtr_t;
typedef std::shared_ptr<InductionStep> InductionStepPtr_t;


class NygaDistribution : public DeterministicSumUnit {
public:

    double min_likelihood_improvement;
    size_t min_samples_per_quantile;
    ContinuousPtr_t variable;

    explicit NygaDistribution(const ContinuousPtr_t &variable, size_t min_samples_per_quantile = 1,
                              double min_likelihood_improvement = 0.1);

    double log_pdf(double value) const {
        return 0;
    }

    std::string distribution_representation() const {
        return "Ny";
    }

    template<typename... Args>
    static NygaDistributionPtr_t make_shared(Args &&... args) {
        return std::make_shared<NygaDistribution>(std::forward<Args>(args)...);
    };

    NygaDistributionPtr_t fit(const DataVectorPtr_t &data_p);

    NygaDistributionPtr_t fit_with_initial_induction_step(const InductionStepPtr_t &initial_induction_step);

};


/**
 *  Class for performing induction in the NygaDistributions.
 */
class InductionStep {
public:

    /**
     * A pointer to the entire vector of sorted and unique data points.
     */
    DataVectorPtr_t data_p;

    /**
     * A pointer to the weights of every unique sample in data.
     */
    WeightsVectorPtr_t weights_p;

    /**
     * The index of the first element of the data vector that is included in this step.
     */
    size_t begin_index;

    /**
     * The index of the first element of the data vector that is not included in this step.
     */
    size_t end_index;

    /**
     * The total number of samples in the data vector before it was made unique. This is bigger or equal to the size of
     * the data vector.
     */
    size_t total_number_of_samples;

    /**
     * The pointer to the Nyga Distribution to mount the quantile distributions into and read the parameters from.
     */
    NygaDistributionPtr_t nyga_distribution_p;


    /**
     * Construct an induction step and calculate the logarithm of the weights.
     * @param data_p The pointer to the data vector.
     * @param weights_p The pointer to the weights vector.
     * @param begin_index The index of the first element of the data vector that is included in this step.
     * @param end_index The index of the first element of the data vector that is not included in this step.
     * @param total_number_of_samples The total number of samples in the data vector before it was made unique.
     * @param nyga_distribution_p The pointer to the Nyga Distribution to mount the quantile distributions into and read the parameters from.
     */
    explicit InductionStep(const DataVectorPtr_t &data_p, const WeightsVectorPtr_t &weights_p, size_t begin_index,
                           size_t end_index, size_t total_number_of_samples,
                           const NygaDistributionPtr_t &nyga_distribution_p) : data_p(data_p),
                                                                               weights_p(weights_p),
                                                                               begin_index(begin_index),
                                                                               end_index(end_index),
                                                                               total_number_of_samples(
                                                                                       total_number_of_samples),
                                                                               nyga_distribution_p(
                                                                                       nyga_distribution_p) {
    }



    /**
     * Calculate the left connecting point given some beginning index.
     * @param index The index of the left datapoint.
     * @return the left connecting point
     */
    double left_connecting_point_from_index(size_t index) const;

    /**
     * @return the left connecting point.
     */
    double left_connecting_point() const;

    /**
    * Calculate the right connecting point given some beginning index.
    * @param index The index of the left datapoint.
    * @return the left connecting point
    */
    double right_connecting_point_from_index(size_t index) const;

    /**
     * @return the right connecting point.
     */
    double right_connecting_point() const;

    double data_size() const {
        return (double) data_p->size();
    }

    /**
     * Create a uniform distribution from the datapoint at `begin_index_` to the datapoint at `end_index_`.
     * @param begin_index_  The index of the first datapoint.
     * @param end_index_ The index of the last datapoint.
     * @return The uniform distribution pointer.
     */
    UniformDistributionPtr_t create_uniform_distribution_from_indices(size_t begin_index_, size_t end_index_) const;

    double log_likelihood_of_split(size_t split_index, double connecting_point) const;

    UniformDistributionPtr_t create_uniform_distribution() const;

    /**
     * Sum the logarithm of the weights from `begin_index_` to `end_index_`.
     * @param begin_index_ The index of the first weight.
     * @param end_index_ The index of the excluded last weight.
     * @return
     */
    double sum_log_weights_from_indices(size_t begin_index_, size_t end_index_) const;

    double sum_log_weights() const;

    /**
    * Sum the weights from `begin_index_` to `end_index_`.
    * @param begin_index_ The index of the first weight.
    * @param end_index_ The index of the excluded last weight.
    * @return
    */
    double sum_weights_from_indices(size_t begin_index_, size_t end_index_) const;

    double sum_weights() const;


    std::tuple<double, int> compute_best_split() const;

    /**
     * Construct the left induction step.
     * @param split_index The index of the split.
     * @return The left induction step.
     */
    InductionStepPtr_t construct_left_induction_step(size_t split_index) const;

    /**
     * Construct the right induction step.
     * @param split_index The index of the split.
     * @return The right induction step.
     */
    InductionStepPtr_t construct_right_induction_step(size_t split_index) const;


    [[maybe_unused]] std::optional<std::pair<InductionStepPtr_t, InductionStepPtr_t>> induce();

    template<typename... Args>
    static InductionStepPtr_t make_shared(Args &&... args) {
        return std::make_shared<InductionStep>(std::forward<Args>(args)...);
    };

};