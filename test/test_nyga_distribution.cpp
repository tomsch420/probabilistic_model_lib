#include <random>
#include "gtest/gtest.h"
#include "nyga_distribution.h"
#include "variable.h"

class NygaDistributionTest : public testing::Test {
public:
    ContinuousPtr_t variable_x = make_shared_continuous("x");

    DataVectorPtr_t data_p = std::make_shared<DataVector>(DataVector{1, 2, 3, 4, 7, 9});
    WeightsVectorPtr_t weights_p = std::make_shared<WeightsVector>(
            WeightsVector{1. / 6., 1. / 6. ,1. / 6. ,1. / 6. ,1. / 6. ,1. / 6.});
    NygaDistributionPtr_t model = NygaDistribution::make_shared(variable_x, 1, 0.01);
    InductionStep induction_step = InductionStep(data_p, weights_p, 0, 6, 6, model);

};

TEST_F(NygaDistributionTest, LeftConnectingPoint){
    ASSERT_EQ(induction_step.left_connecting_point_from_index(3), 3.5);
}

TEST_F(NygaDistributionTest, LeftConnectingPointEdgeCase){
    ASSERT_EQ(induction_step.left_connecting_point(), 1);
}

TEST_F(NygaDistributionTest, RightConnectingPoint){
    ASSERT_EQ(induction_step.right_connecting_point_from_index(5), 8);
}

TEST_F(NygaDistributionTest, RightConnectingPointEdgeCase){
    ASSERT_EQ(induction_step.right_connecting_point(), 9);
}

TEST_F(NygaDistributionTest, CreateUniformDistribution){
    auto uniform = induction_step.create_uniform_distribution_from_indices(3, 5);
    ASSERT_EQ(uniform->support->lower(), 3.5);
    ASSERT_EQ(uniform->support->upper(), 8);
}

TEST_F(NygaDistributionTest, CreateUniformDistributionEdgeCase){
    auto uniform = induction_step.create_uniform_distribution();
    ASSERT_EQ(uniform->support->lower(), 1);
    ASSERT_EQ(uniform->support->upper(), 9);
}

TEST_F(NygaDistributionTest, SumWeights){
    ASSERT_DOUBLE_EQ(induction_step.sum_log_weights(), 6 * log(1. / 6));
    ASSERT_DOUBLE_EQ(induction_step.sum_log_weights_from_indices(3, 5), 2 * log(1. / 6.));
}

TEST_F(NygaDistributionTest, ComputeBestSplit){
    auto [likelihood, split_index] = induction_step.compute_best_split();
    ASSERT_EQ(split_index, 1);
}

TEST_F(NygaDistributionTest, ComputeBestSplitWithoutResult){
    model->min_samples_per_quantile = 4;
    auto [likelihood, split_index] = induction_step.compute_best_split();
    ASSERT_EQ(split_index, -1);
    ASSERT_EQ(likelihood, -std::numeric_limits<double>::infinity());
}

TEST_F(NygaDistributionTest, ComputeBestSplitWithInducedIndices){
    model->min_samples_per_quantile = 3;
    auto [likelihood, split_index] = induction_step.compute_best_split();
    ASSERT_EQ(split_index, 3);
}

TEST_F(NygaDistributionTest, Fit){
    auto normal = std::normal_distribution<double>(0, 1);
    std::default_random_engine generator(69);
    auto data = std::make_shared<DataVector>(DataVector(1000));
    std::generate(data->begin(), data->end(), [&](){return normal(generator);});

    model->min_samples_per_quantile = 20;
    auto result = model->fit(data);
    ASSERT_GE(result->sub_circuits.size(), 1);
}

TEST_F(NygaDistributionTest, FitWithSingularData){
    auto data = std::make_shared<DataVector>(DataVector{1, 1, 1});
    auto result = model->fit(data);
    ASSERT_EQ(result->sub_circuits.size(), 1);
    auto subcircuit = std::static_pointer_cast<DiracDeltaDistribution>(result->sub_circuits[0]);
    ASSERT_EQ(subcircuit->location, 1);
    ASSERT_EQ(result->weights.size(), 1);
    ASSERT_EQ(result->weights[0], 1);
}