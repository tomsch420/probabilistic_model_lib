#include "gtest/gtest.h"
#include "probabilistic_circuit.h"
#include "interval.h"
#include "univariate.h"
#include "variable.h"


class SmoothSumUnitTest : public testing::Test {
public:
    ContinuousPtr_t variable_x;
    SymbolicPtr_t variable_a;
    IntegerPtr_t variable_i;
    SmoothSumUnit model;

    SmoothSumUnitTest() {
        variable_x = make_shared_continuous("x");
        variable_a = make_shared_symbolic(std::make_shared<std::string>("symbolic_a"),
                                 make_shared_all_elements(std::set<std::string>{"symbolic_a", "b", "c"}));
        variable_i = make_shared_integer("i");

        model = SmoothSumUnit();
        auto u1 = UniformDistribution::make_shared(variable_x, closed_open<double>(0, 2));
        auto u2 = UniformDistribution::make_shared(variable_x, closed_open<double>(3, 8));
        model.add_subcircuit(0.5, u1);
        model.add_subcircuit(0.5, u2);
    }
};

TEST_F(SmoothSumUnitTest, Likelihood) {
    auto event = std::make_shared<FullEvidence>(FullEvidence{1.0});
    EXPECT_DOUBLE_EQ(model.likelihood(event), 0.5 * 0.5);
    EXPECT_DOUBLE_EQ(model.log_likelihood(event), log(0.5 * 0.5));

    auto event2 = std::make_shared<FullEvidence>(FullEvidence{4.0});
    EXPECT_DOUBLE_EQ(model.likelihood(event2), 0.5 * 0.2);
    EXPECT_DOUBLE_EQ(model.log_likelihood(event2), log(0.5 * 0.2));

    auto event3 = std::make_shared<FullEvidence>(FullEvidence{2.5});
    EXPECT_DOUBLE_EQ(model.likelihood(event3),0);
    EXPECT_DOUBLE_EQ(model.log_likelihood(event3), log(0));
}

TEST_F(SmoothSumUnitTest, Variables) {
    auto variables = model.get_variables();
    EXPECT_EQ(variables->size(), 1);
    EXPECT_TRUE(variables->find(variable_x) != variables->end());
}

class DecomposableProductUnitTest : public testing::Test {
public:
    ContinuousPtr_t variable_x;
    ContinuousPtr_t variable_y;
    DecomposableProductUnit model;

    DecomposableProductUnitTest() {
        variable_x = make_shared_continuous("x");
        variable_y = make_shared_continuous("y");

        model = DecomposableProductUnit();
        auto u1 = UniformDistribution::make_shared(variable_x, closed_open<double>(0, 2));
        auto u2 = UniformDistribution::make_shared(variable_y, closed_open<double>(0, 1));
        model.add_subcircuit(u1);
        model.add_subcircuit(u2);
    }
};

TEST_F(DecomposableProductUnitTest, Variables) {
    auto variables = model.get_variables();
    EXPECT_EQ(variables->size(), 2);
    EXPECT_TRUE(variables->find(variable_x) != variables->end());
    EXPECT_TRUE(variables->find(variable_y) != variables->end());
}

TEST_F(DecomposableProductUnitTest, Likelihood) {
    auto event = std::make_shared<FullEvidence>(FullEvidence{1.0, 0.5});
    EXPECT_DOUBLE_EQ(model.likelihood(event), 0.5);
    EXPECT_DOUBLE_EQ(model.log_likelihood(event), log(0.5));

    auto event2 = std::make_shared<FullEvidence>(FullEvidence{2.0, 0.5});
    EXPECT_DOUBLE_EQ(model.likelihood(event2), 0);
    EXPECT_DOUBLE_EQ(model.log_likelihood(event2), log(0));
}
