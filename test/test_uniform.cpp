#include "gtest/gtest.h"
#include "probabilistic_model.h"
#include "interval.h"
#include "univariate.h"
#include "variable.h"

auto continuous_x = make_shared_continuous("continuous_x");
auto symbolic_a = make_shared_symbolic(std::make_shared<std::string>("symbolic_a"),
                                       make_shared_all_elements(std::set<std::string>{"symbolic_a", "b", "c"}));
auto integer_i = make_shared_integer("integer_i");


TEST(UniformDistribution, Likelihood) {
    auto d1 = UniformDistribution(continuous_x, closed_open<double>(0, 2));
    EXPECT_EQ(d1.log_pdf(0.5), log(0.5));
    EXPECT_EQ(d1.log_pdf(2), log(0.));
    EXPECT_EQ(d1.log_pdf(3), log(0.));
}

TEST(SymbolicDistribution, Likelihood){
    auto d2 = SymbolicDistribution(symbolic_a, std::map<int, double>{{0, 0.7}, {2, 0.3}});
    EXPECT_EQ(d2.log_likelihood(std::make_shared<FullEvidence>(std::vector<double>{0})), log(0.7));
    EXPECT_EQ(d2.log_likelihood(std::make_shared<FullEvidence>(std::vector<double>{1})), log(0.));
    EXPECT_EQ(d2.log_likelihood(std::make_shared<FullEvidence>(std::vector<double>{2})), log(0.3));
}

TEST(SymbolicDistribution, Support){
    auto d2 = SymbolicDistribution(symbolic_a, std::map<int, double>{{0, 0.7}, {2, 0.3}});
    auto support = std::static_pointer_cast<Set>(d2.get_support());
    EXPECT_EQ(support->all_elements->size(), 3);
    EXPECT_EQ(support->simple_sets->size(), 2);
}

TEST(IntegerDistribution, Support){
    auto d3 = IntegerDistribution(integer_i, std::map<int, double>{{0, 0.7}, {2, 0.3}});
    auto support = std::static_pointer_cast<Interval<int>>(d3.get_support());
    EXPECT_EQ(support->lower(), 0);
    EXPECT_EQ(support->upper(), 2);
    EXPECT_EQ(support->simple_sets->size(), 2);
}

TEST(DiracDeltaDistribution, Likelihood){
    auto d4 = DiracDeltaDistribution(continuous_x, 1, 3);
    EXPECT_EQ(d4.pdf(1), 3);
    EXPECT_EQ(d4.pdf(2), 0);
}