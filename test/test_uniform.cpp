#include "gtest/gtest.h"
#include "probabilistic_model.h"
#include "interval.h"
#include "univariate.h"
#include "variable.h"

auto x = Continuous("x");

TEST(AtomicIntervalCreationTestSuite, SimpleInterval) {
SimpleInterval interval = SimpleInterval();
interval.lower = 0.0;
interval.upper = 1.0;
interval.left = BorderType::OPEN;
interval.right = BorderType::CLOSED;
EXPECT_EQ(interval.lower, 0.0);
EXPECT_EQ(interval.upper, 1.0);
EXPECT_EQ(interval.left, BorderType::OPEN);
EXPECT_EQ(interval.right, BorderType::CLOSED);
}

TEST(UniformDistribution, Likelihood) {

    auto uniform = UniformDistribution();

}