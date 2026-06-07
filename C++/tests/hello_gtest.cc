#include <gtest/gtest.h>

TEST(HelloGTest, SanityCheck) {
    EXPECT_EQ(1 + 1, 2);
    EXPECT_TRUE(true);
    EXPECT_NE(1, 2);
}
