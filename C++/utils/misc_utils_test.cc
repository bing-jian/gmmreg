#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "misc_utils.h"

using namespace gmmreg::utils;

// ── parse_tokens (float overload) ────────────────────────────────────────────

TEST(ParseTokensFloat, SpaceDelimited) {
    char buf[] = "1.5 2.5 3.5";
    std::vector<float> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_FLOAT_EQ(v[0], 1.5f);
    EXPECT_FLOAT_EQ(v[1], 2.5f);
    EXPECT_FLOAT_EQ(v[2], 3.5f);
}

TEST(ParseTokensFloat, CommaDelimited) {
    char buf[] = "0.1,0.2,0.3";
    std::vector<float> v;
    parse_tokens(buf, ",", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_FLOAT_EQ(v[0], 0.1f);
    EXPECT_FLOAT_EQ(v[1], 0.2f);
    EXPECT_FLOAT_EQ(v[2], 0.3f);
}

// Multiple delimiter characters: any char in the set acts as a separator
TEST(ParseTokensFloat, MultipleDelimChars) {
    char buf[] = "1.0 2.0,3.0";
    std::vector<float> v;
    parse_tokens(buf, " ,", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_FLOAT_EQ(v[0], 1.0f);
    EXPECT_FLOAT_EQ(v[1], 2.0f);
    EXPECT_FLOAT_EQ(v[2], 3.0f);
}

// strtok skips leading and consecutive delimiters
TEST(ParseTokensFloat, LeadingAndConsecutiveDelimiters) {
    char buf[] = "  1.0  2.0  ";
    std::vector<float> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 2u);
    EXPECT_FLOAT_EQ(v[0], 1.0f);
    EXPECT_FLOAT_EQ(v[1], 2.0f);
}

TEST(ParseTokensFloat, SingleToken) {
    char buf[] = "3.14";
    std::vector<float> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 1u);
    EXPECT_FLOAT_EQ(v[0], 3.14f);
}

TEST(ParseTokensFloat, NegativeValues) {
    char buf[] = "-1.5 0.0 2.5";
    std::vector<float> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_FLOAT_EQ(v[0], -1.5f);
    EXPECT_FLOAT_EQ(v[1],  0.0f);
    EXPECT_FLOAT_EQ(v[2],  2.5f);
}

// Empty string produces no tokens
TEST(ParseTokensFloat, EmptyString) {
    char buf[] = "";
    std::vector<float> v;
    parse_tokens(buf, " ", v);
    EXPECT_TRUE(v.empty());
}

// Tokens are appended to an already-populated vector
TEST(ParseTokensFloat, AppendsToExistingVector) {
    char buf[] = "3.0 4.0";
    std::vector<float> v = {1.0f, 2.0f};
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 4u);
    EXPECT_FLOAT_EQ(v[0], 1.0f);
    EXPECT_FLOAT_EQ(v[1], 2.0f);
    EXPECT_FLOAT_EQ(v[2], 3.0f);
    EXPECT_FLOAT_EQ(v[3], 4.0f);
}

// ── parse_tokens (int overload) ───────────────────────────────────────────────

TEST(ParseTokensInt, SpaceDelimited) {
    char buf[] = "10 20 30";
    std::vector<int> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_EQ(v[0], 10);
    EXPECT_EQ(v[1], 20);
    EXPECT_EQ(v[2], 30);
}

TEST(ParseTokensInt, CommaDelimited) {
    char buf[] = "1,2,3";
    std::vector<int> v;
    parse_tokens(buf, ",", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
}

TEST(ParseTokensInt, NegativeValues) {
    char buf[] = "-3 0 7";
    std::vector<int> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 3u);
    EXPECT_EQ(v[0], -3);
    EXPECT_EQ(v[1],  0);
    EXPECT_EQ(v[2],  7);
}

TEST(ParseTokensInt, SingleToken) {
    char buf[] = "42";
    std::vector<int> v;
    parse_tokens(buf, " ", v);
    ASSERT_EQ(v.size(), 1u);
    EXPECT_EQ(v[0], 42);
}

TEST(ParseTokensInt, EmptyString) {
    char buf[] = "";
    std::vector<int> v;
    parse_tokens(buf, " ", v);
    EXPECT_TRUE(v.empty());
}
