#include <cstdio>

#include <gtest/gtest.h>

#include "gmmreg_api.h"
#include "utils/io_utils.h"
#include "utils/json_utils.h"

using namespace gmmreg;

class GrbfL2FishTest : public ::testing::Test {
protected:
    void TearDown() override {
        std::remove("output.json");
        std::remove("transformed_model.txt");
    }
};

TEST_F(GrbfL2FishTest, OutputMatchesExpected) {
    ASSERT_EQ(gmmreg_api("fish_full.ini", "grbf_l2"), 0)
        << "gmmreg_api returned non-zero for grbf_l2";

    std::string actual   = ReadAll("output.json");
    std::string expected = ReadAll("expected_output/fish_full/grbf_l2.json");
    ASSERT_FALSE(actual.empty())   << "output.json was not created";
    ASSERT_FALSE(expected.empty()) << "expected_output/fish_full/grbf_l2.json not found";

    auto ap = ParseFlatArray(actual,   "params");
    auto ep = ParseFlatArray(expected, "params");
    ASSERT_EQ(ap.size(), ep.size()) << "params array length mismatch";
    for (size_t i = 0; i < ep.size(); ++i)
        EXPECT_NEAR(ap[i], ep[i], 1e-2) << "params[" << i << "] mismatch";

    auto am = ParseMatrix(actual, "transformed_model");
    auto em = ParseMatrix(expected, "transformed_model");
    ASSERT_EQ(am.size(), em.size()) << "transformed_model element count mismatch";
    for (size_t i = 0; i < em.size(); ++i)
        EXPECT_NEAR(am[i], em[i], 8.0e-4)
            << "transformed_model flat index " << i << " mismatch";
}

// ── face (3-D) ────────────────────────────────────────────────────────────────

class GrbfL2FaceTest : public ::testing::Test {
protected:
    void TearDown() override {
        std::remove("output.json");
        std::remove("transformed_model.txt");
    }
};

TEST_F(GrbfL2FaceTest, OutputMatchesExpected) {
    ASSERT_EQ(gmmreg_api("face.ini", "grbf_l2"), 0)
        << "gmmreg_api returned non-zero for face grbf_l2";

    std::string actual   = ReadAll("output.json");
    std::string expected = ReadAll("expected_output/face/grbf_l2.json");
    ASSERT_FALSE(actual.empty())   << "output.json was not created";
    ASSERT_FALSE(expected.empty()) << "expected_output/face/grbf_l2.json not found";

    auto ap = ParseFlatArray(actual,   "params");
    auto ep = ParseFlatArray(expected, "params");
    ASSERT_EQ(ap.size(), ep.size()) << "params array length mismatch";
    for (size_t i = 0; i < ep.size(); ++i)
        EXPECT_NEAR(ap[i], ep[i], 1e-2) << "params[" << i << "] mismatch";

    auto am = ParseMatrix(actual,   "transformed_model");
    auto em = ParseMatrix(expected, "transformed_model");
    ASSERT_EQ(am.size(), em.size()) << "transformed_model element count mismatch";
    for (size_t i = 0; i < em.size(); ++i)
        EXPECT_NEAR(am[i], em[i], 5.0e-3)
            << "transformed_model flat index " << i << " mismatch";
}
