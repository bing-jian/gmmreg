#include <cstdio>

#include <gtest/gtest.h>

#include "gmmreg_api.h"
#include "utils/io_utils.h"
#include "utils/json_utils.h"

using namespace gmmreg;

class TpsKcFishTest : public ::testing::Test {
protected:
    void TearDown() override {
        std::remove("output.json");
        std::remove("transformed_model.txt");
        std::remove("final_affine.txt");
        std::remove("final_tps.txt");
    }
};

TEST_F(TpsKcFishTest, OutputMatchesExpected) {
    ASSERT_EQ(gmmreg_api("fish_full.ini", "tps_kc"), 0)
        << "gmmreg_api returned non-zero for tps_kc";

    std::string actual   = ReadAll("output.json");
    std::string expected = ReadAll("expected_output/fish_full/tps_kc.json");
    ASSERT_FALSE(actual.empty())   << "output.json was not created";
    ASSERT_FALSE(expected.empty()) << "expected_output/fish_full/tps_kc.json not found";

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

// ── fish_partial (2-D) ───────────────────────────────────────────────────────

class TpsKcFishPartialTest : public ::testing::Test {
protected:
    void TearDown() override {
        std::remove("output.json");
        std::remove("transformed_model.txt");
        std::remove("final_affine.txt");
        std::remove("final_tps.txt");
    }
};

TEST_F(TpsKcFishPartialTest, OutputMatchesExpected) {
    ASSERT_EQ(gmmreg_api("fish_partial.ini", "tps_kc"), 0)
        << "gmmreg_api returned non-zero for fish_partial tps_kc";

    std::string actual   = ReadAll("output.json");
    std::string expected = ReadAll("expected_output/fish_partial/tps_kc.json");
    ASSERT_FALSE(actual.empty())   << "output.json was not created";
    ASSERT_FALSE(expected.empty()) << "expected_output/fish_partial/tps_kc.json not found";

    auto ap = ParseFlatArray(actual,   "params");
    auto ep = ParseFlatArray(expected, "params");
    ASSERT_EQ(ap.size(), ep.size()) << "params array length mismatch";
    for (size_t i = 0; i < ep.size(); ++i)
        EXPECT_NEAR(ap[i], ep[i], 1e-2) << "params[" << i << "] mismatch";

    auto am = ParseMatrix(actual,   "transformed_model");
    auto em = ParseMatrix(expected, "transformed_model");
    ASSERT_EQ(am.size(), em.size()) << "transformed_model element count mismatch";
    for (size_t i = 0; i < em.size(); ++i)
        EXPECT_NEAR(am[i], em[i], 8.0e-4)
            << "transformed_model flat index " << i << " mismatch";
}

// ── face (3-D) ────────────────────────────────────────────────────────────────

class TpsKcFaceTest : public ::testing::Test {
protected:
    void TearDown() override {
        std::remove("output.json");
        std::remove("transformed_model.txt");
        std::remove("final_affine.txt");
        std::remove("final_tps.txt");
    }
};

TEST_F(TpsKcFaceTest, OutputMatchesExpected) {
    ASSERT_EQ(gmmreg_api("face.ini", "tps_kc"), 0)
        << "gmmreg_api returned non-zero for face tps_kc";

    std::string actual   = ReadAll("output.json");
    std::string expected = ReadAll("expected_output/face/tps_kc.json");
    ASSERT_FALSE(actual.empty())   << "output.json was not created";
    ASSERT_FALSE(expected.empty()) << "expected_output/face/tps_kc.json not found";

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
