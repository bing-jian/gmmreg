#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gmmreg_api.h"

// ── minimal JSON helpers ─────────────────────────────────────────────────────

static std::string ReadAll(const char* path) {
    std::ifstream f(path);
    return {std::istreambuf_iterator<char>(f), {}};
}

// Parse a flat JSON array "key": [v, v, …] into a vector of doubles.
static std::vector<double> ParseFlatArray(const std::string& json,
                                          const char* key) {
    std::string needle = std::string("\"") + key + "\": [";
    size_t start = json.find(needle);
    if (start == std::string::npos) return {};
    start += needle.size();
    size_t end = json.find(']', start);
    if (end == std::string::npos) return {};
    std::string seg = json.substr(start, end - start);
    for (char& c : seg) if (c == ',') c = ' ';
    std::istringstream ss(seg);
    std::vector<double> vals;
    double v;
    while (ss >> v) vals.push_back(v);
    return vals;
}

// Flatten all numbers from the "transformed_model" array-of-arrays section.
static std::vector<double> ParseTransformedModel(const std::string& json) {
    const char* marker = "\"transformed_model\": [\n";
    size_t start = json.find(marker);
    if (start == std::string::npos) return {};
    start += strlen(marker);
    size_t end = json.find("\n  ]", start);
    if (end == std::string::npos) return {};
    std::string seg = json.substr(start, end - start);
    for (char& c : seg) if (c == '[' || c == ']' || c == ',') c = ' ';
    std::istringstream ss(seg);
    std::vector<double> vals;
    double v;
    while (ss >> v) vals.push_back(v);
    return vals;
}

// ── integration test ─────────────────────────────────────────────────────────

// This test is run with WORKING_DIRECTORY set to C++/testdata so that all
// relative paths in fish_full.ini (./fish_data/…, ./output.json, …) resolve
// against the testdata directory.
class TpsL2FishTest : public ::testing::Test {
protected:
    void TearDown() override {
        std::remove("output.json");
        std::remove("transformed_model.txt");
        std::remove("final_affine.txt");
        std::remove("final_tps.txt");
    }
};

TEST_F(TpsL2FishTest, OutputMatchesExpected) {
    ASSERT_EQ(gmmreg_api("fish_full.ini", "tps_l2"), 0)
        << "gmmreg_api returned non-zero for tps_l2";

    std::string actual   = ReadAll("output.json");
    std::string expected = ReadAll("expected_output/fish_full/tps_l2.json");
    ASSERT_FALSE(actual.empty())   << "output.json was not created";
    ASSERT_FALSE(expected.empty()) << "expected_output/fish_full/tps_l2.json not found";

    // ── params ───────────────────────────────────────────────────────────────
    auto ap = ParseFlatArray(actual,   "params");
    auto ep = ParseFlatArray(expected, "params");
    ASSERT_EQ(ap.size(), ep.size()) << "params array length mismatch";
    for (size_t i = 0; i < ep.size(); ++i)
        EXPECT_NEAR(ap[i], ep[i], 1e-2) << "params[" << i << "] mismatch";

    // ── transformed_model ────────────────────────────────────────────────────
    auto am = ParseTransformedModel(actual);
    auto em = ParseTransformedModel(expected);
    ASSERT_EQ(am.size(), em.size()) << "transformed_model element count mismatch";
    for (size_t i = 0; i < em.size(); ++i)
        EXPECT_NEAR(am[i], em[i], 8.0 * 1e-4)
            << "transformed_model flat index " << i << " mismatch";
}

// ── tps_kc ───────────────────────────────────────────────────────────────────

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

    // ── params ───────────────────────────────────────────────────────────────
    auto ap = ParseFlatArray(actual,   "params");
    auto ep = ParseFlatArray(expected, "params");
    ASSERT_EQ(ap.size(), ep.size()) << "params array length mismatch";
    for (size_t i = 0; i < ep.size(); ++i)
        EXPECT_NEAR(ap[i], ep[i], 1e-2) << "params[" << i << "] mismatch";

    // ── transformed_model ────────────────────────────────────────────────────
    auto am = ParseTransformedModel(actual);
    auto em = ParseTransformedModel(expected);
    ASSERT_EQ(am.size(), em.size()) << "transformed_model element count mismatch";
    for (size_t i = 0; i < em.size(); ++i)
        EXPECT_NEAR(am[i], em[i], 8.0 * 1e-4)
            << "transformed_model flat index " << i << " mismatch";
}
