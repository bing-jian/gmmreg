#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "utils/io_utils.h"
#include "utils/json_utils.h"

using namespace gmmreg;

struct TempFile {
    const char* path;
    explicit TempFile(const char* p) : path(p) {}
    ~TempFile() { std::remove(path); }
};

// Return the double that follows `"key": ` in content.
static double JsonDouble(const std::string& content, const char* key) {
    std::string needle = std::string("\"") + key + "\": ";
    size_t pos = content.find(needle);
    if (pos == std::string::npos)
        throw std::runtime_error(std::string("key not found: ") + key);
    return std::stod(content.substr(pos + needle.size()));
}

// Count non-overlapping occurrences of needle in content.
static int CountOccurrences(const std::string& content,
                             const std::string& needle) {
    int count = 0;
    size_t pos = 0;
    while ((pos = content.find(needle, pos)) != std::string::npos) {
        ++count;
        pos += needle.size();
    }
    return count;
}

// ── SaveOutputToJson ──────────────────────────────────────────────────────────

TEST(SaveOutputToJson, EmptyFilenameIsNoOp) {
    vnl_vector<double> params(2, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    EXPECT_NO_THROW(SaveOutputToJson("", params, model, 0.0, 0.0));
}

TEST(SaveOutputToJson, AllKeysPresent) {
    const char* path = "/tmp/gmmreg_json_keys.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 1.0, 2.0);
    std::string c = ReadAll(path);
    EXPECT_NE(c.find("\"elapsed_time_in_ms\""),        std::string::npos);
    EXPECT_NE(c.find("\"initialization_time_in_ms\""), std::string::npos);
    EXPECT_NE(c.find("\"params\""),                    std::string::npos);
    EXPECT_NE(c.find("\"transformed_model\""),         std::string::npos);
}

TEST(SaveOutputToJson, OuterBraces) {
    const char* path = "/tmp/gmmreg_json_braces.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    EXPECT_EQ(c.front(), '{');
    EXPECT_NE(c.rfind('}'), std::string::npos);
}

TEST(SaveOutputToJson, ElapsedTimeValue) {
    const char* path = "/tmp/gmmreg_json_elapsed.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 123.5, 0.0);
    std::string c = ReadAll(path);
    EXPECT_NEAR(JsonDouble(c, "elapsed_time_in_ms"), 123.5, 1e-9);
}

TEST(SaveOutputToJson, InitializationTimeValue) {
    const char* path = "/tmp/gmmreg_json_inittime.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 45.25);
    std::string c = ReadAll(path);
    EXPECT_NEAR(JsonDouble(c, "initialization_time_in_ms"), 45.25, 1e-9);
}

TEST(SaveOutputToJson, ParamsCount) {
    const char* path = "/tmp/gmmreg_json_params.txt";
    TempFile tf(path);
    vnl_vector<double> params(4);
    params[0]=1.0; params[1]=2.0; params[2]=3.0; params[3]=4.0;
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    size_t open  = c.find("\"params\": [");
    size_t close = c.find("]", open);
    ASSERT_NE(open, std::string::npos);
    std::string arr = c.substr(open, close - open);
    EXPECT_EQ(CountOccurrences(arr, ","), 3);
}

TEST(SaveOutputToJson, EmptyParamsArray) {
    const char* path = "/tmp/gmmreg_json_emptyparams.txt";
    TempFile tf(path);
    vnl_vector<double> params;
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    EXPECT_NE(c.find("\"params\": []"), std::string::npos);
}

TEST(SaveOutputToJson, ModelRowCount) {
    const char* path = "/tmp/gmmreg_json_rowcount.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(3, 2);
    model.fill(0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    size_t tm_pos = c.find("\"transformed_model\"");
    std::string after_tm = c.substr(tm_pos);
    EXPECT_EQ(CountOccurrences(after_tm, "    ["), 3);
}

TEST(SaveOutputToJson, EmptyModelArray) {
    const char* path = "/tmp/gmmreg_json_emptymodel.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(0, 2);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    size_t tm_pos = c.find("\"transformed_model\"");
    std::string after_tm = c.substr(tm_pos);
    EXPECT_EQ(CountOccurrences(after_tm, "    ["), 0);
}

TEST(SaveOutputToJson, ModelValuesRoundTrip) {
    const char* path = "/tmp/gmmreg_json_modelvals.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(2, 2);
    model(0,0)=1.5; model(0,1)=2.5;
    model(1,0)=3.5; model(1,1)=4.5;
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    EXPECT_NE(c.find("1.5"), std::string::npos);
    EXPECT_NE(c.find("2.5"), std::string::npos);
    EXPECT_NE(c.find("3.5"), std::string::npos);
    EXPECT_NE(c.find("4.5"), std::string::npos);
}

// ── ParseFlatArray ────────────────────────────────────────────────────────────

TEST(ParseFlatArray, ParsesKnownValues) {
    std::string json = "{\n  \"params\": [1.0, 2.5, -3.0]\n}\n";
    auto v = ParseFlatArray(json, "params");
    ASSERT_EQ(v.size(), 3u);
    EXPECT_NEAR(v[0],  1.0, 1e-12);
    EXPECT_NEAR(v[1],  2.5, 1e-12);
    EXPECT_NEAR(v[2], -3.0, 1e-12);
}

TEST(ParseFlatArray, EmptyArray) {
    std::string json = "{\n  \"params\": []\n}\n";
    EXPECT_TRUE(ParseFlatArray(json, "params").empty());
}

TEST(ParseFlatArray, MissingKeyReturnsEmpty) {
    std::string json = "{\n  \"other\": [1.0]\n}\n";
    EXPECT_TRUE(ParseFlatArray(json, "params").empty());
}

TEST(ParseFlatArray, RoundTripWithSaveOutputToJson) {
    const char* path = "/tmp/gmmreg_json_flatarray_rt.txt";
    TempFile tf(path);
    vnl_vector<double> params(3);
    params[0] = 1.5; params[1] = -2.5; params[2] = 3.0;
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    auto v = ParseFlatArray(ReadAll(path), "params");
    ASSERT_EQ(v.size(), 3u);
    EXPECT_NEAR(v[0],  1.5, 1e-9);
    EXPECT_NEAR(v[1], -2.5, 1e-9);
    EXPECT_NEAR(v[2],  3.0, 1e-9);
}

// ── ParseMatrix ─────────────────────────────────────────────────────

TEST(ParseMatrix, RoundTripWithSaveOutputToJson) {
    const char* path = "/tmp/gmmreg_json_model_rt.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(2, 2);
    model(0,0) = 1.5; model(0,1) = 2.5;
    model(1,0) = 3.5; model(1,1) = 4.5;
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    auto v = ParseMatrix(ReadAll(path), "transformed_model");
    ASSERT_EQ(v.size(), 4u);
    EXPECT_NEAR(v[0], 1.5, 1e-9);
    EXPECT_NEAR(v[1], 2.5, 1e-9);
    EXPECT_NEAR(v[2], 3.5, 1e-9);
    EXPECT_NEAR(v[3], 4.5, 1e-9);
}

TEST(ParseMatrix, EmptyModel) {
    const char* path = "/tmp/gmmreg_json_emptymodel_rt.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(0, 2);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    EXPECT_TRUE(ParseMatrix(ReadAll(path), "transformed_model").empty());
}
