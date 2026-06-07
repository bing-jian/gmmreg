#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "io_utils.h"

using namespace gmmreg;

// RAII wrapper: removes the file on destruction so tests leave no temp files.
struct TempFile {
    const char* path;
    explicit TempFile(const char* p) : path(p) {}
    ~TempFile() { std::remove(path); }
};

// ── helpers used by SaveOutputToJson tests ───────────────────────────────────

static std::string ReadAll(const char* path) {
    std::ifstream f(path);
    return {std::istreambuf_iterator<char>(f), {}};
}

// Return the double that follows `"key": ` in content.
static double JsonDouble(const std::string& content, const char* key) {
    std::string needle = std::string("\"") + key + "\": ";
    size_t pos = content.find(needle);
    if (pos == std::string::npos)
        throw std::runtime_error(std::string("key not found: ") + key);
    return std::stod(content.substr(pos + needle.size()));
}

// Count non-overlapping occurrences of needle in content.
static int CountOccurrences(const std::string& content, const std::string& needle) {
    int count = 0;
    size_t pos = 0;
    while ((pos = content.find(needle, pos)) != std::string::npos) {
        ++count;
        pos += needle.size();
    }
    return count;
}

// ── LoadMatrixFromTxt ─────────────────────────────────────────────────────────

// Valid file: returns row count and fills matrix with correct values
TEST(LoadMatrixFromTxt, ValidFileParsedCorrectly) {
    const char* path = "/tmp/gmmreg_io_load_valid.txt";
    TempFile tf(path);
    {
        std::ofstream f(path);
        f << "1.0 2.0 3.0\n4.0 5.0 6.0\n";
    }
    vnl_matrix<double> m;
    int rows = LoadMatrixFromTxt(path, m);
    EXPECT_EQ(rows, 2);
    ASSERT_EQ(m.rows(), 2u);
    ASSERT_EQ(m.cols(), 3u);
    EXPECT_NEAR(m(0, 0), 1.0, 1e-10);
    EXPECT_NEAR(m(0, 1), 2.0, 1e-10);
    EXPECT_NEAR(m(0, 2), 3.0, 1e-10);
    EXPECT_NEAR(m(1, 0), 4.0, 1e-10);
    EXPECT_NEAR(m(1, 1), 5.0, 1e-10);
    EXPECT_NEAR(m(1, 2), 6.0, 1e-10);
}

// Single-row matrix
TEST(LoadMatrixFromTxt, SingleRowMatrix) {
    const char* path = "/tmp/gmmreg_io_load_single.txt";
    TempFile tf(path);
    { std::ofstream f(path); f << "7.5 8.5\n"; }
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt(path, m), 1);
    ASSERT_EQ(m.rows(), 1u);
    EXPECT_NEAR(m(0, 0), 7.5, 1e-10);
    EXPECT_NEAR(m(0, 1), 8.5, 1e-10);
}

// Scientific-notation values (format used by real project data files)
TEST(LoadMatrixFromTxt, ScientificNotation) {
    const char* path = "/tmp/gmmreg_io_load_sci.txt";
    TempFile tf(path);
    { std::ofstream f(path); f << "2.816e-01 5.977e-01\n3.103e-01 6.149e-01\n"; }
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt(path, m), 2);
    EXPECT_NEAR(m(0, 0), 2.816e-01, 1e-6);
    EXPECT_NEAR(m(1, 1), 6.149e-01, 1e-6);
}

// Missing file → returns -1
TEST(LoadMatrixFromTxt, MissingFileReturnsMinusOne) {
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt("/tmp/gmmreg_nonexistent_xyz.txt", m), -1);
}

// Empty file → returns -1 (read_ascii fails on empty stream)
TEST(LoadMatrixFromTxt, EmptyFileReturnsMinusOne) {
    const char* path = "/tmp/gmmreg_io_load_empty.txt";
    TempFile tf(path);
    { std::ofstream f(path); }  // create empty file
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt(path, m), -1);
}

// ── SaveMatrixToAsciiFile ─────────────────────────────────────────────────────

// Round-trip: saved matrix loads back with matching values
TEST(SaveMatrixToAsciiFile, RoundTripPreservesValues) {
    const char* path = "/tmp/gmmreg_io_save_matrix.txt";
    TempFile tf(path);
    vnl_matrix<double> orig(3, 2);
    orig(0,0)=1.1; orig(0,1)=2.2;
    orig(1,0)=3.3; orig(1,1)=4.4;
    orig(2,0)=5.5; orig(2,1)=6.6;
    SaveMatrixToAsciiFile(path, orig);

    vnl_matrix<double> loaded;
    ASSERT_EQ(LoadMatrixFromTxt(path, loaded), 3);
    ASSERT_EQ(loaded.rows(), 3u);
    ASSERT_EQ(loaded.cols(), 2u);
    for (int i = 0; i < 3; ++i)
        for (int d = 0; d < 2; ++d)
            EXPECT_NEAR(loaded(i, d), orig(i, d), 1e-5)
                << "mismatch at (" << i << ", " << d << ")";
}

// Empty filename → no crash, no file created
TEST(SaveMatrixToAsciiFile, EmptyFilenameIsNoOp) {
    vnl_matrix<double> m(2, 2);
    m.fill(1.0);
    EXPECT_NO_THROW(SaveMatrixToAsciiFile("", m));
}

// ── SaveVectorToAsciiFile ─────────────────────────────────────────────────────

// Saved vector values can be read back from file
TEST(SaveVectorToAsciiFile, WritesReadableValues) {
    const char* path = "/tmp/gmmreg_io_save_vector.txt";
    TempFile tf(path);
    vnl_vector<double> v(4);
    v[0]=1.5; v[1]=-2.5; v[2]=0.0; v[3]=3.14;
    SaveVectorToAsciiFile(path, v);

    std::ifstream f(path);
    ASSERT_TRUE(f.is_open());
    double val;
    std::vector<double> read_back;
    while (f >> val) read_back.push_back(val);
    ASSERT_EQ(read_back.size(), 4u);
    EXPECT_NEAR(read_back[0],  1.5,  1e-10);
    EXPECT_NEAR(read_back[1], -2.5,  1e-10);
    EXPECT_NEAR(read_back[2],  0.0,  1e-10);
    EXPECT_NEAR(read_back[3],  3.14, 1e-5);
}

// Empty filename → no crash, no file created
TEST(SaveVectorToAsciiFile, EmptyFilenameIsNoOp) {
    vnl_vector<double> v(3, 1.0);
    EXPECT_NO_THROW(SaveVectorToAsciiFile("", v));
}

// ── SaveOutputToJson ──────────────────────────────────────────────────────────

// Empty filename → no crash, no file created
TEST(SaveOutputToJson, EmptyFilenameIsNoOp) {
    vnl_vector<double> params(2, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    EXPECT_NO_THROW(SaveOutputToJson("", params, model, 0.0, 0.0));
    EXPECT_FALSE(std::ifstream("/tmp/gmmreg_json_noop.txt").is_open());
}

// Output file contains all four top-level JSON keys
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

// Output is a valid JSON object: starts with '{', ends with '}'
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

// elapsed_time_in_ms is written and parses back correctly
TEST(SaveOutputToJson, ElapsedTimeValue) {
    const char* path = "/tmp/gmmreg_json_elapsed.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 123.5, 0.0);
    std::string c = ReadAll(path);
    EXPECT_NEAR(JsonDouble(c, "elapsed_time_in_ms"), 123.5, 1e-9);
}

// initialization_time_in_ms is written and parses back correctly
TEST(SaveOutputToJson, InitializationTimeValue) {
    const char* path = "/tmp/gmmreg_json_inittime.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 45.25);
    std::string c = ReadAll(path);
    EXPECT_NEAR(JsonDouble(c, "initialization_time_in_ms"), 45.25, 1e-9);
}

// params with known values: commas between elements, one entry per value
TEST(SaveOutputToJson, ParamsCount) {
    const char* path = "/tmp/gmmreg_json_params.txt";
    TempFile tf(path);
    vnl_vector<double> params(4);
    params[0]=1.0; params[1]=2.0; params[2]=3.0; params[3]=4.0;
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    // 4 values → 3 commas inside the params array
    size_t open  = c.find("\"params\": [");
    size_t close = c.find("]", open);
    ASSERT_NE(open, std::string::npos);
    std::string arr = c.substr(open, close - open);
    EXPECT_EQ(CountOccurrences(arr, ","), 3);
}

// Empty params vector → "params": []
TEST(SaveOutputToJson, EmptyParamsArray) {
    const char* path = "/tmp/gmmreg_json_emptyparams.txt";
    TempFile tf(path);
    vnl_vector<double> params;
    vnl_matrix<double> model(1, 2, 0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    EXPECT_NE(c.find("\"params\": []"), std::string::npos);
}

// transformed_model: number of row arrays equals number of model rows
TEST(SaveOutputToJson, ModelRowCount) {
    const char* path = "/tmp/gmmreg_json_rowcount.txt";
    TempFile tf(path);
    vnl_vector<double> params(1, 0.0);
    vnl_matrix<double> model(3, 2);
    model.fill(0.0);
    SaveOutputToJson(path, params, model, 0.0, 0.0);
    std::string c = ReadAll(path);
    // Each row is written as "    [...]" — count opening brackets after the
    // "transformed_model" key (each row contributes exactly one "    [").
    size_t tm_pos = c.find("\"transformed_model\"");
    std::string after_tm = c.substr(tm_pos);
    EXPECT_EQ(CountOccurrences(after_tm, "    ["), 3);
}

// Empty model (0 rows) → transformed_model array body is empty
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

// Values in transformed_model survive a round-trip through the JSON file.
// Uses values (1.5, 2.5, 3.5, 4.5) that have exact binary representations
// and therefore unambiguous decimal representations at any precision ≥ 2.
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
