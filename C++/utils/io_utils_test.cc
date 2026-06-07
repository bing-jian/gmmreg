#include <cstdio>
#include <fstream>
#include <string>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "utils/io_utils.h"

using namespace gmmreg;

// RAII wrapper: removes the file on destruction so tests leave no temp files.
struct TempFile {
    const char* path;
    explicit TempFile(const char* p) : path(p) {}
    ~TempFile() { std::remove(path); }
};

// ── LoadMatrixFromTxt ─────────────────────────────────────────────────────────

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

TEST(LoadMatrixFromTxt, ScientificNotation) {
    const char* path = "/tmp/gmmreg_io_load_sci.txt";
    TempFile tf(path);
    { std::ofstream f(path); f << "2.816e-01 5.977e-01\n3.103e-01 6.149e-01\n"; }
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt(path, m), 2);
    EXPECT_NEAR(m(0, 0), 2.816e-01, 1e-6);
    EXPECT_NEAR(m(1, 1), 6.149e-01, 1e-6);
}

TEST(LoadMatrixFromTxt, MissingFileReturnsMinusOne) {
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt("/tmp/gmmreg_nonexistent_xyz.txt", m), -1);
}

TEST(LoadMatrixFromTxt, EmptyFileReturnsMinusOne) {
    const char* path = "/tmp/gmmreg_io_load_empty.txt";
    TempFile tf(path);
    { std::ofstream f(path); }
    vnl_matrix<double> m;
    EXPECT_EQ(LoadMatrixFromTxt(path, m), -1);
}

// ── SaveMatrixToAsciiFile ─────────────────────────────────────────────────────

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

TEST(SaveMatrixToAsciiFile, EmptyFilenameIsNoOp) {
    vnl_matrix<double> m(2, 2);
    m.fill(1.0);
    EXPECT_NO_THROW(SaveMatrixToAsciiFile("", m));
}

// ── SaveVectorToAsciiFile ─────────────────────────────────────────────────────

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

TEST(SaveVectorToAsciiFile, EmptyFilenameIsNoOp) {
    vnl_vector<double> v(3, 1.0);
    EXPECT_NO_THROW(SaveVectorToAsciiFile("", v));
}

// ── ReadAll ───────────────────────────────────────────────────────────────────

TEST(ReadAll, ReadsFileContents) {
    const char* path = "/tmp/gmmreg_io_readall.txt";
    TempFile tf(path);
    { std::ofstream f(path); f << "hello\nworld\n"; }
    EXPECT_EQ(ReadAll(path), "hello\nworld\n");
}

TEST(ReadAll, MissingFileReturnsEmpty) {
    EXPECT_TRUE(ReadAll("/tmp/gmmreg_nonexistent_readall.txt").empty());
}
