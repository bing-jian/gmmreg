// This program can be used to evaluate/test Gauss transform implementation.

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>

#include <vnl/vnl_matrix.h>

#include "fgt_utils.h"
#include "nanoflann.hpp"
#include "utils/gauss_transform.h"

typedef std::chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) \
  std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define time_now() std::chrono::high_resolution_clock::now()

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " pts_1 pts_2 scale [gradient]" << endl;
    cerr << "  pts1, pts2, gradient(optional) -- ascii text files " << endl;
    cerr << "  scale -- numerical scalar value" << endl;
    cerr << "Example: " << argv[0] << " pts_1.txt pts_2.txt 1.0 gradient.txt"
         << endl;
    return -1;
  }
  vnl_matrix<double> A;
  ifstream file1(argv[1]);
  A.read_ascii(file1);
  int n1 = A.rows();
  int d1 = A.cols();
  cout << "Loaded point set A (" << n1 << " x " << d1 << ") from " << argv[1]
       << endl;

  vnl_matrix<double> B;
  ifstream file2(argv[2]);
  B.read_ascii(file2);
  int n2 = B.rows();
  int d2 = B.cols();
  cout << "Loaded point set B (" << n2 << " x " << d2 << ") from " << argv[2]
       << endl;

  vnl_matrix<double> gradient;
  gradient.set_size(n1, d1);
  gradient.fill(0);
  double scale = atof(argv[3]);
  double res = 0;

  {
    cout << "------------ " << endl;
    cout << "Evaluating Gauss Transform using a direct brute-force method ... "
         << endl;
    TimeVar t1 = time_now();
    res = gmmreg::GaussTransform(A, B, scale, gradient);
    cout << "  result: " << res << endl;
    cout << "  time elapsed: " << duration(time_now() - t1) << " microseconds."
         << endl;

    if (argc > 4) {
      ofstream outfile(argv[4], ios_base::out);
      gradient.print(outfile);
    }
  }

  {
    cout << "------------ " << endl;
    cout << "Building k-d tree ... " << endl;
    TimeVar t1 = time_now();
    std::unique_ptr<gmmreg::NanoflannTree<double>> tree;
    tree.reset(new gmmreg::NanoflannTree<double>(B));
    tree->tree.buildIndex();
    cout << "  time elapsed: " << duration(time_now() - t1) << " microseconds."
         << endl;

    cout << "------------ " << endl;
    cout << "Evaluating truncated Gauss Transform using k-d tree ... " << endl;
    TimeVar t2 = time_now();
    res = gmmreg::FastGaussTransform(*(tree.get()), A, scale, gradient);
    cout << "  result: " << res << endl;
    cout << "  time elapsed: " << duration(time_now() - t2) << " microseconds."
         << endl;

    if (argc > 5) {
      ofstream outfile(argv[5], ios_base::out);
      gradient.print(outfile);
    }
  }

  return 0;
}
