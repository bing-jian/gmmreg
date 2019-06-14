#include <fstream>
#include <iostream>
#include <memory>

#include <vnl/vnl_matrix.h>

#include "fgt_utils.h"
#include "utils/downsample_utils.h"

using std::cerr;
using std::cout;
using std::endl;


int main(int argc, char* argv[]) {
  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " InputPtsFile MinDistance OutputPtsFile" << endl;
    return -1;
  }
  std::ifstream infile(argv[1]);
  vnl_matrix<double> input_pts;
  input_pts.read_ascii(infile);

  int n = input_pts.rows();
  int d = input_pts.cols();

  cout << n << " "<< d << "-d points loaded from " << argv[1] << endl;

  std::unique_ptr<gmmreg::NanoflannTree<double>> tree;
  tree.reset(new gmmreg::NanoflannTree<double>(input_pts));
  tree->tree.buildIndex();
  vnl_matrix<double> kept = gmmreg::CreateSparseNodes(*(tree.get()), atof(argv[2]));
  std::ofstream outfile_kept(argv[3], std::ios_base::out);
  kept.print(outfile_kept);
  cout << "Input downsampled to " << kept.rows() << " points and saved to "
       << argv[3] << endl;
  return 0;
}
