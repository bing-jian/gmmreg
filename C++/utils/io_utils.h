#ifndef GMMREG_UTILS_IO_UTILS_H_
#define GMMREG_UTILS_IO_UTILS_H_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

namespace gmmreg {

template<typename T>
int LoadMatrixFromTxt(const char* filename, vnl_matrix<T>& matrix);

template<typename T>
void SaveMatrixToAsciiFile(const char * filename, const vnl_matrix<T>& x);

template<typename T>
void SaveVectorToAsciiFile(const char * filename, const vnl_vector<T>& x);


template<typename T>
int LoadMatrixFromTxt(const char* filename, vnl_matrix<T>& matrix) {
  std::ifstream infile(filename, std::ios_base::in);
  if (infile.is_open()) {
    if (matrix.read_ascii(infile)) {
      return matrix.rows();
    } else {
      std::cerr << "unable to parse input file " << filename
                << " as a matrix." << std::endl;
      return -1;
    }
  } else {
    std::cerr << "unable to open input file " << filename << std::endl;
    return -1;
  }
}

template<typename T>
void SaveMatrixToAsciiFile(const char * filename, const vnl_matrix<T>& x) {
  if (strlen(filename) > 0) {
    std::ofstream outfile(filename, std::ios_base::out);
    x.print(outfile);
  }
}

template<typename T>
void SaveVectorToAsciiFile(const char * filename, const vnl_vector<T>& x) {
  if (strlen(filename) > 0) {
    std::ofstream outfile(filename, std::ios_base::out);
    outfile << x;
  }
}

template<typename T>
void SaveOutputToJson(const char* filename,
                      const vnl_vector<T>& params,
                      const vnl_matrix<T>& transformed_model,
                      T elapsed_time_in_ms,
                      T initialization_time_in_ms) {
  if (strlen(filename) == 0) return;
  std::ofstream f(filename);
  if (!f.is_open()) return;
  f << std::setprecision(std::numeric_limits<T>::digits10 + 1);
  f << "{\n";
  f << "  \"elapsed_time_in_ms\": " << elapsed_time_in_ms << ",\n";
  f << "  \"initialization_time_in_ms\": " << initialization_time_in_ms << ",\n";
  f << "  \"params\": [";
  for (unsigned int i = 0; i < params.size(); ++i) {
    if (i) f << ", ";
    f << params[i];
  }
  f << "],\n";
  f << "  \"transformed_model\": [\n";
  for (unsigned int i = 0; i < transformed_model.rows(); ++i) {
    f << "    [";
    for (unsigned int d = 0; d < transformed_model.cols(); ++d) {
      if (d) f << ", ";
      f << transformed_model(i, d);
    }
    f << "]";
    if (i + 1 < transformed_model.rows()) f << ",";
    f << "\n";
  }
  f << "  ]\n";
  f << "}\n";
}

}  // namespace gmmreg

#endif // GMMREG_UTILS_IO_UTILS_H_
