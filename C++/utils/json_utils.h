#ifndef GMMREG_UTILS_JSON_UTILS_H_
#define GMMREG_UTILS_JSON_UTILS_H_

#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace gmmreg {

// Write registration output as a structured JSON file.
// No-op when filename is empty.
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

// Parse a flat JSON array  "key": [v, v, …]  into a vector of doubles.
// Matches the format produced by SaveOutputToJson.
std::vector<double> ParseFlatArray(const std::string& json, const char* key);

// Parse a JSON array-of-arrays identified by key into a flat vector of doubles
// (row-major order). Matches the format produced by SaveOutputToJson.
std::vector<double> ParseMatrix(const std::string& json, const char* key);

}  // namespace gmmreg

#endif  // GMMREG_UTILS_JSON_UTILS_H_
