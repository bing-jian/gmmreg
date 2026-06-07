#include "json_utils.h"

#include <sstream>

namespace gmmreg {

std::vector<double> ParseFlatArray(const std::string& json, const char* key) {
  std::string needle = std::string("\"") + key + "\"";
  size_t key_pos = json.find(needle);
  if (key_pos == std::string::npos) return {};
  size_t open = json.find('[', key_pos + needle.size());
  if (open == std::string::npos) return {};
  size_t close = json.find(']', open + 1);
  if (close == std::string::npos) return {};
  std::string seg = json.substr(open + 1, close - open - 1);
  for (char& c : seg) if (c == ',') c = ' ';
  std::istringstream ss(seg);
  std::vector<double> vals;
  double v;
  while (ss >> v) vals.push_back(v);
  return vals;
}

std::vector<double> ParseMatrix(const std::string& json, const char* key) {
  std::string needle = std::string("\"") + key + "\"";
  size_t key_pos = json.find(needle);
  if (key_pos == std::string::npos) return {};
  size_t outer_open = json.find('[', key_pos + needle.size());
  if (outer_open == std::string::npos) return {};
  // Walk forward counting brackets to find the matching outer close.
  size_t pos = outer_open + 1;
  int depth = 1;
  while (pos < json.size() && depth > 0) {
    if (json[pos] == '[') ++depth;
    else if (json[pos] == ']') --depth;
    ++pos;
  }
  if (depth != 0) return {};
  // pos now points one past the outer ']'
  std::string seg = json.substr(outer_open + 1, pos - outer_open - 2);
  for (char& c : seg) if (c == '[' || c == ']' || c == ',') c = ' ';
  std::istringstream ss(seg);
  std::vector<double> vals;
  double v;
  while (ss >> v) vals.push_back(v);
  return vals;
}

}  // namespace gmmreg
