#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Read entire file into a string.
inline std::string ReadAll(const char* path) {
    std::ifstream f(path);
    return {std::istreambuf_iterator<char>(f), {}};
}

// Parse a flat JSON array "key": [v, v, …] into a vector of doubles.
inline std::vector<double> ParseFlatArray(const std::string& json,
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
inline std::vector<double> ParseTransformedModel(const std::string& json) {
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
