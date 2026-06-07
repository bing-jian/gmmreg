#include <iostream>
#include <string>
#include <vector>

#include "gmmreg_api.h"

int main(int argc, char* argv[]) {
  bool verbose = false;
  std::vector<const char*> args;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--verbose") {
      verbose = true;
    } else {
      args.push_back(argv[i]);
    }
  }

  if (args.size() < 2) {
    std::cerr << "Usage: " << argv[0]
              << " [--verbose] config_file method" << std::endl;
    print_usage();
    return -1;
  }

  if (verbose) {
    std::cout << "Robust Point Set Registration Using Gaussian Mixture Models"
              << std::endl;
    std::cout << "Compiled on " << __DATE__ << ", " << __TIME__ << std::endl;
    std::cout << "Copyright Bing Jian & Baba C. Vemuri" << std::endl;
  }

  gmmreg_api(args[0], args[1]);
  return 0;
}
