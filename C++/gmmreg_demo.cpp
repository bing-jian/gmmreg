/*=========================================================================
$Author: bingjian $
$Date: 2013-01-04 01:39:25 -0500 (Fri, 04 Jan 2013) $
$Revision: 145 $
=========================================================================*/

#include <iostream>

#include "gmmreg_api.h"

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " config_file method" << std::endl;
    print_usage();
    return -1;
  }
  gmmreg_api(argv[1], argv[2]);
  return 0;
}
