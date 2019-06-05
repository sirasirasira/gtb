#include <iostream>
#include <sstream>
#include <fstream>
#include "getopt.h"
#include "gspan.h"

#define OPT " [-m minsup] [-x maxpat] [-i] graph-file"

int main(int argc, char **argv) {
  unsigned int minsup = 0;
  unsigned int maxpat = -1;
  bool out_instances = false;
  int opt;
  while ((opt = getopt(argc, argv, "m:x:i")) != -1) {
    switch (opt) {
    case 'm':
      minsup = atoi (optarg);
      break;
    case 'x':
      maxpat = atoi (optarg);
      break;
    case 'i':
      out_instances = true;
      break;
    default:
      std::cerr << "Usage: "<< argv[0] << OPT<< std::endl;
      return -1;
    }
  }
	
  if(argc-optind != 1){
    std::cerr << "Usage: "<< argv[0] << OPT<< std::endl;
    return -1;
  }
	
  std::ifstream graph_file(argv[optind++]);
  if(graph_file.fail()){
    std::cerr << "File not found: " << argv[optind-1] << std::endl;
    return -1;
  }
  
  Gspan gspan;
  gspan.minsup = minsup;
  if(maxpat>0){
    gspan.maxpat = maxpat;
  }
  gspan.out_instances = out_instances;
  gspan.set_data(graph_file);
  gspan.run();

  return 0;
}
