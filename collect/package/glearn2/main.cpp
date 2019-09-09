// main.cpp
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include "glearn.h"
#include "lossfunc.h"
#include "getopt.h"

#define OPT " [-1 alpha] [-2 beta] [-i] graph-file"

int main(int argc, char **argv) {
  bool out_instances = false;
  double lm1 = 1.0;
  double lm2 = 0.0;
  double scale = 1.0;
  int maxpat = std::numeric_limits<int>::max();
  int maxitr = 1000;
  
  bool prune = true;
  bool unique = true;
  bool inexact = false;
  bool continuation = false;
  
  std::string out_fname = "model.out";
  
  static struct option long_options[] =
    {
      {"outfile", required_argument, 0, 0}, // index = 0
      {"maxpat", required_argument, 0, 0},  // index = 1
      {"maxitr", required_argument, 0, 0},  // index = 2
      {"noprune", no_argument, 0, 0},       // index = 3
      {"noequiv", no_argument, 0, 0},       // index = 4
      {"inexact", no_argument, 0, 0},       // index = 5
      {"continuation", no_argument, 0, 0}   // index = 6
    };
  
  int option_index = 0;
  int opt;
  
  while ((opt = getopt_long(argc, argv, "m:1:2:s:i",
			    long_options, &option_index)) != -1) {
    switch (opt) {
    case 0:
      if (long_options[option_index].flag != 0)
	break;
      if(option_index==0){
	out_fname = optarg;
      }else if(option_index==1){
	maxpat = atoi(optarg);
	std::cout << "maxpat set to " << maxpat << std::endl;
      }else if(option_index==2){
	maxitr = atoi(optarg);
	std::cout << "maxitr set to " << maxitr << std::endl;
      }else if(option_index==3){
	prune = false;
	std::cout << "No pruning option selected." << std::endl;
      }else if(option_index==4){
	unique = false;
	std::cout << "No processing for equivalent class." << std::endl;
      }else if(option_index==5){
	inexact = true;
	std::cout << "inexact pruning." << std::endl;
      }else if(option_index==6){
	continuation = true;
	std::cout << "continuation (warm start)" << std::endl;
      }
      break;
    case '1':
      lm1 = atof (optarg);
      break;
    case '2':
      lm2 = atof (optarg);
      break;
    case 's':
      scale = atof (optarg);
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
  
  Glearn gl;
  gl.lm1 = lm1;
  gl.lm2 = lm2;
  gl.maxpat= maxpat;
  gl.maxitr= maxitr;
  gl.scaling = scale; // standardize?
  gl.gamma = 0.0;
  gl.sigma = 0.1;
  gl.c     = 0.5;
  gl.alpha = 1.0;
  gl.gs_cnst = 0.9;
  gl.tol = 1e-6;
  //gl.tol = 1e-3;
  gl.prune  = prune;
  gl.unique = unique;
  gl.inexact = inexact;
  gl.continuation = continuation;
  gl.out_instances = out_instances;
  gl.set_data(graph_file);
  
  gl.ofname = out_fname;
  
  gl.run();
  
  std::cout << "output result into a file..." << std::endl;
  
  std::ofstream outfile(out_fname.c_str());
  
  outfile.setf(std::ios::fixed,std::ios::floatfield);
  outfile.precision(12);
  
  outfile << gl.param.theta[0] << std::endl;
  //for(int s=1; s <= (int) gl.param.get_num_vars(); ++s){
  for(set<int>::iterator it = gl.param.nonzero_idx.begin();
      it != gl.param.nonzero_idx.end(); ++it){
    outfile << gl.param.theta[*it]*scale << " " << gl.param.retrieve(*it) << std::endl;
  }
  
  return 0;
}
