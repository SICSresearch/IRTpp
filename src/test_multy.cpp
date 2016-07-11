#include <Rcpp.h>
 using namespace Rcpp;
 using std::cout;
 
 #include <iostream>
 #include "MultiPoly/simulation/simulation.h"
 
 inline void improveIO () {
	std::ios_base::sync_with_stdio(0);
	//std::cin.tie(0); std::cout.tie(0);
}

// [[Rcpp::export]]
int multiTest(std::string path){
	
	cout.setf(std::ios_base::fixed);
	cout.precision(5);

	improveIO();
	simulation sim;

	sim.run_single_dichotomous(3, 1, path, 0.0001);
	return 0;
}