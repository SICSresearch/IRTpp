/*
 * quadrature_points.h
 *
 *  Created on: 22/04/2016
 *      Author: Milder
 */

#ifndef SRC_UTIL_QUADRATUREPOINTS_H_
#define SRC_UTIL_QUADRATUREPOINTS_H_

#include <cstring>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include "matrix.h"
#include <vector>
#include <cmath>

namespace irtpp {
	/**
	 * Functions to compute and save quadrature points and weights
	 * based on gauss.quad() function from statmod library in R
	 * */

	void compute_and_save_quadrature_points(int, int);
	void compute_and_save_weights(int, int);


	/**
	 * Functions to load quadrature points and weights previously computed
	 * by functions above
	 * */

	matrix<double> load_quadrature_points(int);
	std::vector<double> load_weights(int);
}



#endif /* SRC_UTIL_QUADRATUREPOINTS_H_ */
