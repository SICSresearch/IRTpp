/*
 * initial_values.h
 *
 *  Created on: 13/06/2016
 *      Author: Milder
 */

#ifndef UTIL_INITIAL_VALUES_H_
#define UTIL_INITIAL_VALUES_H_

#include <vector>
#include <exception>
#include "matrix.h"

namespace irtpp {

	/**
	 * Computing initial values according to
	 * Andrade, Tavares & Valle (2000),
	 *
	 * */


	/**
	 * Computes initial values for Unidimensional - dichotomous case
	 *
	 * */
	void find_initial_values ( matrix<char> &data, std::vector<double> &alpha, std::vector<double> &gamma );

	/**
	 * Computes biserial correlation
	 * Needed to find a(discrimation) and b(difficulty) initial values
	 *
	 *
	 * PARAMETERS:
	 * 		IN: scores, Y
	 * 		OUT: p  ---> biserial correlation of each item
	 *
	 * */
	void compute_biserial_correlation(matrix<char> &Y, std::vector<double> &p);

	/**
	 * Receives p vector as biserial correlation
	 *
	 * and vector a where alphas will be saved
	 * */
	void compute_alphas(std::vector<double> &p, std::vector<double> &a);

	/**
	 * Receives a and p vectors (alphas and biserial correlation)
	 * and vector d where gammas will be saved
	 * */
	void compute_gammas(matrix<char> &Y, std::vector<double> &a, std::vector<double> &p, std::vector<double> &d);


	/**
	 * Helper functions to find mean of an array and standard deviation
	 * */
	inline double mean ( std::vector<double> &v );
	inline double sd ( std::vector<double> &v );
} //irtpp

#endif /* UTIL_INITIAL_VALUES_H_ */
