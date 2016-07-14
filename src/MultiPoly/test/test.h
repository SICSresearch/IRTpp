/*
 * pimatrixtest.h
 *
 *  Created on: 25/04/2016
 *      Author: Milder
 */

#ifndef TEST_TEST_H_
#define TEST_TEST_H_

#include <vector>
#include "../util/matrix.h"

namespace irtpp {

	/**
	 * Tests the sum if each columns of pi matrix is equals to 1
	 * */
	bool test_pi ( matrix<double> &pi );

	/**
	 * The sum of matrix r has to be the number of examinees multiplied by number of items
	 *
	 * N x p
	 * */
	bool test_r ( std::vector<matrix<double> > &r, int N, int p );
	bool test_r ( matrix<double> &r, int N, int p );
}



#endif /* TEST_TEST_H_ */
