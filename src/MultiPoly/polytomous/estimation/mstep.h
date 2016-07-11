/*
 * mstep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_ESTIMATION_MSTEP_H_
#define POLYTOMOUS_ESTIMATION_MSTEP_H_

#include "../model/model.h"

#include "../../util/matrix.h"

#include "../type/estimationdata.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <set>

//including optimization files from dlib library
#include <dlib/optimization.h>

namespace irtpp {

namespace polytomous {

// Necessary typedef to be able to maximize using dlib
typedef dlib::matrix<double,0,1> column_vector;

/**
 * M step of the EM Algorithm
 *
 * Receives an estimation_data reference that MUST bring all the
 * data needed to run the Mstep
 */
double Mstep(estimation_data&);

/**
 * Log likelihood Function to maximize
 * */
class Qi {
public:
	/**
	 * Receives the number of the current item (i)
	 * and the estimation_data pointer
	 * */
    Qi (int, estimation_data*);
    //Evaluates the function
    double operator() (const column_vector&) const;
private:
    int i;
    estimation_data *data;
};

/**
 * Derivative for Log likelihood
 * */
class Qi_derivative {
public:
	/**
	 * Receives the number of the current item (i)
	 * and the estimation_data pointer
	 * */
	Qi_derivative (int, estimation_data*);
	const column_vector operator() (const column_vector&) const;
private:
    int i;
    estimation_data *data;
};

}

} /* namespace irtpp */

#endif /* ESTIMATION_MSTEP_H_ */
