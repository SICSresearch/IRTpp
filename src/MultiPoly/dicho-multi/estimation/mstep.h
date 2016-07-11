/*
 * mstep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef DICHOMULTI_ESTIMATION_MSTEP_H_
#define DICHOMULTI_ESTIMATION_MSTEP_H_

#include "../../util/matrix.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <set>

#include "../../dicho-multi/model/model.h"
#include "../../dicho-multi/type/estimationdata.h"

namespace irtpp {

namespace dichomulti {

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
    double operator() (const item_parameter&) const;
private:
    int i;
    estimation_data *data;
};

}

} /* namespace irtpp */

#endif /* ESTIMATION_MSTEP_H_ */
