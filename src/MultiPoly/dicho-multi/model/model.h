/*
 * model.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef DICHOMULTI_MODEL_MODEL_H_
#define DICHOMULTI_MODEL_MODEL_H_

#include <vector>
#include <cmath>
#include <cassert>

//including optimization files from dlib library
#include "../../include/dlib/optimization.h"

namespace irtpp {

namespace dichomulti {

const double LOWER_BOUND_ = 1e-08;
const double UPPER_BOUND_ = 0.999999;

typedef dlib::matrix<double,0,1> item_parameter;

/*
 * Model class
 * It represents what is the model approach to use
 * Might be 1PL, 2PL, 3PL
 * */

class model {

public:
	/**
	 * Number of parameters of the model
	 * */
	int parameters;

	model();

	/**
	 * This receives 1, 2 or 3. Depending on the model to use
	 * */
	model(int);
	virtual ~model();

	double P(std::vector<double>&, const item_parameter&);
};

}

} /* namespace irtpp */

#endif /* MODEL_MODEL_H_ */
