/*
 * model.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_MODEL_MODEL_H_
#define POLYTOMOUS_MODEL_MODEL_H_

#include <vector>
#include <cmath>
#include <cassert>

//including optimization files from dlib library
#include <dlib/optimization.h>

namespace irtpp {

const double LOWER_BOUND_ = 1e-08;
const double UPPER_BOUND_ = 0.999999;

namespace polytomous {

/*
 * Model class
 * It represents what is the model approach to use
 * Might be 1PL, 2PL, 3PL
 * */

class model {

typedef dlib::matrix<double,0,1> item_parameter;

public:
	/**
	 * Number of parameters of the model
	 * */
	int parameters;

	std::vector<int> *categories_item;
	int d;

	model();

	/**
	 * This receives 1, 2 or 3. Depending on the model to use
	 * */
	model(int, int, std::vector<int>*);
	virtual ~model();

	/**
	 * Probability in dichotomous case
	 * */
	double Pstar_ik(std::vector<double>&, const item_parameter&, int i, int k);

	/**
	 * This method computes the probability that a response pattern U_l has the category k to item
	 * i, given that its latent trait vector theta, and the item paramters
	 * */
	double Pik(std::vector<double>&, const item_parameter&, int i, int k);
};

}

} /* namespace irtpp */

#endif /* MODEL_MODEL_H_ */
