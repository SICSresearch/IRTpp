/*
 * model.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dicho-multi/model/model.h"

namespace irtpp {

namespace dichomulti {

model::model() {
	// TODO Auto-generated constructor stub

}

model::model(int parameters) {
	this->parameters = parameters;
}

double model::P(std::vector<double> &theta, const item_parameter &parameters) {
	if ( this->parameters == 1 ) {
		/**
		 * 1PL Approach
		 *
		 * */

		//Initialized with gamma_k value
		double eta = parameters(0);

		//Computing dot product
		for ( short i = 0; i < theta.size(); ++i )
			eta += 1 * theta[i]; //no alpha in this model

		double P = 1.0 / (1.0 + std::exp(-eta));
		P = std::max(P, LOWER_BOUND_);
		P = std::min(P, UPPER_BOUND_);
		return P;
	}
	if ( this->parameters == 2 ) {
		/**
		 * 2PL Approach
		 * */

		//Initialized with gamma value
		double eta = parameters(parameters.size() - 1);

		//Computing dot product
		for ( short i = 0; i < theta.size(); ++i )
			eta += parameters(i) * theta[i];

		double P = 1.0 / (1.0 + std::exp(-eta));
		P = std::max(P, LOWER_BOUND_);
		P = std::min(P, UPPER_BOUND_);
		return P;
	}
	/**
	 * 3PL Approach
	 * */
	double c = parameters(parameters.size() - 1);
	c = 1.0 / (1.0 + exp(-c));

	//Initialized with gamma value
	double eta = parameters(parameters.size() - 2);

	//Computing dot product
	for ( short i = 0; i < theta.size(); ++i )
		eta += parameters(i) * theta[i];

	double P = c + (1.0 - c) / (1.0 + exp(-eta));
	P = std::max(P, LOWER_BOUND_);
	P = std::min(P, UPPER_BOUND_);
	return P;
}

model::~model() {
	// TODO Auto-generated destructor stub
}

}

} /* namespace irtpp */
