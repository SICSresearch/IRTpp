/*
 * model.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "model.h"

namespace irtpp {

namespace polytomous {

model::model() {
	// TODO Auto-generated constructor stub

}

model::model(int parameters, int d, std::vector<int> *categories_item) {
	this->parameters = parameters;
	this->categories_item = categories_item;
	this->d = d;
}

double model::Pstar_ik(std::vector<double> &theta, const item_parameter &parameters, int i, int k) {
	int mi = (*categories_item)[i];
	/**
	 * Base cases
	 *
	 * Because of implementation indexes start from 0 and theory from 1
	 * It's necessary subtract 1
	 * */

	if ( k == -1 ) return 1;
	if ( k == mi - 1 ) return 0;

	if ( this->parameters == 1 ) {
		/**
		 * 1PL Approach
		 *
		 * */

		/**
		 * Initialized with gamma_k value
		 * */
		double eta = parameters(k);

		//Computing dot product
		for ( short i = 0; i < theta.size(); ++i )
			eta += 1 * theta[i]; //no alpha in this model

		return 1.0 / (1.0 + std::exp(-eta));

	}
	/**
	 * 2PL Approach
	 *
	 * */

	/**
	 * Initialized with gamma_k value
	 * */
	double eta = parameters(d + k);

	//Computing dot product
	for ( short i = 0; i < theta.size(); ++i )
		eta += parameters(i) * theta[i];

	return 1.0 / (1.0 + std::exp(-eta));
}

double model::Pik(std::vector<double> &theta, const item_parameter &parameters, int i, int k) {
	double P_ik = Pstar_ik(theta, parameters, i, k - 1) - Pstar_ik(theta, parameters, i, k);

	P_ik = std::max(P_ik, LOWER_BOUND_);
	P_ik = std::min(P_ik, UPPER_BOUND_);
	return P_ik;
}

model::~model() {
	// TODO Auto-generated destructor stub
}

}

} /* namespace irtpp */
