/*
 * estimationdata.h
 *
 *  Created on: 20/06/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_UTIL_ESTIMATIONDATA_H_
#define POLYTOMOUS_UTIL_ESTIMATIONDATA_H_
#include <vector>
#include <set>
#include "../model/model.h"
#include "../../util/matrix.h"

#include <dlib/optimization.h>

namespace irtpp {

namespace polytomous {

typedef dlib::matrix<double,0,1> item_parameter;

/**
 * Contains the information needed to execute the estimation
 * */
class estimation_data {
public:
	//Matrix of answers
	matrix<char> *dataset;
	//Dimension
	int d;
	//Matrix of response patterns
	matrix<char> Y;
	//Frequencies of each response pattern
	std::vector<int> nl;
	//Number of examines
	int N;
	//Number of response patterns
	int s;
	//Number of items
	int p;
	//Number of quadrature points
	int G;
	//Number of categories for each item
	std::vector<int> categories_item;
	//Latent traits vectors
	matrix<double> theta;
	//Weights
	std::vector<double> w;
	//Matrix r
	std::vector<matrix<double> > r;
	/**
	 * Probability matrix P
	 *
	 * P_gik means the probability that an individual has selected the category k
	 * to item i and belongs to group g
	 */
	std::vector<matrix<double> > P;
	//Matrix pi
	matrix<double> pi;
	//Pinned items (won't be estimated)
	std::set<int> pinned_items;
	//Vector or item parameters
	std::vector<item_parameter> zeta;
	//Model to use
	model m;

	estimation_data(int);
	estimation_data();
	virtual ~estimation_data();
};

}

} /* namespace irtpp */

#endif /* UTIL_ESTIMATIONDATA_H_ */
