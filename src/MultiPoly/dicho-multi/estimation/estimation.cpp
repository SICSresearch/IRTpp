/*
 * estimation.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dicho-multi/estimation/estimation.h"

namespace irtpp {

namespace dichomulti {

estimation::estimation(int themodel, matrix<char> &dataset, short d,
					   double convergence_difference) {
	/**
	 * Object to allocate all data needed in estimation process
	 * */
	data = estimation_data(d);
	data.dataset = &dataset;


	//-------------------------------------------------------------------------------------

	//Model to be used
	model &m = data.m;

	//Number of examinees
	int &N = data.N;

	//Number of items
	int &p = data.p;

	//Number of response patterns (s <= N)
	int &s = data.s;

	//Matrix of response patterns. Its size is s x p
	matrix<char> &Y = data.Y;

	//Frequency of each pattern
	std::vector<int> &nl = data.nl;

	//Number of quadrature points
	int &G = data.G;

	//Latent trait vectors
	matrix<double> &theta = data.theta;

	//Weights
	std::vector<double> &w = data.w;

	//Matrix r. Needed in Estep and Mstep
	matrix<double> &r = data.r;

	/**
	 * Probability matrix P
	 *
	 * P_gi
	 *
	 * P_gi means the probability that an individual has selected the correct answer
	 *
	 *
	 * The purpose of this matrix is to allocate the value of P_gi
	 * to avoid recompute them while matrix Pi and r are computed in EStep
	 * */
	matrix<double> &P = data.P;

	//Matrix pi
	matrix<double> &pi = data.pi;
	//f
	std::vector<double> &f = data.f;

	//-------------------------------------------------------------------------------------


	//Matrix of response patterns and their frequency
	std::map<std::vector<char>, int> freq;
	for ( int i = 0; i < dataset.rows(); ++i )
		++freq[dataset.get_row(i)];

	Y = matrix<char>();
	nl = std::vector<int>();

	std::map<std::vector<char>, int>::iterator it;
	for ( it = freq.begin(); it != freq.end(); ++it ) {
		Y.add_row(it->first);
		nl.push_back(it->second);
	}

	N = dataset.rows();
	s = Y.rows();
	p = Y.columns(0);

	/**
	 * Number of quadrature points (G) is computed based on
	 * MAX_NUMBER_OF_QUADRATURE_POINTS and dimension of the problem, in this way
	 *
	 *
	 * G will be in 1dimension = 40 ---> 40^1 = 40
	 * 				2dimension = 20 ---> 20^2 = 400
	 * 				3dimension = 10 ---> 10^3 = 1000
	 * 				> 4dimension = 5 ---> 5^d
	 * */
	G = MAX_NUMBER_OF_QUADRATURE_POINTS / (std::min(1 << (d - 1), 8));

	// Latent trait vectors loaded from file
	theta = load_quadrature_points(d);

	// Weights loaded from file
	w = load_weights(d);

	G = theta.rows();

	//Builds r and P matrixes
	P = matrix<double>(G, p);
	pi = matrix<double>(G, s);
	r = matrix<double>(G, p);
	f = std::vector<double>(G);

	//Configurations for the estimation
	m = model(themodel);
	this->convergence_difference = convergence_difference;
	this->iterations = 0;
}

estimation::estimation(int themodel, matrix<char> &dataset, short d,
					   double convergence_difference, std::vector<int> &number_of_items) {

	estimation(themodel, dataset, d, convergence_difference);

	//Pinned items in multidimensional case (the first of each dimension)
	std::set<int> &pinned_items = data.pinned_items;

	int before = 0;
	pinned_items.insert(0);
	for ( unsigned int i = 0; i < number_of_items.size() - 1; ++i ) {
		before += number_of_items[i];
		pinned_items.insert(before);
	}

	this->convergence_difference = convergence_difference;
	this->iterations = 0;
}

void estimation::initial_values() {
	//Parameters of the items
	std::vector<item_parameter> &zeta = data.zeta;
	//Dimension
	int &d = data.d;
	//Number of examinees
	int &N = data.N;
	//Number of items
	int &p = data.p;
	//Model used in the problem
	model &m = data.m;
	//Matrix of answers of the examinees
	matrix<char> &dataset = *data.dataset;

	zeta = std::vector<item_parameter>(p);
	int total_parameters = m.parameters == 1 ? 1 : m.parameters - 1 + d;

	for ( int i = 0; i < p; ++i ) {
		zeta[i] = item_parameter(total_parameters);
		for ( int j = 0; j < total_parameters; ++j )
			zeta[i](j) = 1.0;
	}

	if ( d == 1 ) {
		std::vector<double> alpha, gamma;
		find_initial_values(dataset, alpha, gamma);

		for ( int i = 0; i < p; ++i ) {
			item_parameter &item_i = zeta[i];

			if ( m.parameters > 1 ) {
				item_i(0) = alpha[i];
				item_i(1) = gamma[i];
				if ( m.parameters == 3 ) item_i(2) = -1.1;
			} else {
				item_i(0) = gamma[i];
			}
		}
	} else {
		std::vector<double> alpha, gamma;
		find_initial_values(dataset, alpha, gamma);

		for ( int i = 0; i < p; ++i ) {
			item_parameter &item_i = zeta[i];

			if ( m.parameters < 3 ) item_i(item_i.size() - 1) = gamma[i];
			else {
				item_i(item_i.size() - 2) = gamma[i];
				item_i(item_i.size() - 1) = -1.1;
			}
		}

		//Items that will not be estimated
		std::set<int> &pinned_items = data.pinned_items;

		if ( pinned_items.empty() ) {
			int items_for_dimension = p / d;
			for ( int i = 0, j = 0; i < p; i += items_for_dimension, ++j ) {
				item_parameter &item = zeta[i];
				pinned_items.insert(i);
				for ( int h = 0; h < d; ++h )
					item(h) = 0;
				item(j) = 1;
			}
		}
	}
}

void estimation::EMAlgortihm() {
	initial_values();
	double dif = 0.0;
	do {
		Estep(data);
		dif = Mstep(data);
		++iterations;
		//std::cout << "Iteration: " << iterations << " \tMax-Change: " << dif << std::endl;
	} while ( dif > convergence_difference && iterations < MAX_ITERATIONS );
}

void estimation::print_results ( ) {
	std::vector<item_parameter> &zeta = data.zeta;
	int &p = data.p;
	model &m = data.m;

	std::cout << "Finished after " << iterations << " iterations.\n";

	bool guessing_parameter = m.parameters == 3;
	for ( int i = 0; i < p; ++i ) {
		std::cout << "Item " << i + 1 << '\n';
		for ( int j = 0; j < zeta[i].size() - guessing_parameter; ++j )
			std::cout << zeta[i](j) << ' ';
		if ( guessing_parameter ) {
			double c = zeta[i](zeta[i].size() - 1);
			std::cout << 1.0 / (1.0 + exp(-c));
		}
		std::cout << '\n';
	}
}

void estimation::print_results ( std::ofstream &fout, int elapsed ) {
	std::vector<item_parameter> &zeta = data.zeta;
	int &d = data.d;
	int &p = data.p;
	model &m = data.m;

	bool guessing_parameter = m.parameters == 3;
	for ( int i = 0; i < p; ++i ) {
		for ( int j = 0; j < zeta[i].size() - guessing_parameter; ++j ) {
			if ( j ) fout << ';';
			fout << zeta[i](j);
		}
		if ( guessing_parameter ) {
			double c = zeta[i](zeta[i].size() - 1);
			fout << 1.0 / (1.0 + exp(-c));
		}
		fout << ';' << elapsed << '\n';
	}
}

estimation::~estimation() {

}

} /* namespace dichomulti */

} /* namespace irtpp */
