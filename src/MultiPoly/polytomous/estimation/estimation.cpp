/*
 * estimation.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "estimation.h"

namespace irtpp {

namespace polytomous {

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

	//Number of categories by item
	std::vector<int> &categories_item = data.categories_item;

	//Latent trait vectors
	matrix<double> &theta = data.theta;

	//Weights
	std::vector<double> &w = data.w;

	//Matrix r. Needed in Estep and Mstep
	std::vector<matrix<double> > &r = data.r;

	/**
	 * Probability matrix P
	 *
	 * P_gik
	 *
	 * P_gik means the probability that an individual has selected the category k
	 * to item i and belongs to group g
	 *
	 *
	 * The purpose of this matrix is to allocate the value of P_gik
	 * to avoid recompute them while numerators and denominators in Estep are computed
	 * */
	std::vector<matrix<double> > &P = data.P;

	/**
	 * Matrix of probabilities pi, denominators vector and matrix of numerators
	 * needed in Estep
	 * */
	matrix<double> &pi = data.pi;

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

	//Number of categories of each item
	categories_item = std::vector<int>(p);
	for ( int j = 0; j < p; ++j ) {
		int max_category = -1;
		for ( int i = 0; i < s; ++i ) {
			//Number of categories of an item is defined as the max category found in the answers
			if ( Y(i, j) > max_category )
				max_category = Y(i, j);
		}
		categories_item[j] = max_category;
	}
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
	r = std::vector<matrix<double> >(G);
	P = std::vector<matrix<double> >(G);
	for ( int g = 0; g < G; ++g ) {
		r[g] = matrix<double>();
		P[g] = matrix<double>();
		for ( int i = 0; i < p; ++i ) {
			r[g].add_row(categories_item[i]);
			P[g].add_row(categories_item[i]);
		}
	}

	//Matrixes needed in Estep
	pi = matrix<double>(G, s);

	//Configurations for the estimation
	m = model(themodel, d, &categories_item);
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
	//Number of categories of each item
	std::vector<int> &categories_item = data.categories_item;
	//Model used in the problem
	model &m = data.m;
	//Matrix of answers of the examinees
	matrix<char> &dataset = *data.dataset;

	zeta = std::vector<item_parameter>(p);

	for ( int i = 0; i < p; ++i ) {
		int total_parameters = m.parameters == 1 ? categories_item[i] - 1 : categories_item[i] - 1 + d;
		zeta[i] = item_parameter(total_parameters);
		for ( int j = 0; j < total_parameters; ++j )
			zeta[i](j) = 1.0;
	}

	if ( d == 1 ) {
		/**
		 * Here, it is necessary find dichotomous items for each polytomous item
		 * */

		for ( int i = 0; i < p; ++i ) {
			item_parameter &item_i = zeta[i];
			int mi = categories_item[i];

			matrix<char> data_dicho(N, mi - 1);
			for ( int k = 1; k < mi; ++k ) {
				for ( int j = 0; j < N; ++j )
					data_dicho(j, k - 1) = dataset(j, i) >= k + 1;
			}

			//std::cout << data << std::endl;
			//std::cout << data_dicho << std::endl;

			std::vector<double> alpha, gamma;
			find_initial_values(data_dicho, alpha, gamma);

			/**
			 * Real alpha for this item will be the average among all alphas computed
			 * */

			double a = mean(alpha);

			if ( m.parameters > 1 ) {
				item_i(0) = a;
				for ( int k = 0; k < mi - 1; ++k )
					item_i(k + 1) = gamma[k];
			} else {
				for ( int k = 0; k < mi - 1; ++k )
					item_i(k) = gamma[k];
			}
		}

		//Here one item (with maximum number of categories) is pinned
		int item_to_pin = 0, max_categories = 0;
		for ( int i = 0; i < p; ++i ) {
			if ( data.categories_item[i] > max_categories ) {
				max_categories = data.categories_item[i];
				item_to_pin = i;
			}
		}

		data.pinned_items.insert(item_to_pin);
	}

	else {

		//TODO Compute Alphas

		/**
		 * Multidimensional case
		 * */

		int alphas = m.parameters == 2 ? d : 0;
		/**
		 * Polytomous case
		 *
		 * Here, it is necessary find dichotomous items for each polytomous item
		 * */

		for ( int i = 0; i < p; ++i ) {
			item_parameter &item_i = zeta[i];
			int mi = categories_item[i];

			matrix<char> data_dicho(N, mi - 1);
			for ( int k = 1; k < mi; ++k ) {
				for ( int j = 0; j < N; ++j )
					data_dicho(j, k - 1) = dataset(j, i) >= k + 1;
			}

			std::vector<double> alpha, gamma;
			find_initial_values(data_dicho, alpha, gamma);

			//As there is more than one gamma, it is necessary iterate over the number of categories
			for ( int k = 0; k < mi - 1; ++k )
				item_i(alphas + k) = gamma[k];
		}

		//Items that will not be estimated
		std::set<int> &pinned_items = data.pinned_items;

		/**
		 * It is supposed that there are p / d items for each dimension
		 * if the user does not specify them
		 *
		 *
		 * */

		if ( pinned_items.empty() ) {
			int items_for_dimension = p / d;
			for ( int i = 0, j = 0; i < p; i += items_for_dimension, ++j ) {
				item_parameter &item = zeta[i];
				pinned_items.insert(i);
				for ( int h = 0; h < alphas; ++h )
					item(h) = 0.0;
				item(j) = 1.0;
			}
		}
	}
}

void estimation::EMAlgortihm() {
	initial_values();
	double dif;
	do {
		Estep(data);
		dif = Mstep(data);
		++iterations;
		std::cout << "Iteration: " << iterations << " \tMax-Change: " << dif << std::endl;
	} while ( dif > convergence_difference && iterations < MAX_ITERATIONS );
}

void estimation::print_results ( ) {
	std::vector<item_parameter> &zeta = data.zeta;
	int &p = data.p;
	model &m = data.m;

	std::cout << "Finished after " << iterations << " iterations.\n";
	for ( int i = 0; i < p; ++i ) {
		std::cout << "Item " << i + 1 << '\n';
		for ( int j = 0; j < zeta[i].size(); ++j )
			std::cout << zeta[i](j) << ' ';
		std::cout << '\n';
	}
}

void estimation::print_results ( std::ofstream &fout, int elapsed ) {
	std::vector<item_parameter> &zeta = data.zeta;
	int &d = data.d;
	int &p = data.p;
	model &m = data.m;

	for ( int i = 0; i < p; ++i ) {
		for ( int j = 0; j < zeta[i].size(); ++j ) {
			if ( j ) fout << ';';
			fout << zeta[i](j);
		}
		fout << ';' << elapsed << '\n';
	}
}

estimation::~estimation() {

}

}

} /* namespace irtpp */
