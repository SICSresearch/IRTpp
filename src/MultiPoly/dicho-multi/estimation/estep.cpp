/*
 * estep.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dicho-multi/estimation/estep.h"

namespace irtpp {

namespace dichomulti {

void Estep ( estimation_data &data ) {
	//Number of items
	int &p = data.p;
	//Number of response patterns
	int &s = data.s;
	//Number of quadrature points
	int &G = data.G;
	//Model used
	model &m = data.m;
	//Matrix of response patterns
	matrix<char> &Y = data.Y;
	//Frequency of each pattern
	std::vector<int> &nl = data.nl;
	//Latent trait vectors
	matrix<double> &theta = data.theta;
	//Weights
	std::vector<double> &w = data.w;
	//Vector of parameters of the items
	std::vector<item_parameter> &zeta = data.zeta;
	//f
	std::vector<double> &f = data.f;
	f.assign(f.size(), 0);

	//pi matrix
	matrix<double> &pi = data.pi;


	// Probability matrix P
	matrix<double> &P = data.P;

	//r matrix
	matrix<double> &r = data.r;
	r.reset();

	/**
	 * Computing each element of matrix P
	 * P_gi
	 * */
	for ( int g = 0; g < G; ++g ) {
		std::vector<double> &theta_g = *theta.get_pointer_row(g);
		for ( int i = 0; i < p; ++i )
			P(g, i) = m.P(theta_g, zeta[i]);
	}

	std::vector<int> correct(p);
	int correct_size;
	//Patterns
	for ( int l = 0; l < s; ++l ) {
		double denonimator_l = 0;
		//Quadrature points
		for ( int g = 0; g < G; ++g ) {
			double &pi_gl = pi(g, l);
			pi_gl = w[g];
			correct_size = 0;
			//Items
			for ( int i = 0; i < p; ++i ) {
				if ( Y(l, i) ) {
					pi_gl *= P(g, i);
					correct[correct_size++] = i;
				}
				else
					pi_gl *= 1 - P(g, i);
			}
			/**
			 * As denominator for a response pattern l is the summation over the latent traits
			 * here pi(g, l) is added to denominator_l
			 * */
			denonimator_l += pi_gl;
		}

		for ( int g = 0; g < G; ++g ) {
			double &pi_gl = pi(g, l);
			pi_gl /= denonimator_l;

			f[g] += nl[l] * pi_gl;
			for ( int i = 0; i < correct_size; ++i )
				r(g, correct[i]) += nl[l] * pi_gl;
		}
	}

	//Asserting pi correctness
//	bool pi_ok = test_pi(pi);
//	assert(("Each column of pi matrix must sum 1.0", pi_ok));

	//Asserting r correctness
//	bool r_ok = test_r(r, data.N, p);
//	assert(("Sum of elements in r must be N x p", r_ok));
}

}

} /* namespace irtpp */
