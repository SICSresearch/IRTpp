/*
 * estimation.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef DICHOMULTI_ESTIMATION_ESTIMATION_H_
#define DICHOMULTI_ESTIMATION_ESTIMATION_H_

#include "../../util/initial_values.h"
#include "../../util/matrix.h"
#include "../../util/quadraturepoints.h"

#include <map>
#include <cmath>

#include "../estimation/estep.h"
#include "../estimation/mstep.h"
#include "../model/model.h"
#include "../type/estimationdata.h"

namespace irtpp {

namespace dichomulti {
/**
 * Max number of iterations of EMAlgorithm
 * */
const int MAX_ITERATIONS = 500;

const int MAX_NUMBER_OF_QUADRATURE_POINTS = 40;

/*
 * Class to set up and run the estimation process
 *
 * The main method is EMAlgorithm
 * */

class estimation {
	private:

		/**
		 * Counts the actual number of iterations
		 * */
		short iterations;

		/**
		 * Epsilon to stop the EMAlgorithm
		 * */
		double convergence_difference;

		/**
		 * Saves all data need in the estimation process
		 * */
		estimation_data data;

	public:
		/**
		 * Receives:
		 * 	1. A specific model to use -> [1, 3] that means 1PL, 2PL, or 3PL
		 *  2. A matrix containing the answers of examinees
		 *  3. The dimension of the problem
		 *	4. The epsilon (convergence difference) that the algoritm will use
		 *		as a criterion to stop
		 *
		 *
		 *  Then it sets up all the data needed to start the estimation process
		 *
		 * */
		estimation(int, matrix<char>&, short, double);


		/**
		 * Constructor of multidimensional case
		 *
		 * Receives:
		 * 	1. A specific model to use -> [1, 3] that means 1PL, 2PL, or 3PL
		 *  2. A matrix containing the answers of examinees
		 *  3. The dimension of the problem
		 *	4. The epsilon (convergences difference) that the algoritm will use
		 *		as a criterion to stop
		 *	5. A vector that stores the number of items that each dimension has
		 *
		 *	   Vector's size MUST BE equal to the dimension
		 *
		 *
		 *  Then it sets up all the data needed to start the estimation process
		 *
		 * */
		estimation(int, matrix<char>&, short, double, std::vector<int>&);

		virtual ~estimation();

		/**
		 * Finds the initial values for every parameter of the items to start the estimation
		 * */
		void initial_values();
		/*
		 * Runs the EMAlgorithm to find out the parameters
		 * */
		void EMAlgortihm();

		/**
		 * Prints the results
		 * Values of the estimated parameters
		 * */
		void print_results();

		/**
		 * Prints the results to a specific output stream,
		 * including time elapsed in ms
		 * */
		void print_results(std::ofstream &, int);
};

}

} /* namespace irtpp */

#endif /* ESTIMATION_ESTIMATION_H_ */
