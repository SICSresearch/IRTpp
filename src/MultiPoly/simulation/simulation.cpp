/*
 * simulation.cpp
 *
 *  Created on: 1/06/2016
 *      Author: Milder
 */

#include "simulation.h"
#define START_CLOCK		clock_t start = clock();
#define END_CLOCK		clock_t stop = clock();
#define REPORT_TIME     double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;\
						std::cout << "Time elapsed: " << elapsed << " ms." << '\n';


namespace irtpp {

simulation::simulation() {
	// TODO Auto-generated constructor stub

}

void simulation::simulate ( int model, int d, int start, int end, std::string folder,
							std::string name, double dif, bool dicho ) {
	std::ofstream report_parameters;
	std::stringstream ss;
	ss << folder << "/estimation-" << name << '-' << start << '-' << end << ".csv";
	std::string parameters = ss.str();
	report_parameters.open(parameters.c_str());
	report_parameters.precision(4);

	ss.str("");
	ss << folder << "/" << name;
	const std::string base_name = ss.str();
	for ( int i = start; i <= end; ++i ) {
		matrix<char> Y;
		input<char> in(';');

		std::stringstream ss;
		ss << base_name << i << ".csv";

		std::string file_name = ss.str();
		in.importData(file_name, Y);
		std::cout << file_name << " imported" << std::endl;

		if ( dicho ) {
			START_CLOCK

			dichomulti::estimation e(model, Y, d, dif);
			e.EMAlgortihm();

			END_CLOCK
			REPORT_TIME
			e.print_results(report_parameters, elapsed);
		} else {
			START_CLOCK

			polytomous::estimation e(model, Y, d, dif);
			e.EMAlgortihm();

			END_CLOCK
			REPORT_TIME
			e.print_results(report_parameters, elapsed);
		}
	}

	report_parameters.close();
}

void simulation::simulate ( int model, int d, int iterations, std::string folder,
							std::string name, int interval, double dif, bool dicho ) {
	for ( int i = 1; i <= iterations; i += interval ) {
		simulate(model, d, i, i + interval - 1, folder, name, dif, dicho);
	}
}

void simulation::run_single ( int model, int d, std::string filename, double dif, bool dicho ) {
	if ( dicho ) run_single_dichotomous(model, d, filename, dif);
	else		 run_single_polytomous(model, d, filename, dif);
}

void simulation::run_single_polytomous ( int model, int d, std::string filename, double dif ) {
	matrix<char> Y;
	input<char> in(';');
	in.importData(filename, Y);
	std::cout << "Data imported" << std::endl;

	START_CLOCK

	polytomous::estimation e(model, Y, d, dif);
	e.EMAlgortihm();

	END_CLOCK

	e.print_results();
	REPORT_TIME
}

void simulation::run_single_dichotomous ( int model, int d, std::string filename, double dif ) {
	matrix<char> Y;
	input<char> in(';');
	in.importData(filename, Y);
	std::cout << "Data imported" << std::endl;

	START_CLOCK

	dichomulti::estimation e(model, Y, d, dif);
	e.EMAlgortihm();

	END_CLOCK

	e.print_results();
	REPORT_TIME
}

simulation::~simulation() {
	// TODO Auto-generated destructor stub
}

} /* namespace irtpp */
