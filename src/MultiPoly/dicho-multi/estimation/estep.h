/*
 * estep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef DICHOMULTI_ESTIMATION_ESTEP_H_
#define DICHOMULTI_ESTIMATION_ESTEP_H_

#include "../../util/matrix.h"
#include "../../test/test.h"
#include <ctime>

#include <iostream>

#include "../../dicho-multi/model/model.h"
#include "../../dicho-multi/type/estimationdata.h"

namespace irtpp {

namespace dichomulti {

typedef dlib::matrix<double,0,1> item_parameter;

/**
 * Estep of the EMAlgortihm
 *
 * Receives an estimation_data reference that MUST bring all the
 * data needed to run the Estep
 * */
void Estep(estimation_data&);

} /* namespace irtpp */

}

#endif /* ESTIMATION_ESTEP_H_ */
