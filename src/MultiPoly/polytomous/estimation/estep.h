/*
 * estep.h
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#ifndef POLYTOMOUS_ESTIMATION_ESTEP_H_
#define POLYTOMOUS_ESTIMATION_ESTEP_H_

#include "../model/model.h"

#include "../../util/matrix.h"
#include "../type/estimationdata.h"

#include "../../test/test.h"
#include <ctime>

#include <iostream>

namespace irtpp {

namespace polytomous {


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
