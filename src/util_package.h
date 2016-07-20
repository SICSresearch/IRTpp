#ifndef UTILP_H
#define UTILP_H

#include "UniDico/type/Matrix.h"
#include "UniDico/type/parameter.h"
#include "UniDico/type/dataset.h"
#include "UniDico/type/ghquads.h"

#include "UniDico/utils/andrade.h"
#include "UniDico/utils/asa111.h"
#include "UniDico/utils/Input.h"
#include "UniDico/utils/ramsay.h"

#include "UniDico/model/model.h"
#include "UniDico/model/onepl.h"
#include "UniDico/model/twopl.h"
#include "UniDico/model/threepl.h"

#include "UniDico/estimation/emestimation.h"
#include "UniDico/estimation/estep.h"
#include "UniDico/estimation/mstep.h"
#include "UniDico/estimation/LatentTraitEstimation.h"

#include <Rcpp.h>


using namespace std;

Rcpp::List transformParameterOutput(void *);

Rcpp::List transformAbilityOutput(void *);

#endif
