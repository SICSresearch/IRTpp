#ifndef MSTEP_H_
#define MSTEP_H_

#include <Optimization/Unrestricted/Multivariable/BFGS.hpp>
#include <type/parameter.h>
#include <model/model.h>
#include <vector>

namespace irtpp
{
  typedef double* (*Func)(double * p, ll_parameter info);

  double mstep(model * m, Matrix<double> & z, m_parameter param);

}

#endif
