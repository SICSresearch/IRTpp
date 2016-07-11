#ifndef RAMSAY_H_
#define RAMSAY_H_

#include <type/Matrix.h>
#include <iostream>

using namespace std;

namespace irtpp
{

  inline void ramsay(Matrix<double>** args_hist)
  {
    double dX,
           dX2,
           d2X2,
           accel,
           numerator =   0.0,
           denominator = 0.0;

    for (int i = 0; i < args_hist[0]->nR(); i++)
    {
      for (int j = 0; j < args_hist[0]->nC(); j++)
      {
        dX = (*args_hist[2])(i,j) - (*args_hist[1])(i,j);
        dX2 = (*args_hist[1])(i,j) - (*args_hist[0])(i,j);
        d2X2 = dX - dX2;

        numerator += dX * dX;
        denominator += d2X2 * d2X2;
      }
    }

    accel = 1 - sqrt(numerator / denominator);

    if (accel < -5.0)
    {
      accel = -5;
    }

    for (int i = 0; i < args_hist[0]->nR(); i++)
    {
      for (int j = 0; j < args_hist[0]->nC(); j++)
      {
        (*args_hist[2])(i,j) = (1 - accel) * (*args_hist[2])(i,j) + accel * (*args_hist[1])(i,j);
      }
    }
  }
}

#endif