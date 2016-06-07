#ifndef PARAM_H_
#define PARAM_H_

#include <type/Matrix.h>
#include <type/dataset.h>

namespace irtpp
{
  // Definition of a probability function
  // Name P_Function
  // 1st parameter, the theta value
  // 2nd parameter, the z parameters
  // return the result of the probability according to the model
  typedef double (*P_Function)(double, double*);
  typedef void (*Boundary_Function)(double*);

  struct ll_parameter
  {
    //Matrix<double>*   theta;
    Matrix<double>*   r;
    Matrix<double>*   f;
    double*           gradient;
    double*           sum;
    P_Function        probability;
    Boundary_Function boundary;
    int               index;
    double      LL;
  };

  // Definition of a probability function
  // Name G_Function
  // 1st parameter, the z parameters
  // 2nd parameter, the ll_parameters
  // return the result of the probability according to the model
  typedef double* (*G_Function)(double*, ll_parameter);

  struct e_parameter
  {
    Matrix<double>* f;
    Matrix<double>* r;
    Matrix<double>* weight;
    Matrix<double>* probability;
    dataset*        d;
    double*         faux;
    int*            counter_temp;
  };

  struct m_parameter
  {
    Matrix<double>* f;
    Matrix<double>* r;
    Matrix<double>* weight;
    Matrix<double>* theta;
    double*         gradient;
    double*         sum;
    dataset*        d;
    int             items;
    int             param_size;
    double LL;
  };

}
#endif
