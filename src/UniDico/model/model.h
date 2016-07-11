#ifndef MODEL_H_
#define MODEL_H_

#include <type/Matrix.h>
#include <type/parameter.h>
#include <cmath>
#include <vector>
#include <utils/andrade.h>
#include <utils/asa111.h>
#include <type/ghquads.h>

namespace irtpp
{
  class model
  {
    public:

      model(){}

      // To obtain a pointer to the static probability function
      virtual P_Function getP_Function() = 0;
      // To obtain a pointer to the static gradient function
      virtual G_Function getGrad_Function() = 0;

      virtual Boundary_Function getBoundary_Function() = 0;

      virtual Matrix<double>* getZ(int) = 0;

      virtual int getParamSize() = 0;

      virtual void transform(Matrix<double>*) = 0;

      virtual void untransform(Matrix<double>*) = 0;

      virtual void setInitialValues(Matrix<double>*, dataset* data) = 0;

      virtual void calculateError(double& max_diff, Matrix<double>* z, Matrix<double>* z_temp, int size) = 0;

      virtual void savePrevValues(Matrix<double>* z, Matrix<double>* z_temp, int size) = 0;

      static double* loglikelihood(double* z, ll_parameter param)
      {
        double tp = 0,
               tq = 0;
               param.sum[0] = 0;

        //std::cout << "b: " << z[0] << std::endl;
        param.boundary(z);
        //std::cout << "b: " << z[0] << std::endl;

        for (int k = 0; k < 40; ++k)
        {
          tp = param.probability(quads(40)[k],z);

          if (tp < 1e-08)
          {
            tp = 1e-08;
          }
          tq = 1 - tp;
          if (tq < 1e-08)
          {
            tq=1e-08;
          }

          param.sum[0] += (((*(param.r))(k,param.index))*log(tp))+(((*(param.f))(k,0))
                           - ((*(param.r))(k,param.index)))*log(tq);
        }

        param.sum[0] = -param.sum[0];
        param.LL = param.sum[0];
        return (param.sum);
      }

      virtual ~model() { delete probability; }

      Matrix<double>* probability;
      int qnodes;
  };

}
#endif
