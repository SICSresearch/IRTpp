#ifndef TWOPL_H_
#define TWOPL_H_

#include <model/model.h>

namespace irtpp
{

  class twopl : public model
  {
    public:
      static void boundary(double* z)
      {
        // if(abs(z[0]) > 5)
        // {
        //   z[0] = 0.851;
        // }
        // double dd = -z[1]/z[0];
        // if(abs(dd)>5)
        // {
        //   z[1] = 0;
        // }
      }

      Boundary_Function getBoundary_Function()
      {
        return boundary;
      }

      void transform(Matrix<double>*){}
      void untransform(Matrix<double>* z)
      {
        for (int i = 0; i < z->nR(); ++i)
        {
          (*z)(i, 1) = -(*z)(i, 1)/(*z)(i, 0);
        }
      }

      void setInitialValues(Matrix<double>* z, dataset* data)
      {
        double * result = Andrade(data);
        int ifault;

        for (int i = 0; i < data->size; i++)
        {
          (*z)(i, 0) = std::sqrt((result[1] * result[1]) / (1.0 - result[1] * result[1]));
          (*z)(i, 1) = -(ppnd(result[0], &ifault)) / result[1];
        }

        delete[] result;
      }

      static double probability(double theta, double* z)
      {
        double exponential = (z[0] * theta) + z[1];

        if (exponential > 35) { exponential = 35; }
        else if (exponential < -35) { exponential = -35; }

        return (1 / (1 + exp(-exponential)));
      }

      P_Function getP_Function()
      {
        return probability;
      }

      G_Function getGrad_Function()
      {
        return gradient;
      }

      static double* gradient(double* z, ll_parameter param)
      {
        double p;
        double factor;

        param.gradient[0] = 0;
        param.gradient[1] = 0;

        for (int k = 0; k < 40; k++)
        {
          p = probability(quads(40)[k], z);
          factor = (((*(param.r))(k,param.index)) - ((*(param.f))(k,0))*(p));

          param.gradient[0] -= factor * quads(40)[k];
          param.gradient[1] -= factor;
        }

        return param.gradient;
      }

      Matrix<double>* getZ(int items)
      {
        return new Matrix<double>(items, 2);
      }

      int getParamSize()
      {
        return 2;
      }

      void calculateError(double& max_diff, Matrix<double>* z, Matrix<double>* z_temp, int size)
      {
        max_diff = -1;
        double t;

        for(int i = 0; i < size; i++)
        {
          t =        fabs((*z_temp)(i, 0) - (*z)(i, 0));
          max_diff = t > max_diff ? t : max_diff;
          t =        fabs((*z_temp)(i, 1) - (*z)(i, 1));
          max_diff = t > max_diff ? t : max_diff;
        }
      }

      void savePrevValues(Matrix<double>* z, Matrix<double>* z_temp, int size)
      {
        for(int i = 0; i < size; i++)
        {
          (*z_temp)(i, 0) = (*z)(i, 0);
          (*z_temp)(i, 1) = (*z)(i, 1);
        }
      }
  };

}

#endif
