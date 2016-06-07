#ifndef ONEPL_H_
#define ONEPL_H_

#include <model/model.h>

namespace irtpp
{

  class onepl : public model
  {
    public:
      
      static void boundary(double* z)
      {
        if(abs(z[0]) > 10)
	{
          z[0] = 0;
        }
      }

      void transform(Matrix<double>*){}
      void untransform(Matrix<double>*){}

      void setInitialValues(Matrix<double>* z, dataset* data)
      {
        double * result = Andrade(data);
        int ifault;
      
        for (int i = 0; i < data->size; i++)
	{
          (*z)(i, 0) = -(ppnd(result[0], &ifault)) / result[1];
	}

        delete[] result;
      }

      Boundary_Function getBoundary_Function()
      {
        return boundary;
      }
    
      static double probability(double theta, double* z)
      {
        double exponential = ((theta) - (z[0]));

        if (exponential > 35.0) { exponential = 35.0; }
        else if (exponential < -35.0) { exponential = -35.0; }

        return (1 / (1.0 + exp(-exponential)));
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

        param.gradient[0] = 0;

        for (int k = 0; k < 40; k++)
        {
          p = probability(quads(40)[k], z);
          param.gradient[0] += (((*(param.r))(k,param.index)) - ((*(param.f))(k,0))*(p));
        }

        return param.gradient;
      }

      Matrix<double>* getZ(int items)
      {
        return new Matrix<double>(items, 1);
      }

      int getParamSize()
      {
        return 1;
      }

      void calculateError(double& max_diff, Matrix<double>* z, Matrix<double>* z_temp, int size)
      {
        max_diff = -1;
        double t;

        for(int i = 0; i < size; i++)
        {
          t =        fabs((*z_temp)(i, 0) - (*z)(i, 0));
          max_diff = t > max_diff ? t : max_diff;
        }
      }

      void savePrevValues(Matrix<double>* z, Matrix<double>* z_temp, int size)
      {
        for(int i = 0; i < size; i++)
        {
          (*z_temp)(i, 0) = (*z)(i, 0);
        }
      }
  };

}

#endif