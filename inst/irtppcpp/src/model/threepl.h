#ifndef THREEPL_H_
#define THREEPL_H_

#include <model/model.h>

namespace irtpp
{

  class threepl : public model
  {
    public:
      static void boundary(double* z)
      {
        if(abs(z[0]) > 5)
        {
          z[0] = 0.851;
        }
        if(z[0] < 0){
                z[0] = 0.5;
        }
        if(abs(-z[1]/z[0]) > 5)
        {
          z[1] = 0;
        }
        if(abs(z[2]) > 5)
        {
          z[2] = -1.3;
        }
      }

      Boundary_Function getBoundary_Function()
      {
        return boundary;
      }

      void transform(Matrix<double>* z)
      {
        for (int i = 0; i < z->nR(); ++i)
        {
          double qc = (*z)(i, 2);
          (*z)(i, 2) = log(qc / (1 - qc));
        }
      }

      void untransform(Matrix<double>* z)
      {
        for (int i = 0; i < z->nR(); ++i)
        {
          double qa = (*z)(i, 0);
          double qb = (*z)(i, 1);
          double qc = (*z)(i, 2);
          double ec = exp(qc);

          (*z)(i, 2) = ec / (1 + ec);
          (*z)(i, 1) = -qb / qa; //Transformacion del B   d=-b/a
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
          (*z)(i, 2) = 0.2;
        }

        delete[] result;
      }

      static double probability(double theta, double* z)
      {
        double exponential = (z[0]*theta+z[1]);

        if ( exponential > 35 )
          exponential = 35;

        else if ( exponential < -35 )
          exponential = -35;

        exponential = exp(-exponential) ;
        double ec = exp(z[2]);

        return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
      }

      P_Function getP_Function()
      {
        return probability;
      }

      G_Function getGrad_Function()
      {
        return gradient;
      }

      static double * gradient(double* z, ll_parameter param)
      {
        double p, P_Star;
        double factor;
        double ec;      // e^c_i
        double ecp1i;   // 1 / (e^c_i + 1)
        double W;

        param.gradient[0] = 0;
        param.gradient[1] = 0;
        param.gradient[2] = 0;

        ecp1i=1/(1+exp(z[2]));
        ec=exp(z[2]);

        for (int k = 0; k < 40; k++)
        {
          p = probability(quads(40)[k], z);
          P_Star = 1/(1+exp(-(z[0]*quads(40)[k]+z[1])));

          W = P_Star * ( 1 - P_Star ); // Numerator
          W /= p * ( 1 - p );// Denominator

          factor = (((*(param.r))(k,param.index)) - ((*(param.f))(k,0))*(p)) * W;

          param.gradient[0] -= factor * quads(40)[k] * ecp1i;
          param.gradient[1] -= factor * ecp1i;
          param.gradient[2] -= factor * (ec * (ecp1i*ecp1i) / P_Star);
        }

        return param.gradient;
      }

      Matrix<double>* getZ(int items)
      {
        return new Matrix<double>(items, 3);
      }

      int getParamSize()
      {
        return 3;
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
          t =        fabs((*z_temp)(i, 2) - (*z)(i, 2));
          max_diff = t > max_diff ? t : max_diff;
        }
      }

      void savePrevValues(Matrix<double>* z, Matrix<double>* z_temp, int size)
      {
        for(int i = 0; i < size; i++)
        {
          (*z_temp)(i, 0) = (*z)(i, 0);
          (*z_temp)(i, 1) = (*z)(i, 1);
          (*z_temp)(i, 2) = (*z)(i, 2);
        }
      }
  };

}

#endif
