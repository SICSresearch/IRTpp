/*
 * LatentTraitEstimation.h
 *
 *      Author: cesandovalp
 */

#ifndef LATENTTRAITESTIMATION_H_
#define LATENTTRAITESTIMATION_H_
#include <model/model.h>
#include <model/onepl.h>
#include <model/twopl.h>
#include <model/threepl.h>
#include <type/Matrix.h>
#include <type/dataset.h>
#include <type/parameter.h>
#include <estimation/estep.h>
#include <estimation/mstep.h>
#include <type/ghquads.h>
#include <type/LatentTraits.h>
#include <Types/Generic.hpp>
#include <Optimization/Unrestricted/Singlevariable/Brent.hpp>
#include <Optimization/Optim.hpp>

#include <sstream>
#include <cstdlib>
#include <iostream>

#define _DELTA 0.0001220703


namespace irtpp
{
  class LatentTraitEstimation
  {

    public:

      model*        m;
      LatentTraits* lt;

      typedef struct
      {
        bool*  pattern;
        int    node;
        model* m;
        Matrix<double>* z;
      }
      Parameter_logLP;

      typedef double (*Function)(double, Any*);

      LatentTraitEstimation() {}

      LatentTraitEstimation(dataset* d) { lt = new LatentTraits(d); }

      virtual ~LatentTraitEstimation() { delete lt; }

      inline double probabilities(bool* pattern, int size, int node)
      {
        double p = 1;

        for (int i = 0; i < size; i++)
        {
          if (pattern[i]) { p *= (*(m->probability))(node, i); }
          else { p *= 1 - (*(m->probability))(node, i); }
        }

        return (p);
      }

      static double probabilities(double theta, bool* pattern, int node, model* t_model, Matrix<double>* z)
      {
        double p;
        double temp_zita[3];

        p = 1;

        for (int i = 0; i < z->nR(); i++)
        {
          temp_zita[0] = (*z)(i, 0);
          temp_zita[1] = (*z)(i, 1);
          temp_zita[2] = (*z)(i, 2);

          if (pattern[i]) { p *= t_model->getP_Function()(theta, temp_zita); }
          else { p *= 1 - t_model->getP_Function()(theta, temp_zita); }
        }

        return (p);
      }

      LatentTraits* getLatentTraits(){ return (lt); }

      void setLatentTraits(LatentTraits* ltt) { lt = ltt; }

      void setModel(model* m) { this->m = m; }

      void estimateLatentTraitsEAP(Matrix<double>* z)
      {
        Matrix<char>* list;
        int           size,
                      counter;
        double        sum_num,
                      sum_den,
                      pp,
                      temp_zita[3];

        for (int i = 0; i < lt->pm->size; ++i)
        {
          double qc  = (*z)(i, 2);
          (*z)(i, 2) = log(qc / (1 - qc));
          (*z)(i, 1) = -(*z)(i, 1) * (*z)(i, 0);
        }

        list    = lt->pm->getBitsetList();
        size    = lt->pm->matrix.size();
        counter = 0;

        for (int pattern = 0; pattern < size; pattern++, ++counter)
        {
          bool* pattern_bool;

          sum_num      = 0;
          sum_den      = 0;
          pattern_bool = new bool[lt->pm->size];

          for(int idx = 0; idx < lt->pm->size; idx++)
          {
            pattern_bool[idx] = (*list)(pattern, idx);
          }

          for (int i = 0; i < m->qnodes; ++i)
          {
            pp = probabilities(quads(m->qnodes)[i], pattern_bool, i, this->m, z);
  
            sum_num += quads(m->qnodes)[i] * weights(m->qnodes)[i] * pp;
            sum_den += weights(m->qnodes)[i] * pp;
          }

          (*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;

          delete[] pattern_bool;
        }
      }

      void estimateLatentTraitsEAP()
      {
        Matrix<char>* list;
        int           size,
                      counter;
        double        sum_num,
                      sum_den,
                      pp;

        counter = 0;

        list = lt->pm->getBitsetList();
        size = lt->pm->matrix.size();

        for (int pattern = 0; pattern < size; pattern++, ++counter)
        {
          bool* pattern_bool;

          sum_num      = 0;
          sum_den      = 0;
          pattern_bool = new bool[lt->pm->size];

          for(int idx = 0; idx < lt->pm->size; idx++)
          {
            pattern_bool[idx] = (*list)(pattern, idx);
          }

          for (int i = 0; i < m->qnodes; ++i)
          {
            pp = probabilities(pattern_bool,
                               lt->pm->size,
                               i);

            sum_num += quads(m->qnodes)[i] * weights(m->qnodes)[i] * pp;
            sum_den += weights(m->qnodes)[i] * pp;
          }

          (*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;
        }
      }

      static double logLP(double theta, Any* param_temp)
      {
        Parameter_logLP param = param_temp->getData<Parameter_logLP>();
        double p = probabilities(theta, param.pattern, param.node, param.m, param.z);

        return (-(log(p) - ((theta * theta) / 2)));
      }

      void estimateLatentTraitsMAP(Matrix<double>* z)
      {
        Matrix<char>*    list;
        int              counter;
        vector<double>   bounds; // Define the vector that hold the bounds
        vector<Function> f;      // Define the vector that hold the function

        counter = 0;
        list    = lt->pm->getBitsetList();
        
        // Store the function
        f.push_back(logLP);
        // Store the bounds
        bounds.push_back(-5);
        bounds.push_back( 5);

        for (int index = 0; index < lt->pm->matrix.size(); index++, ++counter)
        {
          vector<Any*>     info;   // Define the auxiliar values (bounds)
          bool* pattern_bool = new bool[lt->pm->size];

          for(int idx = 0; idx < lt->pm->size; idx++)
          {
            pattern_bool[idx] = (bool)((*list)(index, idx));
          }

          spgo::Optim<Function, double, vector<Any*>>* test;
          // Define the parameter
          double p;

          Parameter_logLP f_param;

          {
            f_param.pattern = pattern_bool;
            f_param.node    = counter;
            f_param.m       = this->m;
            f_param.z       = z;
          }

          // Store the constants and bounds vectors
          info.push_back(new spgo::Generic<Parameter_logLP>(f_param));
          info.push_back(new spgo::Generic<vector<double>>(bounds));

          // Assign the Brent optimizer
          test = new spgo::Brent<Function, double, vector<Any*>>();
    
          test->run(f, p, info);

          (*lt->traits)(counter, lt->dim - 1) = p;

          delete   test;
          delete[] pattern_bool;
          delete   info[0];
          delete   info[1];
        }
      }
  };
}

#endif /* LATENTTRAITESTIMATION_H_ */
