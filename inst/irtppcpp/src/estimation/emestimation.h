#ifndef EMESTIMATION_H_
#define EMESTIMATION_H_

#include <model/onepl.h>
#include <model/twopl.h>
#include <model/threepl.h>
#include <type/Matrix.h>
#include <type/dataset.h>
#include <type/parameter.h>
#include <estimation/estep.h>
#include <estimation/mstep.h>
#include <type/ghquads.h>
#include <utils/ramsay.h>
#include <utils/Input.h>

namespace irtpp
{
  class emestimation
  {
    public:
      ~emestimation()
      {
      	delete m;
        delete f;
        delete r;
        delete z;
        delete z_temp;

        delete[] gradient;
        delete[] sum;
        delete[] faux;
        delete[] counter_temp;
      }

      emestimation(model* m, dataset* d)
      {
        iterations   = 0;
        m->qnodes       = 40;
        //Matrix<double> cuad(qnodes, 2);

        this->d      = d;
        m->probability  = new Matrix<double>(m->qnodes, d->size);

        items        = d->size;
        this->m      = m;
        param_size   = m->getParamSize();

        f            = new Matrix<double>(m->qnodes, 1);
        r            = new Matrix<double>(m->qnodes, items);

        z            = m->getZ(items);
        z_temp       = m->getZ(items);
        z->reset();

        gradient     = new double[param_size]{0};
        sum          = new double[1]{0};
        faux         = new double[m->qnodes];
        counter_temp = new int[d->countItems()];

        p1.f            = f;
        p1.r            = r;
        p1.probability  = m->probability;
        p1.d            = d;
        p1.faux         = faux;
        p1.counter_temp = counter_temp;

        p2.f            = f;
        p2.r            = r;
        p2.d            = d;
        p2.items        = items;
        p2.param_size   = param_size;
        p2.gradient     = gradient;
        p2.sum          = sum;
      }

      Matrix<double>* coef()
      {
        Matrix<double>* result = new Matrix<double>(items, z->nC());

        for(int i = 0; i < z->nR(); i++)
        {
          for(int j = 0; j < z->nC(); j++)
          {
            (*result)(i, j) = (*z)(i, j);
          }
        }

        return result;
      }
      
      Matrix<double>* getF(){
        Matrix<double>* result = new Matrix<double>(m->qnodes, 1);
        for(int i = 0 ; i < m->qnodes ; i++){
          (*result)(i,0) = (*f)(i,0);
        }
        return result;
      }
      
      Matrix<double>* getR(){
        Matrix<double>* result = new Matrix<double>(m->qnodes, items);
        for(int i = 0; i < m->qnodes; i++)
        {
          for(int j = 0; j < items; j++)
          {
            (*result)(i, j) = (*r)(i, j);
          }
        }
        return result;
      }

      double LogLik()
      {
        return LL;
      }

      void updateProbabilityMatrix()
      {
        for (int k = 0; k < m->qnodes; k++)
        {
          for (int i = 0; i < items; i++)
          {
            (*m->probability)(k, i) = m->getP_Function()(quads(40)[k], z->memory[i]);
          }
        }
      }

      void ** estimate()
      {
        void ** return_list = new void*[5];

        Matrix<double>* state[3];

        state[0] = m->getZ(items);
        state[1] = m->getZ(items);
        state[2] = m->getZ(items);

        state[0]->reset();
        state[1]->reset();
        state[2]->reset();

        bool convergenceSignal = false;

        m->setInitialValues(z, d);

        m->transform(z);
        for (;!(iterations++ > 500 || convergenceSignal);)
        {
          // Ramsay setup
          delete state[0];
          state[0] = state[1];
          state[1] = state[2];
          state[2] = m->getZ(items);
          m->savePrevValues(z, state[2], d->size);
          // Ramsay setup finished

          m->savePrevValues(z, z_temp, d->size);
          updateProbabilityMatrix();
          /**/
          estep(p1);
          LL = mstep(m, *z, p2);
          if(iterations > 5 && (iterations) % 3 == 0)
          {
            ramsay(state);
            m->savePrevValues(state[2], z, d->size);
          }

          /**/
          m->calculateError(max_diff, z, z_temp, d->size);
          convergenceSignal = max_diff <=  0.0001 ? true : false;
        }

        m->untransform(z);

        return_list[0] = new int(iterations);
        return_list[1] = new bool(convergenceSignal);
        return_list[2] = m->probability;
        return_list[3] = f;
        return_list[4] = r;

        delete state[0];
        delete state[1];
        delete state[2];

        return return_list;
      }


private:
      int iterations;
      double LL;
      Input input;
      Matrix<double>* f;
      Matrix<double>* r;
      Matrix<double>* z;
      Matrix<double>* z_temp;
      e_parameter     p1;
      m_parameter     p2;
      model*          m;
      dataset*        d;
      int             items;
      int             param_size;
      double          max_diff;
      double*         gradient;
      double*         sum;
      double*         faux;
      int*            counter_temp;
  };
}

#endif
