/*
 * LatentTraits.h
 *
 *  Created on: Feb 2, 2015
 *      Author: cesandovalp
 */

#ifndef TYPE_LATENTTRAITS_H_
#define TYPE_LATENTTRAITS_H_
#include <type/dataset.h>

namespace irtpp
{
  class LatentTraits
  {

  public:

    int             dim;
    dataset*        pm;
    Matrix<double>* traits;

    LatentTraits(dataset* p, const int dims = 1)
    {
      int rows;
        
      pm     = p;
      rows   = pm->matrix.size();
      traits = new Matrix<double>(rows, dims);
      dim    = dims;
    }

    double** getListPatternTheta()
    {
      double** result;
      Matrix<char>* pattern_list;
      int      items;

      pattern_list = pm->getBitsetList();
      result       = new double*[pm->matrix.size()];
      items        = pm->countItems();

      for(int i = 0; i < pm->matrix.size(); i++)
      {
        result[i] = new double[items + 1];

        for(int j = 0; j < items; j++)
        {
          result[i][j] = (bool)(*pattern_list)(i, j);
        }

        result[i][items] = (*traits)(i,0);
      }

      return result;
    }

    void deleteListPatternTheta(double** p)
    {
      for(int i = 0; i < pm->matrix.size(); i++) { delete[] p[i]; }

      delete[] p;
    }

    virtual ~LatentTraits() { delete traits; };

  };
}

#endif /* TYPE_LATENTTRAITS_H_ */
