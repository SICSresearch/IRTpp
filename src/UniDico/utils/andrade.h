#ifndef ANDRADE_H_
#define ANDRADE_H_

#include <type/dataset.h>
#include <cmath>

namespace irtpp
{
  static double* Andrade(dataset* data)
  {
    int pSize, index;
    double Ni, frequencyV, PII, corr, mT;
    double sdT, sdU, covar, mUU, mTU, mU;
    double* T;
    double* U;
    double* TU;
    double* UU;
    double* Tm;
    double* Um;
    double* result;
    Matrix<int>* frequency_list;
    Matrix<char>* bitset_list;

    pSize = data->matrix.size();
    Ni = data->countIndividuals();
    frequency_list = data->getFrequencyList();
    bitset_list = data->getBitsetList();

    T = new double[pSize];
    U = new double[pSize];
    TU = new double[pSize];
    UU = new double[pSize];
    Tm = new double[pSize];
    Um = new double[pSize];

    PII = corr = 0;

    for (int i = 0; i < data->size; i++)
    {
      PII = mT = mU = mTU = mUU = 0.0;
      for (index = 0; index < pSize; index++)
      {
        frequencyV = (*frequency_list)(index,0);
        T[index] = 0;
        T[index] = data->countBitSet(bitset_list, index);
        PII += frequencyV * (*bitset_list)(index,i);
        U[index] = (*bitset_list)(index,i);
        TU[index] = T[index] * U[index];
        UU[index] = U[index] * U[index];
        mT += frequencyV * T[index];
        mU += frequencyV * U[index];
        mTU += frequencyV * TU[index];
        mUU += frequencyV * UU[index];
      }

      PII /= Ni;
      mT /= Ni;
      mU /= Ni;
      mTU /= Ni;
      mUU /= Ni;
      covar = mTU - mU * mT;
      sdT = sdU = 0.0;

      for (index = 0; index < pSize; index++)
      {
        frequencyV = (*frequency_list)(index,0);
        Tm[index] = T[index] - mT;
        Um[index] = U[index] - mU;
        sdT += frequencyV * Tm[index] * Tm[index];
        sdU += frequencyV * Um[index] * Um[index];
      }

      sdT = std::sqrt(sdT / (Ni - 1.0));
      sdU = std::sqrt(sdU / (Ni - 1.0));
      corr = covar / (sdT * sdU);
    }

    result = new double[2]{PII, corr};

    delete[] T;
    delete[] U;
    delete[] TU;
    delete[] UU;
    delete[] Tm;
    delete[] Um;

    return (result);
  }
}

#endif