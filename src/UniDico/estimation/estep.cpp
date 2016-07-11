#include <estimation/estep.h>
#include <type/parameter.h>
#include <type/ghquads.h>

namespace irtpp
{

  void estep(e_parameter param)
  {
    double sum;

    int i,
        k,
        items,
        counter_set;

    Matrix<char>* bitset_list;
    Matrix<int>*  frequency_list;

    items          = param.d->countItems();
    bitset_list    = param.d->getBitsetList();
    frequency_list = param.d->getFrequencyList();

    param.r->reset();
    param.f->reset();

    for (int pattern = 0; pattern < param.d->matrix.size(); pattern++)
    {
      sum = 0.0;
      //Calculate g*(k) for all the k's
      //first calculate the P for each k and store it in the array f aux
      for (k = 0; k < 40; k++)
      {
        param.faux[k] = weights(40)[k];
        //Calculate the p (iterate over the items in the productory)
        counter_set   = 0;
        for (i = 0; i < items; i++)
        {
          if ((*bitset_list)(pattern,i))
          {
            param.counter_temp[counter_set++] = i + 1;
            param.faux[k] *= (*(param.probability))(k,i);
          }
          else
          {
            param.faux[k] *= 1 - (*(param.probability))(k,i);
          }
        }
        //At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
        //Now multiply by the weight
        sum += param.faux[k];
      }

      for (k = 0; k < 40; k++)
      {
        param.faux[k] *= ((*frequency_list)(pattern,0)) / sum; //This is g*_j_k
        (*(param.f))(k,0) += param.faux[k];

        for (i = 0; i < counter_set; i++)
        {
          (*(param.r))(k, param.counter_temp[i] - 1) += param.faux[k];
        }
      }
    }
  }
}
