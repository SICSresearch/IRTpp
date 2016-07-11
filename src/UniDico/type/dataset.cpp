/*
* dataset.cpp
*
*  Created on: May 30, 2014
*      Author: jcliberatol
*      Updated by: cesandovalp
update again by : jcliberatol
*/

#include "dataset.h"

namespace irtpp{

  dataset::dataset(int size)
  {
    begin = matrix.begin();
    end = matrix.end();
    this->size = size;
    count_set_bits = NULL;
    bitset_list = NULL;
    frequency_list = NULL;
  }

  dataset::~dataset()
  {
    delete bitset_list;
    delete frequency_list;
    delete[] count_set_bits;
  }

  int dataset::countItems() const { return ((matrix.empty()) ? 0 : size); }

  int dataset::freq(std::vector<char> bitset) { return (matrix[bitset]); }

  int dataset::countIndividuals() const
  {
    std::map<std::vector<char>, int>::const_iterator it;
    int counter = 0;

    for (it = matrix.begin(); it != matrix.end(); ++it)
    counter += it->second;

    return (counter);
  }

  void dataset::push(std::vector<char> n) { matrix[n]++; }

  void dataset::push(std::vector<char> n, int k) { matrix[n] = matrix[n] + k; }

  void dataset::flush() { matrix.clear(); }

  void dataset::print()
  {
    for (iterator = matrix.begin(); iterator != matrix.end(); ++iterator)
    {
      for (int var = 0; var < size; ++var)
      {
        //std::cout << iterator->first[var];
      }

      //std::cout << " " << iterator->second << std::endl;
    }
  }

  int dataset::countBitSet(Matrix<char>* bitset, int index)
  {
    //When the pointer is not set yet, set it by flipping all to -1, this means they are not yet counted
    if (count_set_bits == NULL)
    {
      count_set_bits = new int[matrix.size()];
      for (int i = 0; i < matrix.size(); i++)
      {
        count_set_bits[i] = -1;
      }
    }

    //So, the pointer was initialized but the index is not counted for this entry
    //Then count it with the function.
    //TODO : What happens when the matrix changes = ?
    if (count_set_bits[index] == -1)
    {
      count_set_bits[index] = 0;
      for (int i = 0; i < size; i++)
      {
        if ((*bitset)(index,i))
        {
          count_set_bits[index]++;
        }
      }
    }

    return (count_set_bits[index]);
  }

  Matrix<char>* dataset::getBitsetList()
  {
    if (bitset_list == NULL)
    {
      std::map<std::vector<char>, int>::const_iterator it;
      std::map<std::vector<char>, int>::const_iterator begin = matrix.begin();
      std::map<std::vector<char>, int>::const_iterator end = matrix.end();

      int msize = matrix.size();
      
      bitset_list = new Matrix<char>(msize,size);
      frequency_list = new Matrix<int>(msize,1);

      int counter = 0;
      int j;

      for (it = begin; it != end; ++it, ++counter)
      {
        j = 0;

        for (auto it2 = it->first.begin(); it2 != it->first.end(); ++it2, ++j)
        {
          (*bitset_list)(counter, j) = *it2;
        }

        (*frequency_list)(counter,0) = it->second;
      }
    }

    return (bitset_list);
  }

  Matrix<int> * dataset::getFrequencyList()
  {
    if (bitset_list == NULL)
    {
      getBitsetList();
    }

    return (frequency_list);
  }

}
