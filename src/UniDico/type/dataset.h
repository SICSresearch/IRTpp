#ifndef DATASET_H_
#define DATASET_H_

#include <map>
#include <vector>
#include <type/Matrix.h>

namespace irtpp
{

class dataset
{

public:

    int size;
    std::map<std::vector<char>, int> matrix;
    Matrix<char>* bitset_list;
    Matrix<int> * frequency_list;
    int * count_set_bits;

    //Constructor
    dataset(int size);

    // Methods
    std::map<std::vector<char>, int>::const_iterator iterator; /**use this when reading in order*/
    std::map<std::vector<char>, int>::const_iterator begin;
    std::map<std::vector<char>, int>::const_iterator end;
    inline void resetIterator(){iterator = matrix.begin();}
    inline char checkEnd(){return (iterator==matrix.end());}
    inline void iterate(){++iterator;}
    inline std::vector<char> getCurrentBitSet(){return (iterator->first);}
    inline int getCurrentFrequency(){return (iterator->second);}
    int countBitSet(Matrix<char> * bitset, int index);
    void push(std::vector<char>);/**Use this to fill the pattern matrix*/
    void push(std::vector<char>,int);/**Use this to fill the pattern matrix many times with a pattern*/
    int freq(std::vector<char>);//Frequency
    void flush();/**Use this to clean matrix*/
    void print();

    friend std::ostream& operator<< (std::ostream &, dataset &);/**Output operator*/

    //DataSet implementations
    int countItems () const;
    int countIndividuals () const;
    Matrix<char>* getBitsetList();
    Matrix<int>* getFrequencyList();

    //Destructor
    ~dataset();
};

}

#endif
