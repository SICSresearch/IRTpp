#ifndef ANY_HPP
#define ANY_HPP

#include <string>

using namespace std;

class Any
{

  public:

    Any(){}
    virtual ~Any(){};
    template<typename T>
    T getData() { return ( *(T*) data ); }

  protected:

    void* data;
};


#endif
