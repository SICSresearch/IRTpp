#ifndef FUNCTION_HPP
#define FUNCTION_HPP

namespace spgo
{

template <typename Input, typename Output>
class Function
{
    public:
        
        Function(){}
        
        virtual Output run(Input) = 0;
        
        virtual ~Function(){}
};

}

#endif