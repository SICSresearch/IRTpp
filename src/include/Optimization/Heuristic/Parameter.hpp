#ifndef PARAMETER_HPP
#define PARAMETER_HPP

#include <Types/Generic.hpp>
#include <vector>

namespace spgo
{
    class Parameter
    {
        public:
            Parameter(){}
            virtual ~Parameter()
            {
                //ToDo;
            }
            virtual Parameter * getNeighbor() = 0;

            std::vector<Any*> parameters;
    };
}

#endif