#ifndef HC_HPP
#define HC_HPP

#include <Optimization/Optim.hpp>
#include <Optimization/Heuristic/Parameter.hpp>
#include <cmath>

namespace spgo
{
    template <typename Function, typename Info>
    class HillClimbing : public Optim<Function, Parameter, Info>
    {
        public:
            HillClimbing(){}
            ~HillClimbing(){}
            void testttt(){}
            
            int run(std::vector<Function> f, Parameter &p, Info info)
            {
                double eval;
                Parameter * new_p;
                double new_eval;
                double change = 1000;

                eval = (f[0])(&p, info)[0];
                new_eval = eval + 100;
                
                while(change > 0.001)
                {
                    new_p = p.getNeighbor();
                    new_eval = (f[0])(new_p, info)[0];

                    if(new_eval < eval)
                    {
                        change = abs(eval - new_eval);
                        delete p.parameters[0];
                        delete p.parameters[1];
                        p.parameters[0] = new_p->parameters[0];
                        p.parameters[1] = new_p->parameters[1];
                        eval = new_eval;
                    }

                    delete new_p;
                }

                return 0;
            }

        private:
            ;
    };
}

#endif