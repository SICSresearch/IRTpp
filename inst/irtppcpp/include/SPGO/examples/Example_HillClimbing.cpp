#include <Optimization/Heuristic/HillClimbing.hpp>
#include <Optimization/Heuristic/Parameter.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <math.h>

using namespace spgo;
using namespace std;

class ParameterTest : public Parameter
{
    public:
        ParameterTest(double a1, double a2)
        {
        	this->parameters.push_back(new Generic<double>(a1));
        	this->parameters.push_back(new Generic<double>(a2));
        }

        ParameterTest(){}

        virtual ~ParameterTest(){}

        Parameter * getNeighbor()
        {
            std::random_device rd;
            std::mt19937 e2(rd());
            std::normal_distribution<> dist(0, 1);
            
            Parameter* result;
            result = new ParameterTest();

            result->parameters.push_back(new Generic<double>(this->parameters[0]->getData<double>() + dist(e2) - .5));
            result->parameters.push_back(new Generic<double>(this->parameters[1]->getData<double>() + dist(e2) - .5));

            return result;
        }
};

typedef void* Info;
typedef vector<double> (*Function)(Parameter *p, Info info);

vector<double> function(Parameter *p, Info info)
{    
    vector<double> result;
    double x, y;
    double aux1, aux2;

    x = p->parameters[0]->getData<double>();
    y = p->parameters[1]->getData<double>();

    aux1 = (x+(2*y)-7);
    aux2 = ((2*x) + y - 5);

    result.push_back((aux1*aux1) + (aux2*aux2));

    return result;
}

int main()
{
    // Define the optimizer
    HillClimbing<Function, Info> * test;
    // Define the vector that hold the function and the gradient
    vector<Function> f;
    // Define the parameters
    ParameterTest p(100,100);

    // Store the function
    f.push_back(function);

    // Assign the BFGS optimizer
    test = new HillClimbing<Function, Info>();

    // Print the state before the optimization
    cout << "p = (" << p.parameters[0]->getData<double>() << ", " << p.parameters[1]->getData<double>() << ")" << endl;
    cout << "f(p) = f(" << p.parameters[0]->getData<double>() << ", " << p.parameters[1]->getData<double>() << ") = "
                        << f[0](&p, 0)[0] << endl;

    // Print the status of the optimizer
    cout << "Test status: " << ((test->run(f, p, 0) == 0) ? "Success" : "-") << endl;

    // Print the state after the optimization
    cout << "p = (" << p.parameters[0]->getData<double>() << ", " << p.parameters[1]->getData<double>() << ")" << endl;
    cout << "f(p) = f(" << p.parameters[0]->getData<double>() << ", " << p.parameters[1]->getData<double>() << ") = "
                        << f[0](&p, 0)[0] << endl;
    
    delete test;

    return 0;
}