#include <Optimization/Unrestricted/Singlevariable/Brent.hpp>
#include <iostream>

using namespace spgo;
using namespace std;

typedef double Parameter;
typedef double (*Function)(Parameter a, vector<double> info);
typedef vector< vector<double> > Info;

double function(Parameter a, vector<double> info) { return ((a+5)*(a+5)); }

int main()
{
    // Define the optimizer
    Optim<Function, Parameter, Info> * test;
    // Define the vector that hold the function
    vector<Function> f;
    // Define the parameter
    Parameter p;
    // Define the auxiliar values (bounds)
    Info info;
    // Define the vector that hold the bounds
    vector<double> bounds;
    // Define the vector that hold the constants (void)
    vector<double> constants;

    // Store the function
    f.push_back(function);
    // Store the bounds
    bounds.push_back(-10000);
    bounds.push_back(10000);
    // Store the constants and bounds vectors
    info.push_back(constants);
    info.push_back(bounds);

    // Assign the Brent optimizer
    test = new Brent<Function, Parameter, Info>();
    
    // Print the status of the optimizer
    cout << "Test status: " << ((test->run(f, p, info) == 0) ? "Success" : "-") << endl;

    // Print the state after the optimization
    cout << "p = " << p << endl;
    cout << "f(p) = f(" << p << ") = " << f[0](p, constants) << endl;

    delete test;

    return 0;
}