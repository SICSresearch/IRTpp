#include <Optimization/Unrestricted/Multivariable/BFGS.hpp>
#include <iostream>
#include <vector>

using namespace spgo;
using namespace std;

typedef vector<double> Parameter;
typedef void* Info;
typedef vector<double> (*Function)(Parameter p, Info info);

vector<double> function(Parameter p, Info info)
{	
	vector<double> result;
	double x, y;
	double aux1, aux2;

	x = p[0];
	y = p[1];

	aux1 = (x+(2*y)-7);
	aux2 = ((2*x) + y - 5);

	result.push_back((aux1*aux1) + (aux2*aux2));

	return result;
}

vector<double> gradient(Parameter p, Info info)
{
	vector<double> result;
	double x, y;

	x = p[0];
	y = p[1];

	result.push_back((10*x) + (8*y) - 34);
	result.push_back((8*x) + (10*y) - 38);
	
	return result;
}

int main()
{
	// Define the optimizer
	Optim<Function, Parameter, Info> * test;
	// Define the vector that hold the function and the gradient
	vector<Function> f;
	// Define the parameters
	Parameter p;

	// Store the function
	f.push_back(function);
	// Store the gradient
	f.push_back(gradient);

	// Store the initial values
	p.push_back(100);
	p.push_back(100);

	// Assign the BFGS optimizer
	test = new BFGS<Function, Parameter, Info>();

	// Print the state before the optimization
	cout << "p = (" << p[0] << ", " << p[1] << ")" << endl;
	cout << "f(p) = f(" << p[0] << ", " << p[1] << ") = " << f[0](p, 0)[0] << endl;
	cout << "g(p) = g(" << p[0] << ", " << p[1] << ") = (" << f[1](p, 0)[0] << ", " << f[1](p, 0)[1] << ")" << endl;

	// Print the status of the optimizer
	cout << "Test status: " << ((test->run(f, p, 0) == 0) ? "Success" : "-") << endl;

	// Print the state after the optimization
	cout << "p = (" << p[0] << ", " << p[1] << ")" << endl;
	cout << "f(p) = f(" << p[0] << ", " << p[1] << ") = " << f[0](p, 0)[0] << endl;
	cout << "g(p) = g(" << p[0] << ", " << p[1] << ") = (" << f[1](p, 0)[0] << ", " << f[1](p, 0)[1] << ")" << endl;
	
	delete test;

	return 0;
}