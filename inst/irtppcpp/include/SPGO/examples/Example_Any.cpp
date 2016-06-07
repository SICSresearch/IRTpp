#include <iostream>
#include <vector>
#include <map>
#include <Types/Generic.hpp>

using namespace std;
using namespace spgo;

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

int main()
{
    map<string, Any*> test;
    Parameter v;
    Info i;

    test["1st"] = (new Generic<int>(1015446513));
    test["2nd"] = (new Generic<double>(3.141592654));
    test["3rd"] = (new Generic<string>("lel"));
    test["4th"] = (new Generic<Parameter>(v));
    test["5th"] = (new Generic<Info>(i));
    test["6th"] = (new Generic<Function>(function));

    cout << test["1st"]->getData<int>() << endl;
    cout << test["2nd"]->getData<double>() << endl;
    cout << test["3rd"]->getData<string>() << endl;

    delete test["1st"];
    delete test["2nd"];
    delete test["3rd"];
    delete test["4th"];
    delete test["5th"];
    delete test["6th"];

    return 0;
}