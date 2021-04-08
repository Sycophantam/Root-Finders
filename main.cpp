#include <iostream>
#include "polynomial.h"
using namespace std;

int main()
{
    Polynomial e(vector<double>{1, 2, 10, -20}, true);
    cout << "Bisection" << endl;
    e.bisection(2,1);
    cout << endl << endl;
    cout << "False-position" << endl;
    e.falsePosition(2,1);
    cout << endl << endl;
    cout << "Newton-Raphson" << endl;
    e.newtonRaphson(2,1);
    cout << endl;
    cout << "Secant" << endl;
    e.secant(2, 1);
    cout << endl;
    cout << "Modified Secant" << endl;
    e.ModifiedSecant(2, 1);
    cout << endl;

//    cout << "Second equation" << endl;
//    Polynomial e2;
//    e2.bisection();
//    cout << endl;
//    cout << "False-position" << endl;
//    e2.falsePosition();
//    cout << endl;
//    cout << "Newton-Raphson" << endl;
//    e2.newtonRaphson();
//    cout << endl;
//    cout << "Secant" << endl;
//    e2.secant();
//    cout << endl;
//    cout << "Modified Secant" << endl;
//    e2.ModifiedSecant();
//    cout << endl;
    return 0;
}
