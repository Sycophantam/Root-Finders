#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
class Polynomial
{
public:
    Polynomial();
    Polynomial(vector<double> coef, bool leading);

    void print();
    double solveWithVectors(vector<double> x);
    double solveWithValue(double val);
    bool getApp(double exact, double app, double tolerance);
    vector<double> bisection(double x0=0, double x1=0);
    vector<double> falsePosition(double x0=0, double x1=0);
    vector<double> newtonRaphson(double x0=0, double x1=0);
    vector<double> secant(double x0=0, double x1=0);
    vector<double> ModifiedSecant(double x0 = 0, double x1 = 0);

    vector<double> coefficients;    //Polynomial coefficients
    bool extra = false;             //True if there is one constant

private:
    double getNormOfVector(vector<double> v);
    vector< pair<float, float> > getBoundaries(short int stop, float step);
    double findMidpoint(double beg, double end);
    void printTable(vector<double> entries);
    vector< pair<float, float> > getPossibleRoots();

};

#endif // POLYNOMIAL_H
