#include "polynomial.h"

Polynomial::Polynomial()
{

}

Polynomial::Polynomial(vector<double> coef, bool leading)
{
    extra = leading;
    coefficients = coef;
}

double Polynomial::solveWithVectors(vector<double> x)
{
    double total = 0;
    std::vector<double>::size_type to = coefficients.size();
    if(extra)
    {
        total += coefficients.back();
        to = coefficients.size() - 1;
    }
    for(unsigned u = 0; u < to; u++)
    {
        total += coefficients.at(u) * x.at(u);
    }

    return total;
}

double Polynomial::solveWithValue(double val)
{
    double sum = coefficients.at(0);
    for(unsigned i = 1; i < coefficients.size(); i++)
    {
//        cout << sum << "*" << val << " = " << sum*val << endl;
//        cout << sum*val << " + " << coefficients.at(i) << " = " << sum*val + coefficients.at(i) << endl;
        sum = sum*val + coefficients.at(i);
    }
    //This extra bit accounts for whether or not the final term is a constant
    //or if it's multiplied by an x
    if(!extra)
        sum*=val;
    return sum;
}

void Polynomial::print()
{
    cout << coefficients.at(0) << "x1";
    for(unsigned u = 1; u < coefficients.size() - 1; u++)
    {
        cout << " + " << coefficients.at(u) << "x" << u + 1;
    }
    cout << " + " << coefficients.at(coefficients.size() - 1);
    if(!extra)
        cout << "x" << coefficients.size();
    cout << endl;
}

bool Polynomial::getApp(double exact, double app, double tolerance)
{
    return abs(exact - app)/exact < tolerance;
}

vector<double> Polynomial::bisection(double x0, double x1)
{
    //Stopping conditions:
    //  If there is no root found between -1000 and 1000
    //  If more than 100 iterations pass without finding the root
    //

    //Choosing the boundaries
    vector< pair<float, float> > boundaries;
    boundaries = getPossibleRoots();
    vector<double> roots;

    //The number of roots discovered is the size of boundaries
    for(unsigned u = 0; u < boundaries.size(); u++)
    {
        pair<float, float> curr = boundaries.at(u);
        cout << "Root is between " << curr.first << " and "
             << curr.second << endl;
        //This means that we found a root naturally
        if(abs(static_cast<double>(curr.first - curr.second)) < .0001)
        {
            roots.push_back(static_cast<double>(curr.first));
        }
        else
        {
            double epsilon = 1;
            cout << "|\tn\t" << "|\ta\t" << "|\tb\t" << "|\tc\t" << "|\tf(a)\t"
                 << "|\tf(b)\t" << "|\tf(c)\t" << "|\te\t|"
                 << endl;
            double a = static_cast<double>(curr.first);
            double b = static_cast<double>(curr.second);
            double c = 0;
            double c2 = 0;
            int count = 0;
            unsigned iter = 0;
            while(abs(epsilon) > .001 && count < 100)
            {                
                c2 = static_cast<double>(findMidpoint(a,b));
                double fa, fb, fc2 = 0;

                //To test the second equation,
                //we want to make an empty Polynomial
                if(coefficients.size() != 0)
                {
                    fa = solveWithValue(a);
                    fb = solveWithValue(b);
                    fc2 = solveWithValue(c2);
                }
                else
                {
                    fa = a + 10 - a * cosh(50/a);
                    fb = b + 10 - b * cosh(50/b);
                    fc2 = c2 + 10 - c2 * cosh(50/c2);
                }

                epsilon = abs(c2 - c)/c2;
                cout << setprecision(0);
                cout << "|\t" << iter << "\t|" << setw(8) << setprecision(6)
                     << fixed << a << "\t|" << setw(8) << setprecision(6)
                     << fixed << b << "\t|" << setw(8) << setprecision(6)
                     << fixed << c2 << "\t|" << setw(8) << setprecision(6)
                     << fixed << fa << "\t|" << setw(8) << setprecision(6)
                     << fixed << fb << "\t|" << setw(8) << setprecision(6)
                     << fixed << fc2 << "\t|" << setw(8) << setprecision(6)
                     << fixed << epsilon << "\t|" << endl;
                if((fc2 * fa) < 0)
                {
                    b = c2;
                }
                else if((fc2 * fa) > 0)
                {
                    a = c2;
                }
                else if (abs(fc2 * fa) < .0001)
                {
                    break;
                }
                c = c2;
                count++;
                iter++;
            }
            cout << "Root is at: " << c2 << endl;
            roots.push_back(c2);
        }
        cout << endl;
    }
    return roots;
}

vector<double> Polynomial::falsePosition(double x0, double x1)
{
    //Stopping conditions:
    //  If there is no root found between -1000 and 1000
    //  If more than 100 iterations pass without finding the root
    //

    //Choosing the boundaries
    vector< pair<float, float> > boundaries;
    boundaries = getPossibleRoots();
    if(x0 == x1 && x0 == 0)
    {
        //Choosing the boundaries

        boundaries = getPossibleRoots();

    }
    else
    {
        boundaries.push_back({x0,x1});
    }
    vector<double> roots;

    //The number of roots discovered is the size of boundaries
    for(unsigned u = 0; u < boundaries.size(); u++)
    {
        pair<float, float> curr = boundaries.at(u);
        cout << "Root is between " << curr.first << " and "
             << curr.second << endl;
        //This means that we found a root naturally
        if(abs(static_cast<double>(curr.first - curr.second)) < .0001)
        {
            roots.push_back(static_cast<double>(curr.first));
        }
        else
        {
            double epsilon = 1;
            cout << "|\tn\t" << "|\ta\t" << "|\tb\t" << "|\tf(a)\t"
                 << "|\tf(b)\t"<< "|\tc\t" << "|\tf(c)\t" << "|\te\t|" << endl;
            double a = static_cast<double>(curr.first);
            double b = static_cast<double>(curr.second);
            double c = 0;
            double c2 = 0;
            int count = 0;
            unsigned iter = 0;
            while(abs(epsilon) > .001 && count < 100)
            {
                double fa, fb, fc2 = 0;

                //To test the second equation,
                //we want to make an empty Polynomial
                if(coefficients.size() != 0)
                {
                    fa = solveWithValue(a);
                    fb = solveWithValue(b);
                    c2 = ((a * fb) - (b * fa))/(fb - fa);
                    fc2 = solveWithValue(c2);
                }
                else
                {
                    fa = a + 10 - a * cosh(50/a);
                    fb = b + 10 - b * cosh(50/b);
                    c2 = ((a * fb) - (b * fa))/(fb - fa);
                    fc2 = c2 + 10 - c2 * cosh(50/c2);
                }

                epsilon = abs(c2 - c)/c2;

                cout << setprecision(0);
                cout << "|\t" << iter << "\t|" << setw(8) << setprecision(6)
                     << fixed << a << "\t|" << setw(8) << setprecision(6)
                     << fixed << b << "\t|" << setw(8) << setprecision(6)
                     << fixed << fa << "\t|" << setw(8) << setprecision(6)
                     << fixed << fb << "\t|" << setw(8) << setprecision(6)
                     << fixed << c2 << "\t|" << setw(8) << setprecision(6)
                     << fixed << fc2 << "\t|" << setw(8) << setprecision(6)
                     << fixed << epsilon << "\t|" << endl;

                if(fc2 < 0)
                {
                    a = c2;
                }
                else if(fc2 > 0)
                {
                    b = c2;
                }
                else if (abs(fc2) < .0001)
                {
                    break;
                }
                c = c2;
                count++;
                iter++;
            }
            cout << "Root is at: " << c2 << endl;
            roots.push_back(c2);
        }
        cout << endl;
    }
    return roots;
}

vector<double> Polynomial::newtonRaphson(double x0, double x1)
{
    //Stopping conditions:
    //  If there is no root found between -1000 and 1000
    //  If more than 100 iterations pass without finding the root
    //

    //Choosing the boundaries
    vector< pair<float, float> > boundaries;
    if(x0 == x1 && x0 == 0)
    {
        //Choosing the boundaries
        boundaries = getPossibleRoots();
    }
    else
    {
        boundaries.push_back({x0,x1});
    }
    vector<double> roots;

    //The number of roots discovered is the size of boundaries
    for(unsigned u = 0; u < boundaries.size(); u++)
    {
        pair<float, float> curr = boundaries.at(u);
        //This means that we found a root naturally
        double epsilon = 1;
        cout << "|\tn\t" << "|\tx_i\t" << "|\tf(x_i)\t" << "|\tf'(x_i)\t"
             << "|\tx_i+1\t" << "|\te\t|" << endl;
        double a = static_cast<double>(curr.first);
        double b = static_cast<double>(curr.second);

        //Initial guess
        double xi = a;
        double xip1 = 0;

        //If the method fails, we want to try a new starting point. Our next
        //point wants to be in the opposite direction of where we started.
        double start = xi;
        //Variation variable. If the method fails, we start again with a
        //slight offset
        double var = .1;
        int count = 0;
        unsigned iter = 0;
        bool again = true;
        while(again)
        {
            while(abs(epsilon) > .001 && count < 100)
            {
                again = false;
                double fxi, dfxi, fxip1 = 0;

                //To test the second equation,
                //we want to make an empty Polynomial
                if(coefficients.size() != 0)
                {
                    fxi = solveWithValue(xi);
                    dfxi = 3 * pow(xi,2) + 4 * xi + 10;
                    xip1 = xi - fxi/dfxi;

                }
                else
                {
                    fxi = xi + 10 - xi * cosh(50/xi);
                    dfxi = ((50 * sinh(50/xi))/xi) - cosh(50/xi) + 1;
                    xip1 = xi - fxi/dfxi;
                }
                //If the guess goes outside the boundary
//                if(xip1 > b || xip1 < a)
//                {
//                    again = true;
//                    if(xip1 < start)
//                        xi = start + var;
//                    else
//                        xi = start - var;
//                    var +=.1;
//                    break;
//                }



                epsilon = abs(xip1 - xi)/xip1;

                cout << setprecision(0);
                cout << "|\t" << iter << "\t|" << setw(8) << setprecision(5)
                     << fixed << xi << "\t|" << setw(8) << setprecision(5)
                     << fixed << fxi << setw(8) << "\t|" << setw(8)
                     << setprecision(5) << fixed << dfxi << "\t|" << setw(8)
                     << setprecision(5) << fixed << xip1 << "\t|" << setw(8)
                     << setprecision(5) << fixed << epsilon << "\t|"
                     << endl;

                if(coefficients.size() != 0)
                {
                    fxip1 = solveWithValue(xip1);
                }
                else
                {
                    fxip1 = xip1 + 10 - xip1 * cosh(50/xip1);
                }

                //Testing whether we actually got to the root we wanted
                if(abs(epsilon) < .001)
                {
                    cout << "Infinite loop?" << endl;
                    //If the error is small but f(xip1) is not zero, we
                    //are oscillating around a point of zero slope
                    if(fxip1 > .0001)
                    {
                        again = true;
                        cout << "Oscillating detected. Going again.\n";

                        //Finding out which way to go to find another start
                        //point

                        //If the position is less than where we started, we
                        //want to go to the right
                        if(xip1 <= start)
                            xi = start + var;
                        else
                            xi = start - var;
                        var += .1;
                        break;
                    }
                }
                //cout << "e" << endl;
                xi = xip1;
                count++;
                iter++;
                if(count == 100)
                {
                    again = false;
                    break;
                }
            }
        }
        cout << "Root is at: " << xip1 << endl;
        roots.push_back(xip1);

        cout << endl;
    }
    return roots;
}

vector<double> Polynomial::secant(double x0, double x1)
{
    //Stopping conditions:
    //  If there is no root found between -1000 and 1000
    //  If more than 100 iterations pass without finding the root
    //

    vector< pair<float, float> > boundaries;
    if(x0 == x1 && x0 == 0)
    {
        //Choosing the boundaries

        boundaries = getPossibleRoots();

    }
    else
    {
        boundaries.push_back({x0,x1});
    }
    vector<double> roots;

    //The number of roots discovered is the size of boundaries
    for(unsigned u = 0; u < boundaries.size(); u++)
    {
        pair<float, float> curr = boundaries.at(u);
        double epsilon = 1;
        cout << "|\tn\t" << "|\tx_i-1\t" << "|\tf(x_i-1)\t" << "|\tx_i\t"
             << "|\tf(x_i)\t" << "|\tx_i+1\t" << "|\te\t|" << endl;
        double a = static_cast<double>(curr.first);
        double b = static_cast<double>(curr.second);

        //Initial guess
        double xim1 = a;
        double xi = b;
        double xip1 = 0;

        //If the method fails, we want to try a new starting point. Our next
        //point wants to be in the opposite direction of where we started.
        double start1 = xim1;
        double start2 = xi;
        //Variation variable. If the method fails, we start again with a
        //slight offset
        double var = .1;
        int count = 0;
        unsigned iter = 0;
        bool again = true;
        while(again)
        {
            epsilon = 1;
            iter = 0;
            count = 0;
            while(abs(epsilon) > .001 && count < 100)
            {
                again = false;
                double fxim1, fxi, fxip1 = 0;

                //To test the second equation,
                //we want to make an empty Polynomial
                if(coefficients.size() != 0)
                {
                    fxim1 = solveWithValue(xim1);
                    fxi = solveWithValue(xi);
                }
                else
                {
                    fxim1 = xim1 + 10 - xim1 * cosh(50/xim1);
                    fxi = xi + 10 - xi * cosh(50/xi);
                }
                xip1 = xi - (((xi - xim1) / (fxi - fxim1)) * fxi);
                //If the guess goes outside the boundary
//                if(xip1 > b || xip1 < a)
//                {
//                    cout << "Infinite loop?" << endl;
//                    again = true;
//                    cout << "xim1: " << xim1 << endl;
//                    cout << "xi: " << xi << endl;
//                    cout << "xip1: " << xip1 << endl;
//                    cout << "var: " << var << endl;
//                    if(xip1 <= start1)
//                    {
//                        xim1 = start1 + var;
//                        xi = start2 + var;
//                    }
//                    else if(xip1 >= start2)
//                    {
//                        xim1 = start1 - var;
//                        xi = start2 - var;
//                    }
//                    var +=.1;
//                    break;
//                }



                epsilon = abs(xip1 - xi)/xip1;

                cout << setprecision(0);
                cout << "|\t" << iter << "\t|" << setw(8) << setprecision(5)
                     << fixed << xim1 << "\t|\t" << setw(8) << setprecision(5)
                     << fixed << fxim1 << setw(8) << "\t|" << setw(8)
                     << setprecision(5) << fixed << xi << "\t|" << setw(8)
                     << setprecision(5) << fixed << fxi << "\t|" << setw(8)
                     << setprecision(5) << fixed << xip1 << "\t|" << setw(8)
                     << setprecision(5) << fixed << epsilon << "\t|"
                     << endl;

                if(coefficients.size() != 0)
                {
                    fxip1 = solveWithValue(xip1);
                }
                else
                {
                    fxip1 = xip1 + 10 - xip1 * cosh(50/xip1);
                }

                //Testing whether we actually got to the root we wanted
                if(abs(epsilon) < .001)
                {
                    //If the error is small but f(xip1) is not zero, we
                    //are oscillating around a point of zero slope
                    if(fxip1 > .0001)
                    {
                        again = true;
                        cout << "Oscillating detected. Going again.\n";

                        //Finding out which way to go to find another start
                        //point

                        //If the position is less than where we started, we
                        //want to go to the right
                        if(xip1 <= start1)
                        {
                            xim1 = start1 + var;
                            xi = start2 + var;
                        }
                        else if(xip1 >= start2)
                        {
                            xim1 = start1 - var;
                            xi = start2 - var;
                        }
                        else
                        {
                            xim1 = start1 - var;
                            xi = start2 + var;
                        }
                        var += .1;
                        break;
                    }
                }
                //cout << "e" << endl;
                xim1 = xi;
                xi = xip1;

                count++;
                iter++;
            }
        }
        cout << "Root is at: " << xip1 << endl;
        roots.push_back(xip1);

        cout << endl;
    }
    return roots;
}

vector<double> Polynomial::ModifiedSecant(double x0, double x1)
{
    float delta = static_cast<float>(.01);
    vector< pair<float, float> > boundaries;
    if(x0 == x1 && x0 == 0)
    {
        //Choosing the boundaries
        boundaries = getPossibleRoots();
    }
    else
    {
        boundaries.push_back({x0,x1});
    }
    vector<double> roots;

    //The number of roots discovered is the size of boundaries
    for(unsigned u = 0; u < boundaries.size(); u++)
    {
        pair<float, float> curr = boundaries.at(u);
        double epsilon = 1;
        cout << "|\tn\t" << "|\tdelta\t" << "|\tx_i\t\t" << "|\tf(x_i)\t"
             << "|\tx_i+1\t" << "|\te\t|" << endl;
        double a = static_cast<double>(curr.first);

        //Initial guess
        double xi = a;
        double xip1 = 0;

        int count = 0;
        unsigned iter = 0;
        bool again = true;
        while(again)
        {
            epsilon = 1;
            iter = 0;
            count = 0;
            while(abs(epsilon) > .001 && count < 100)
            {
                again = false;
                double xd = delta * xi;
                double fxi, fxd = 0;

                //To test the second equation,
                //we want to make an empty Polynomial
                if(coefficients.size() != 0)
                {
                    fxi = solveWithValue(xi);
                    fxd = solveWithValue((xi + xd));

                }
                else
                {
                    fxi = xi + 10 - xi * cosh(50/xi);
                    fxd = (xi + xd) + 10 - (xi + xd) * cosh(50/(xi + xd));
                }
                xip1 = xi - (fxi * (xd/(fxd - fxi)));

                epsilon = abs(xip1 - xi)/xip1;

                cout << setprecision(0);
                cout << "|\t" << iter << "\t|" << setw(8) << setprecision(2)
                     << fixed << delta << "\t|\t" << setw(8) << setprecision(5)
                     << fixed << xi << setw(8) << "\t|" << setw(8)
                     << setprecision(5) << fixed << fxi << "\t|" << setw(8)
                     << setprecision(5) << fixed << xip1 << "\t|" << setw(8)
                     << setprecision(5) << fixed << epsilon << "\t|" << endl;

                xi = xip1;
                count++;
                iter++;
            }
        }
        cout << "Root is at: " << xip1 << endl;
        roots.push_back(xip1);

        cout << endl;
    }
    return roots;

}
double Polynomial::getNormOfVector(vector<double> v)
{
    double total = 0;
    for(unsigned u = 0; u < v.size(); u++)
    {
        total += pow(v.at(u), 2);
    }
    return sqrt(total);
}

vector< pair<float, float> > Polynomial::
getBoundaries(short int stop, float step)
{
    pair<float, float> brackets;
    vector< pair<float, float> > boundaries;
    for(float i = -stop; i < stop; i+=step)
    {
        if(coefficients.size() != 0)
        {
            //Close enough to zero for our purposes
            if(abs(solveWithValue(static_cast<double>(i))) < 0.00001)
            {
                brackets.first = brackets.second = i;
                boundaries.push_back(brackets);
            }
            else if(solveWithValue(static_cast<double>(i)) *
                    solveWithValue(static_cast<double>(i + step)) < 0)
            {
                brackets.first = i;
                brackets.second = i + step;
                boundaries.push_back(brackets);
            }
        }
        else
        {
            float value = i + 10 - i * cosh(50/i);
            //cout << "Value is: " << value << endl;

            //Close enough to zero for our purposes
            if(abs(static_cast<double>(value)) < .00001)
            {
                brackets.first = brackets.second = i;
                boundaries.push_back(brackets);
            }
            else
            {
                float value2 = (i + step) + 10 - (i + step) *
                        cosh(50/(i + step));
                if(value * value2 < 0)
                {
                    brackets.first = i;
                    brackets.second = i + step;
                    boundaries.push_back(brackets);
                }
            }
        }

    }
    return boundaries;
}

double Polynomial::findMidpoint(double beg, double end)
{
    return (beg + end)/2;
}

void Polynomial::printTable(vector<double> entries)
{
    cout << "|\t" << entries.at(0);
    cout << "\t|\t" << setprecision(4) << fixed
         << entries.at(1) << "\t|\t" << setprecision(4)
         << fixed << entries.at(2);

        cout << "\t|\t" << setprecision(4) << fixed << entries.at(3);
    cout << "\t|\t" << setprecision(4) << fixed << entries.at(4) << "\t|\t"
         << setprecision(4) << fixed << entries.at(5) << "\t|\t"
         << setprecision(4) << fixed << entries.at(6);
    cout << "\t|\t" << setprecision(4) << fixed << entries.at(7)
         << "\t|" << endl;
}

vector< pair<float, float> > Polynomial::getPossibleRoots()
{
    vector< pair<float, float> > boundaries;
    //If we're testing an actual polynomial and not the second option
    if(coefficients.size() != 0)
    {
        short int n = 10;
        unsigned long long degree = coefficients.size();
        if(extra)
            degree--;
        while(n <= 1000)
        {
            boundaries = getBoundaries(n, 1);
            //The amount of roots a polynomial has won't be greater
            //than its degree
            if(boundaries.size() == degree)
                break;
            if(n < 100)
                n += 10;
            else
                n += 100;
        }

        //If there are fewer roots than the degree of the polynomial, it's
        //possible that there were more that one root in between the values
        //tested, so we test at -10, -9.5, -9, etc... instead of -10, -9, -8
        if(boundaries.size() < degree)
        {
            n = 10;
            while(n <= 1000)
            {
                boundaries = getBoundaries(n, .5);
                if(boundaries.size() == degree)
                    break;
                if(n < 100)
                    n += 10;
                else
                    n += 100;
            }
        }
    }
    //If we're testing the second function, we just find all roots between
    //-1000 and 1000 with a step size of .5
    else
    {
        boundaries = getBoundaries(1000, 1);
    }
    return boundaries;
}
