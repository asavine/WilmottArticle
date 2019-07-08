#pragma once

#include <vector>
#include <iostream>
using namespace std;

#include "gaussians.h"
#include "BlackScholes.h"

struct Node
{
    int     numArg;       //  number of arguments: 0, 1 or 2
    int     idx1;         //  index of first argument on tape
    int     idx2;         //  index of second argument on tape
    double  der1;         //  partial derivative to first argument
    double  der2;         //  partial derivative to second argument
};

//  The tape, declared as a global variable
vector<Node> tape;

struct Number
{
    double  value;
    int     idx;

    //  default constructor does nothing
    Number() {}

    //  constructs with a value and record
    Number(const double& x) : value(x)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  reference record on tape
        idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 0;
    }

    Number operator +() const { return *this; }
    Number operator -() const { return Number(0.0) - *this; }

    Number& operator +=(const Number& rhs) { *this = *this + rhs; return *this; }
    Number& operator -=(const Number& rhs) { *this = *this - rhs; return *this; }
    Number& operator *=(const Number& rhs) { *this = *this * rhs; return *this; }
    Number& operator /=(const Number& rhs) { *this = *this / rhs; return *this; }

    friend Number operator+(const Number& lhs, const Number& rhs)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = lhs.value + rhs.value; //  calling double overload 

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 2;
        node.idx1 = lhs.idx;
        node.idx2 = rhs.idx;

        //  compute derivatives, both derivatives of addition are 1
        node.der1 = 1;
        node.der2 = 1;

        return result;
    }

    friend Number operator-(const Number& lhs, const Number& rhs)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = lhs.value - rhs.value; //  calling double overload 

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 2;
        node.idx1 = lhs.idx;
        node.idx2 = rhs.idx;

        //  compute derivatives, both derivatives of addition are 1
        node.der1 = 1;
        node.der2 = -1;

        return result;
    }

    friend Number operator*(const Number& lhs, const Number& rhs)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = lhs.value * rhs.value; //  Different value here

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 2;
        node.idx1 = lhs.idx;
        node.idx2 = rhs.idx;

        //  compute derivatives, both derivatives of addition are 1
        node.der1 = rhs.value;                  //  Different derivative here
        node.der2 = lhs.value;                  //  Different derivative here

        return result;
    }

    friend Number operator/(const Number& lhs, const Number& rhs)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = lhs.value / rhs.value; //  Different value here

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 2;
        node.idx1 = lhs.idx;
        node.idx2 = rhs.idx;

        //  compute derivatives, both derivatives of addition are 1
        node.der1 = 1.0 / rhs.value;;
        node.der2 = -lhs.value / (rhs.value * rhs.value);

        return result;
    }

    friend Number log(const Number& arg)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = log(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 1;
        node.idx1 = arg.idx;

        //  compute derivative
        node.der1 = 1.0 / arg.value;

        return result;
    }

    friend Number exp(const Number& arg)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = exp(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 1;
        node.idx1 = arg.idx;

        //  compute derivative
        node.der1 = result.value;

        return result;
    }

    friend Number sqrt(const Number& arg)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = sqrt(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 1;
        node.idx1 = arg.idx;

        //  compute derivative
        node.der1 = 0.5 / result.value;

        return result;
    }

    friend Number normalDens(const Number& arg)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = normalDens(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 1;
        node.idx1 = arg.idx;

        //  compute derivative
        node.der1 = -result.value * arg.value;

        return result;
    }

    friend Number normalCdf(const Number& arg)
    {
        //  create a new record on tape
        tape.push_back(Node());
        Node& node = tape.back();

        //  compute result
        Number result;
        result.value = normalCdf(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        node.numArg = 1;
        node.idx1 = arg.idx;

        //  compute derivative
        node.der1 = normalDens(arg.value);

        return result;
    }
};

vector<double> calculateAdjoints(Number& result)
{
    //  initialization
    vector<double> adjoints(tape.size(), 0.0);  //  initialize all to 0
    int N = result.idx;                         //  find N
    adjoints[N] = 1.0;                          //  seed aN = 1
    
    //  backward propagation
    for(int j=N; j>0; --j)  //  iterate backwards over tape
    {
        if (tape[j].numArg > 0)
        {
            //  propagate first argument
            adjoints[tape[j].idx1] += adjoints[j] * tape[j].der1;       
            if (tape[j].numArg > 1)
            {
                //  propagate second argument
                adjoints[tape[j].idx2] += adjoints[j] * tape[j].der2;   
            }
        }
    }

    return adjoints;
}

void differentiateBlackScholes()
{
    // initializes and records inputs
    Number spot = 100, 
        rate = 0.02, 
        yield = 0.05, 
        vol = 0.2, 
        strike = 110, 
        mat = 2; 
    // evaluates and records operations
    auto result = blackScholes(spot, rate, yield, vol, strike, mat);                
    cout << "Value = " << result.value << endl;   //  5.03705
 
    //  propagate adjoints
    vector<double> adjoints = calculateAdjoints(result);

    //  show derivatives
    cout << "Derivative to spot (delta) = " 
        << adjoints[spot.idx] << endl;          
    //  0.309
    cout << "Derivative to rate (rho) = " 
        << adjoints[rate.idx] << endl;
    //  51.772
    cout << "Derivative to dividend yield = " 
        << adjoints[yield.idx] << endl;       
    //  -61.846
    cout << "Derivative to volatility (vega) = " 
        << adjoints[vol.idx] << endl;      
    //  46.980
    cout << "Derivative to strike (-digital) = " 
        << adjoints[strike.idx] << endl;   
    //  -0.235
    cout << "Derivative to maturity (-theta) = " 
        << adjoints[mat.idx] << endl;      
    //  1.321

    //  clear
    tape.clear();
}