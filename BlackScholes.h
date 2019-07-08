#pragma once

#include "gaussians.h"

/*
inline double blackScholes(
const double spot, 
const double rate, 
const double yield, 
const double vol, 
const double strike, 
const double mat)
{
    double df = exp(-rate * mat), 
        fwd = spot * exp((rate - yield) * mat), 
        std = vol * sqrt(mat);
    double d = log(fwd / strike) / std;
    double d1 = d + 0.5 * std, d2 = d - 0.5 * std;
    double p1 = normalCdf(d1), p2 = normalCdf(d2);
    return df * (fwd * p1 - strike * p2);
}
*/

template <class T>
inline T blackScholes(
const T spot, 
const T rate, 
const T yield, 
const T vol, 
const T strike, 
const T mat)
{
    auto df = exp(-rate * mat), 
        fwd = spot * exp((rate - yield) * mat), 
        std = vol * sqrt(mat);
    auto d = log(fwd / strike) / std;
    auto d1 = d + 0.5 * std, d2 = d - 0.5 * std;
    auto p1 = normalCdf(d1), p2 = normalCdf(d2);
    return df * (fwd * p1 - strike * p2);
}