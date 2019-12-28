
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

//  Excel export wrappers to functions in main.h

#pragma warning(disable:4996)

#include "dupireBarrier.h"
#include "sobol.h"
#include "mrg32k3a.h"

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "xlcall.h"
#include "xlframework.h"
#include "xlOper.h"

//  Helpers
struct NumericalParam
{
    bool              useSobol;
    int               numPath;
    int               seed1 = 12345;
    int               seed2 = 1234;
};
NumericalParam xl2num(
    const double              useSobol,
    const double              seed1,
    const double              seed2,
    const double              numPath)
{
    NumericalParam num;

    num.numPath = static_cast<int>(numPath + EPS);
	if (seed1 >= 1)
	{
		num.seed1 = static_cast<int>(seed1 + EPS);
	}
	else
	{
		num.seed1 = 1234;
	}

	if (seed2 >= 1)
	{
		num.seed2 = static_cast<int>(seed2 + EPS);
	}
	else
	{
		num.seed2 = num.seed1 + 1;
	}

	num.useSobol = useSobol > EPS;

    return num;
}

//	Wrappers

extern "C" __declspec(dllexport)
double xDupireBarrierPricer(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,
    double              mat,
    double              strike,
    double              barrier,
    double              paths,
	double				batchSize,
    double              steps,
    double              epsilon,
    double              useSobol,
    double              seed1,
    double              seed2,
    double              parallel)
{
    FreeAllTempMemory();

    //  Make sure we have paths and steps
    if (paths <= 0.0 || steps <= 0.0) return -1;

    //  Unpack

    if (spots->rows * spots->columns * times->rows * times->columns != vols->rows * vols->columns)
    {
        return -1;
    }

    vector<double> vspots = to_vector(spots);
    vector<double> vtimes = to_vector(times);
    matrix<double> vvols = to_matrix(vols);

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (useSobol > 0.5) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1 > 0.5 ? int(seed1): 12345, seed2 > 0.5? int(seed2): 123456);
	rng->init(int(steps));

    //  Call and return
    return parallel > 0.5 
        ? dupireBarrierPricerMT(spot, vspots, vtimes, vvols, mat, strike, barrier, int(paths), int(batchSize), int(steps), spot*epsilon, *rng)
        : dupireBarrierPricer(spot, vspots, vtimes, vvols, mat, strike, barrier, int(paths), int(batchSize), int(steps), spot*epsilon, *rng);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xDupireBarrierRisks(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,
    double              mat,
    double              strike,
    double              barrier,
    double              paths,
    double              batchSize,
    double              steps,
    double              epsilon,
    double              useSobol,
    double              seed1,
    double              seed2,
    double              parallel)
{
    FreeAllTempMemory();

    //  Make sure we have paths and steps
    if (paths <= 0.0 || steps <= 0.0) TempErr12(xlerrNA);

    //  Unpack

    if (spots->rows * spots->columns * times->rows * times->columns != vols->rows * vols->columns)
    {
        return TempErr12(xlerrNA);
    }

    vector<double> vspots = to_vector(spots);
    vector<double> vtimes = to_vector(times);
    matrix<double> vvols = to_matrix(vols);

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (useSobol > 0.5) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1 > 0.5 ? int(seed1): 12345, seed2 > 0.5? int(seed2): 123456);
	rng->init(int(steps));

    //  Call 
	double price, delta;
	matrix<double> vegas(vvols.rows(), vvols.cols());
    
    if (parallel > 0.5)
    {
        dupireBarrierRisksMT(spot, vspots, vtimes, vvols, mat, strike, barrier, int(paths), int(batchSize), int(steps), spot*epsilon, *rng,
		    price, delta, vegas);        
    }
    else
    {
        dupireBarrierRisks(spot, vspots, vtimes, vvols, mat, strike, barrier, int(paths), int(batchSize), int(steps), spot*epsilon, *rng,
		    price, delta, vegas);    
    }

	//	Pack and return
	LPXLOPER12 results = TempXLOPER12();
	resize(results, 2 + vegas.rows()*vegas.cols(), 1);
	setNum(results, price, 0, 0);
	setNum(results, delta, 1, 0);
	size_t r = 2;
	for (int i=0; i<vegas.rows(); ++i) for (int j=0; j<vegas.cols(); ++j)
	{
		setNum(results, vegas[i][j], r, 0);
		++r;
	}

	return results;
}

//	Registers

extern "C" __declspec(dllexport) int xlAutoOpen(void)
{
	XLOPER12 xDLL;

	Excel12f(xlGetName, &xDLL, 0);

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireBarrierPricer"),
        (LPXLOPER12)TempStr12(L"BBK%K%K%BBBBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xDupireBarrierPricer"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, mat, strike, barrier, paths, batchSize, steps, epsilon, useSobol, [seed1], [seed2], [parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Toy Dupire Barrier MC"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireBarrierRisks"),
        (LPXLOPER12)TempStr12(L"QBK%K%K%BBBBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xDupireBarrierRisks"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, mat, strike, barrier, paths, batchSize, steps, epsilon, useSobol, [seed1], [seed2], [parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Toy Dupire Barrier MC AAD risks"),
        (LPXLOPER12)TempStr12(L""));

	/* Free the XLL filename */
	Excel12f(xlFree, 0, 1, (LPXLOPER12)&xDLL);

	return 1;
}

extern "C" __declspec(dllexport) int xlAutoClose(void)
{

    return 1;
}