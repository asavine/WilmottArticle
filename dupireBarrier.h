#pragma once

#include "matrix.h"
#include "random.h"
#include "interp.h"
#include "AAD.h"

#include <numeric>

template <class T>
inline T dupireBarrierMCBatch(
    //  Spot
    const T            S0,
    //  Local volatility
    const vector<T>&   spots,
    const vector<T>&   times,
    const matrix<T>&   vols,
    //  Product parameters
    const T            maturity,
    const T            strike,
    const T            barrier,
    //  First and last path
    const int          firstPath,
	const int		   lastPath,
	//	Time steps
    const int          Nt,
    //  Smoothing
    const T            epsilon,
    //  Random number generator
    RNG&               random)
{
    //  Initialize
    T result = 0;
    vector<double> gaussianIncrements(Nt); 
	//	double because the RNG is not templated
	//  	(and correctly so, see chapter 12)
	
	//	Set RNG state to the first path in the batch
	random.skipTo(firstPath);

    //  Loop over paths
    const T dt = maturity / Nt, sdt = sqrt(dt);
    for (int i = firstPath; i < lastPath; ++i)
    {
        //  Generate Nt Gaussian Numbers
        random.nextG(gaussianIncrements);

		//	Inntialize path
        T spot = S0, time = 0;
        T notionalAlive = 1.0; 
        
		//  Step by step
		for (int j = 0; j < Nt; ++j)
        {
            //  Interpolate volatility
            const T vol = interp2D(spots, times, vols, spot, time);
            //  Simulate return
            spot *= exp(-0.5 * vol * vol * dt + vol * sdt * gaussianIncrements[j]);
			//	Increase time
			time += dt;

            //  Monitor barrier
            if (spot > barrier + epsilon) { notionalAlive = 0.0; break; }       //   definitely dead
            else if (spot < barrier - epsilon) { /* do nothing */ }     //   definitely alive
            else /* in between, interpolate */ notionalAlive *= 1.0 - (spot - barrier + epsilon) / (2 * epsilon);
        }

        //  Payoff
        if (spot > strike) result += notionalAlive * (spot - strike); // pay on surviving notional
    }   

    return result / (lastPath - firstPath);
}

inline double dupireBarrierPricer(
    //  Spot
    const double			S0,
    //  Local volatility
    const vector<double>&   spots,
    const vector<double>&   times,
    const matrix<double>&   vols,
    //  Product parameters
    const double            maturity,
    const double            strike,
    const double            barrier,
    //  Number of paths
    const int				Np,
	//	Number of simulations in every batch
	const int				Nb,
	//	Time steps
    const int				Nt,
    //  Smoothing
    const double            epsilon,
    //  Random number generator
    RNG&					random)
{	
	random.init(Nt);
	double result = 0.0;
	int firstPath = 0;

	while (firstPath < Np)
	{
		int lastPath = firstPath + Nb;
		lastPath = min(lastPath, Np);
		double batchPrice = dupireBarrierMCBatch(
			S0, 
			spots, 
			times, 
			vols, 
			maturity, 
			strike, 
			barrier, 
			firstPath, 
			lastPath, 
			Nt, 
			epsilon, 
			random);
		result += batchPrice * (lastPath-firstPath) / Np;
		firstPath = lastPath;
	}

	return result;
}

inline void dupireBarrierRisks(
    //  Spot
    const double			S0,
    //  Local volatility
    const vector<double>&   spots,
    const vector<double>&   times,
    const matrix<double>&   vols,
    //  Product parameters
    const double            maturity,
    const double            strike,
    const double            barrier,
    //  Number of paths
    const int				Np,
	//	Number of simulations in every batch
	const int				Nb,
	//	Time steps
    const int				Nt,
    //  Smoothing
    const double            epsilon,
    //  Random number generator
    RNG&					random,
	//	Results
	double&					price,
	double&					delta,
	matrix<double>&			vegas)
{	
	//	Allocate and initialize results
	price = 0.0;
	delta = 0.0;
	vegas.resize(spots.size(), times.size());
	for (auto& vega : vegas) vega = 0.0;

	//	Temporary storage for batch risks
	double batchPrice, batchDelta;
	matrix<double> batchVegas(vegas.rows(), vegas.cols());

	//	Temporary storage for parameters as Number types
	Number nS0, nMaturity, nStrike, nBarrier, nEpsilon;
	vector<Number> nSpots(spots.size()), nTimes(times.size());
	matrix<Number> nVols(vols.rows(), vols.cols());

	//	Initialize the RNG
	random.init(Nt);

	//	Loop over batches
	int firstPath = 0;
	while (firstPath < Np)
	{
		int lastPath = firstPath + Nb;
		lastPath = min(lastPath, Np);

		//	Wipe the tape
		tape.clear();
		
		//	Put parameters on tape by initialization of Number types
		//	Note this is inefficient in many ways but it keeps things simple 
		//	and doesn't matter too much as long as the number of batches remains low
		nS0 = S0; 
		nMaturity = maturity;
		nStrike = strike; 
		nBarrier = barrier;
		nEpsilon = epsilon;
		copy(spots.begin(), spots.end(), nSpots.begin());
		copy(times.begin(), times.end(), nTimes.begin());
		copy(vols.begin(), vols.end(), nVols.begin());
	
		//	Compute the batch
		Number nBatchPrice = dupireBarrierMCBatch(
			nS0, 
			nSpots, 
			nTimes, 
			nVols, 
			nMaturity, 
			nStrike, 
			nBarrier, 
			firstPath, 
			lastPath, 
			Nt, 
			nEpsilon, 
			random);

		//	Back-propagate derivatives
		vector<double> adjoints = calculateAdjoints(nBatchPrice);

		//	Pick results
		batchPrice = nBatchPrice.value;
		batchDelta = adjoints[nS0.idx];
		transform(nVols.begin(), nVols.end(), batchVegas.begin(),
			[&](const Number& vol) { return adjoints[vol.idx]; });

		//	Accumulate
		int paths = lastPath - firstPath;
		double w = double(paths) / Np;
		price += batchPrice * w;
		delta += batchDelta * w;
		transform(vegas.begin(), vegas.end(), batchVegas.begin(), vegas.begin(),
			[&](const double& vega, const double& batchVega) { return vega + batchVega * w; });

		//	Next batch
		firstPath = lastPath;
	}
}

inline double dupireBarrierPricerMT(
    //  Spot
    const double			S0,
    //  Local volatility
    const vector<double>&   spots,
    const vector<double>&   times,
    const matrix<double>&   vols,
    //  Product parameters
    const double            maturity,
    const double            strike,
    const double            barrier,
    //  Number of paths
    const int				Np,
	//	Number of simulations in every batch
	const int				Nb,
	//	Time steps
    const int				Nt,
    //  Smoothing
    const double            epsilon,
    //  Random number generator
    RNG&					random)
{	
	//  Memory for the storage of batch-wise results
    const int numBatches = int((Np - 1) / Nb) + 1;
    vector<double> batchResults(numBatches);

	//	Initialize the RNG
	random.init(Nt);

	//	Iterate over batches, in parallel
	#pragma omp parallel for
	for(int batch=0; batch<numBatches; ++batch)
	{
		const int firstPath = batch * Nb;
		const int lastPath = min(firstPath + Nb, Np);

        //  Make a copy of the (mutable) RNG
        auto cRandom = random.clone();

        //  Process the batch
        batchResults[batch] = (lastPath - firstPath) 
            * dupireBarrierMCBatch(
                S0,
                spots,
                times,
                vols,
                maturity,
                strike,
                barrier,
                firstPath,
                lastPath,
                Nt,
                epsilon,
                *cRandom);   //  call with own copy of RNG
	}
    
    //  Average results over batches
    return accumulate(batchResults.begin(), batchResults.end(), 0.0) / Np;
}

inline void dupireBarrierRisksMT(
    //  Spot
    const double			S0,
    //  Local volatility
    const vector<double>&   spots,
    const vector<double>&   times,
    const matrix<double>&   vols,
    //  Product parameters
    const double            maturity,
    const double            strike,
    const double            barrier,
    //  Number of paths
    const int				Np,
	//	Number of simulations in every batch
	const int				Nb,
	//	Time steps
    const int				Nt,
    //  Smoothing
    const double            epsilon,
    //  Random number generator
    RNG&					random,
	//	Results
	double&					price,
	double&					delta,
	matrix<double>&			vegas)
{	
	//  Memory for the storage of batch-wise results
    int numBatches = int((Np - 1) / Nb) + 1;
    vector<double> batchPrices(numBatches);
    vector<double> batchDeltas(numBatches);
    vector<matrix<double>> batchVegas(numBatches);
    for (auto& batchVega : batchVegas)
    {
	    batchVega.resize(spots.size(), times.size());
    }

	//	Initialize the RNG
	random.init(Nt);
	 
	//	Iterate over batches, in parallel
	#pragma omp parallel for
	for(int batch=0; batch<numBatches; ++batch)
	{
		const int firstPath = batch * Nb;
		const int lastPath = min(firstPath + Nb, Np);

        //  Make a copy of the (mutable) RNG
        auto cRandom = random.clone();

        //	Wipe the tape
        tape.clear();

        //	Put parameters on tape by initialization of Number types
        
		//	Working memory
		static thread_local Number nS0, nMaturity, nStrike, nBarrier, nEpsilon;
        static thread_local vector<Number> nSpots, nTimes;
        static thread_local matrix<Number> nVols;

		//	Allocate
		nSpots.resize(spots.size()); 
		nTimes.resize(times.size());
		nVols.resize(vols.rows(), vols.cols());

		//	Initialize and put on tape
		nS0 = S0;
        nMaturity = maturity;
        nStrike = strike;
        nBarrier = barrier;
        nEpsilon = epsilon;
        copy(spots.begin(), spots.end(), nSpots.begin());
        copy(times.begin(), times.end(), nTimes.begin());
        copy(vols.begin(), vols.end(), nVols.begin());

        //	Process the batch
        Number nBatchPrice = dupireBarrierMCBatch(
            nS0,
            nSpots,
            nTimes,
            nVols,
            nMaturity,
            nStrike,
            nBarrier,
            firstPath,
            lastPath,
            Nt,
            nEpsilon,
            *cRandom);

        //	Back-propagate derivatives
        vector<double> adjoints = calculateAdjoints(nBatchPrice);

        //	Pick results
        int paths = lastPath - firstPath;
        batchPrices[batch] = nBatchPrice.value * paths;
        batchDeltas[batch] = adjoints[nS0.idx] * paths;
        transform(nVols.begin(), nVols.end(), batchVegas[batch].begin(),
            [&](const Number& vol) { return adjoints[vol.idx] * paths; });

	}

    //  Average results over batches
    price = accumulate(batchPrices.begin(), batchPrices.end(), 0.0) / Np;
    delta = accumulate(batchDeltas.begin(), batchDeltas.end(), 0.0) / Np;    
    vegas.resize(spots.size(), times.size());
    for (auto& vega : vegas) vega = 0.0;
	//	Average vegas in parallel
	#pragma omp parallel for
	for (int i = 0; i < vegas.rows(); ++i)
	{
		for (int j = 0; j < vegas.cols(); ++j)
		{
			for (int k = 0; k < numBatches; ++k)
			{
				vegas[i][j] += batchVegas[k][i][j];
			}
			vegas[i][j] /= Np;
		}
	}
}
