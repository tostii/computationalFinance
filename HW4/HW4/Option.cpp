#include "Option.h"
double Option::closedEuropeanCall()
{
	double d1 = (log(spot / strike) + (interest - yield + pow(volatility, 2) / 2)*expiry) / (volatility * sqrt(expiry));
	double d2 = d1 - volatility*sqrt(expiry);
	return spot*exp(-yield*expiry)*CumulativeNormal(d1) - strike*exp(-interest*expiry)*CumulativeNormal(d2);
}


double Option::monteCarloGenerator(double numPaths, RandomBase& generator)
{
	generator.ResetDimensionality(1);
	double variance = volatility * volatility * expiry;
	double rootVariance = sqrt(variance);
	double itoCorrection = -0.5 * variance;
	double movedSpot = spot * exp(interest * expiry + itoCorrection);
	double thisSpot;
	double runningSum = 0;
	MJArray VariateArray(2);
	for (unsigned long i = 0; i < numPaths; i++) {
		generator.GetGaussians(VariateArray);
		thisSpot = movedSpot*exp(rootVariance*VariateArray[0]);
		double thisPayoff = thisSpot - strike;
		thisPayoff = thisPayoff > 0 ? thisPayoff : 0;
		runningSum += thisPayoff;
	}
	double mean = runningSum / numPaths;
	mean *= exp(-interest*expiry);
	return mean;
}