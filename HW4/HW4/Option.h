#ifndef OPTION_H
#define OPTION_H

#include <math.h>
#include <string>
#include "Normals.h"
#include "Random2.h"

class Option {
public:
	double spot, strike, interest, expiry, volatility, yield, barrier;
	double closedEuropeanCall();
	double monteCarloGenerator(double numPaths, RandomBase & generator);
};
#endif