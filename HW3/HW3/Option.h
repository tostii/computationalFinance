#ifndef OPTION_H
#define OPTION_H

#include <math.h>
#include "Normals.h"
#include "Random2.h"
#include <string>
class Option {
public:
	double spot, strike, interest, expiry, volatility, yield, barrier;
	double closedEuropeanCall();
	double monteCarloGenerator(double numPaths, RandomBase& generator);
	double monteCarloGaussian(double numPaths, std::string type);
};
#endif
