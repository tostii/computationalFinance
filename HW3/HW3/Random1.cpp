#include "Random1.h"
#include <cstdlib>
#include <cmath>

#include <iostream>
// the basic math functions should be in namespace
// std but aren't in VCPP6

#if !defined(_MSC_VER)
using namespace std;
#endif
ParkMiller pmGenerator;
LEcuyer leGenerator(1,2,3,5,7,11);

double GetOneGaussianBySummation()
{
	double result = 0;
	for (unsigned long j = 0; j < 12; j++)
		result += rand() / static_cast<double>(RAND_MAX);

	result -= 6.0;
	return result;
}


double GetOneGaussianByBoxMullerPark()
{
	double result;

	double x;
	double y;

	double sizeSquared;
	do {
		x = 2.0*pmGenerator.GetOneRandomInteger() / 2147483647 - 1;
		y = 2.0*pmGenerator.GetOneRandomInteger() / 2147483647 - 1;
		sizeSquared = x*x + y*y;
	} while
		(sizeSquared >= 1.0);

	result = x*sqrt(-2 * log(sizeSquared) / sizeSquared);

	return result;
}

double GetOneGaussianByInverseNormal() {
	double x = rand() / static_cast<double>(RAND_MAX);
	return InverseCumulativeNormal(x);
}

double GetOneGaussianByFishman() {
	double U1, U2, U3;
	double x;
	do {
		U1 = pmGenerator.GetOneRandomInteger() / static_cast<double>(2147483647);
		U2 = pmGenerator.GetOneRandomInteger() / static_cast<double>(2147483647);
		U3 = pmGenerator.GetOneRandomInteger() / static_cast<double>(2147483647);
		x = -log(U1);
	} while (U2 > exp(-0.5*((x - 1)*(x - 1))));
	return (U3 <= 0.5 ? -x : x);
}