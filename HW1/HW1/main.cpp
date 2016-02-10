// requires Random1.cpp

#include "Random1.h"
#include "Normals.h"
#include <iostream>
#include <cmath>
#include <string>
#include <time.h>


//generic monte carlo simulation that takes the method used as a parameter
double optionMC(double spot, double strike, double r, double expiry, double vol, double yield, double numPath, double numStep, std::string type, double (*f)(double, double, double, double, double, double, double)){	
	double runningSum = 0;
	for (unsigned long j = 0; j < numPath; j++) {
		double iterateSpot = spot;
		for (int i = 1; i <= numStep; i++) {
			iterateSpot = f(iterateSpot, strike, r, expiry, vol, yield, numStep);
		}
		double thisPayoff;
		if (type == "put") {
			thisPayoff = strike - iterateSpot;
		}
		else if (type == "call") {
			thisPayoff = iterateSpot - strike;
		}
		else {
			return -1;
		}
		thisPayoff = thisPayoff > 0 ? thisPayoff : 0;
		runningSum += thisPayoff;
	}
	double mean = runningSum / numPath;

	mean *= exp(-r*expiry);
	return mean;
}

//closed-form
double closedPut(double spot, double strike, double r, double expiry, double vol, double yield, std::string type) {
	double d1 = (log(spot / strike) + (r - yield + pow(vol, 2) / 2)*expiry) / (vol * sqrt(expiry));
	double d2 = d1 - vol*sqrt(expiry);
	if (type == "put") {
		return strike*exp(-r*expiry)*(1 - CumulativeNormal(d2)) - spot*exp(-yield*expiry)*(1 - CumulativeNormal(d1));
	}
	else if (type == "call") {
		return spot*exp(-yield*expiry)*CumulativeNormal(d1) - strike*exp(-r*expiry)*CumulativeNormal(d2);
	}
	else {
		return -1;
	}
}

//analytical one step
double analytical(double spot, double strike, double r, double expiry, double vol, double yield, double numStep) {
	return spot*exp(vol*expiry*GetOneGaussianByBoxMuller());
}

//euler
double euler(double spot, double strike, double r, double expiry, double vol, double yield, double numStep) {
	double dt = 1 / numStep;
	double dW = sqrt(dt) * GetOneGaussianByBoxMuller();
	double dS = spot*((r - yield)*dt + vol*dW);
	return spot + dS;
}

//euler log
double eulerLog(double spot, double strike, double r, double expiry, double vol, double yield, double numStep) {
	static int count = 1;
	double dt = 1 / numStep;
	double dW = sqrt(dt) * GetOneGaussianByBoxMuller();
	double dS = ((r - yield - pow(vol, 2) / 2)*dt + vol*dW);
	if (count == numStep) {
		count = 1;
		return exp(spot + dS);

	}
	else {
		count++;
		return spot + dS;
	}
}

//milstein
double milstein(double spot, double strike, double r, double expiry, double vol, double yield, double numStep) {
	double dt = 1 / numStep;
	double dW = sqrt(dt) * GetOneGaussianByBoxMuller();
	double dS = spot*((r - yield)*dt + vol*dW + 0.5*vol*vol*(dW*dW - dt));
	return spot + dS;
}


//main function
int main() {
	//begin input
	const double spot = 100;
	const double strike = 110;
	const double r = 0.05;
	const double vol = 0.3;
	const double yield = 0.02;
	const double expiry = 1;
	const double numStep = 252;
	const double numPath = 10000;
	srand(time(NULL));
	//end input

	std::cout << "Assuming that the stock price is geometric Brownian motion with volatility " << vol << ", initial asset price S(0) = " << spot << ", constant risk-free interest rate r = " << r << ", dividend yield = " << yield << ", strike price K = " << strike << ", and maturity T = " << expiry << std::endl;
	std::cout << std::endl;

	std::cout << "European Call Option" << std::endl;
	std::cout << "Option price using closed-form formula: " << closedPut(spot, strike, r, expiry, vol, yield, "call") << std::endl;
	std::cout << "Option price using single-step exact SDE solution: " << optionMC(spot*exp(r*expiry - 0.5*vol*vol*expiry - yield), strike, r, expiry, vol, yield, numPath, 1, "call", (*analytical)) << std::endl;
	std::cout << "Option price using Euler numerical solution of SDE for spot: " << optionMC(spot, strike, r, expiry, vol, yield, numPath, numStep, "call", (*euler)) << std::endl;
	std::cout << "Option price using Euler numerical solutionn of SDE for log spot: " << optionMC(log(spot), strike, r, expiry, vol, yield, numPath, numStep, "call", (*eulerLog)) << std::endl;
	std::cout << "Option price using Milstein numerical solution of SDE for spot: " << optionMC(spot, strike, r, expiry, vol, yield, numPath, numStep, "call", (*milstein)) << std::endl;
	std::cout << std::endl;

	std::cout << "European Put Option" << std::endl;
	std::cout << "Option price using closed-form formula: " << closedPut(spot, strike, r, expiry, vol, yield, "put") << std::endl;
	std::cout << "Option price using single-step exact SDE solution: " << optionMC(spot*exp(r*expiry - 0.5*vol*vol*expiry - yield), strike, r, expiry, vol, yield, numPath, 1, "put", (*analytical)) << std::endl;
	std::cout << "Option price using Euler numerical solution of SDE for spot: " << optionMC(spot, strike, r, expiry, vol, yield, numPath, numStep, "put", (*euler)) << std::endl;
	std::cout << "Option price using Euler numerical solutionn of SDE for log spot: " << optionMC(log(spot), strike, r, expiry, vol, yield, numPath, numStep, "put", (*eulerLog)) << std::endl;
	std::cout << "Option price using Milstein numerical solution of SDE for spot: " << optionMC(spot, strike, r, expiry, vol, yield, numPath, numStep, "put", (*milstein)) << std::endl;
	std::cin.get();
	return 0;

}
