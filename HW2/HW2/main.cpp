#include <iostream>
#include <time.h>
#include "Normals.h"
#include "Ziggurat.h"
#include "Zigg.h"
#include "Random1.h"
#include <vector>
#include <numeric>
#include <algorithm>
#include "runningStatistics.h"
#include "runningRegression.h"


struct Option {
	double spot;
	double strike;
	double interest;
	double expiry;
	double volatility;
	double yield;
	double barrier;
};

double euler(Option o, double spot, double numStep, double norm) {
	double dt = 1 / numStep;
	double dW = sqrt(dt) * norm;
	double dS = spot*((o.interest - o.yield)*dt + o.volatility*dW);
	return spot + dS;
}


double closedVanilla(Option o) {
	double d1 = (log(o.spot / o.strike) + (o.interest - o.yield + pow(o.volatility, 2) / 2) *o.expiry) / (o.volatility * sqrt(o.expiry));
	double d2 = d1 - o.volatility*sqrt(o.expiry);
	return o.spot*exp(-o.yield*o.expiry)*CumulativeNormal(d1) - o.strike*exp(-o.interest*o.expiry)*CumulativeNormal(d2);
}

std::vector<double> generatePair(Option o, double numStep, double(*f)(Option o, double, double, double)) {
	std::vector<double> pair;
	double barrierSpot = o.spot;
	double vanillaSpot = o.spot;
	double barrierPayoff = 0;
	double vanillaPayoff = 0;
	for (int i = 1; i <= numStep; i++){
		double myNorm = GetOneGaussianByBoxMuller();
		barrierSpot = f(o, barrierSpot, numStep, myNorm);
		vanillaSpot = f(o, vanillaSpot, numStep, myNorm);
		if (barrierSpot > o.barrier) {
			break;
		}
	}
	if (barrierSpot < o.barrier) {
		barrierPayoff = barrierSpot - o.strike;
		barrierPayoff = barrierPayoff > 0 ? barrierPayoff : 0;
		vanillaPayoff = barrierPayoff;
	}
	if (barrierSpot > o.barrier) {
		vanillaPayoff = vanillaSpot - o.strike;
		vanillaPayoff = vanillaPayoff > 0 ? vanillaPayoff : 0;
		barrierPayoff = 0;
	}
	pair.push_back(vanillaPayoff);
	pair.push_back(barrierPayoff);
	return pair;
}

double controlMC(Option o, double numPath, double numStep, double numPilotSim, double(*f)(Option o, double, double, double)) {
	runningRegression lineReg;
	runningStatistics statistics;
	runningStatistics statisticsx;
	runningStatistics statisticsy;
	time_t time1 = clock();
	for (int i = 1; i <= numPilotSim; i++) {
		std::vector<double> pair = generatePair(o, numStep, (*euler));
		lineReg.Push(pair.at(0), pair.at(1));
	}
	double b_optimal = lineReg.Slope();
	for (int j = 1; j < numPath; j++) {
		std::vector<double> pair1 = generatePair(o, numStep, (*euler));
		double thisPayoff = pair1.at(1)- (b_optimal * (pair1.at(0) - closedVanilla(o)));
		statisticsx.Push(pair1.at(0));
		statisticsy.Push(pair1.at(1));
		statistics.Push(thisPayoff);
	}
	time_t time2 = clock();
	double time_elapsed = double((time2 - time1) / static_cast<double>(CLOCKS_PER_SEC));
	double controlMC = statistics.Mean() * exp(-o.interest * o.expiry);
	std::cout << controlMC;
	std::cout << ", std error = " << sqrt(statistics.Variance()) / sqrt(numPath);
	std::cout << ", run time = " << time_elapsed;
	std::cout << ", control coefficient = " << b_optimal;
	std::cout << ", correlation coefficient = " << b_optimal / sqrt(statisticsy.Variance()/statisticsx.Variance()) << std::endl;
	return statistics.Mean();
}

double optionMC(Option o, double numPath, double numStep, double(*f)(Option o, double, double, double)) {
	runningStatistics statistics;
	double dt = 1 / numStep;
	time_t time1 = clock();

	for (unsigned long j = 0; j < numPath; j++) {
		Ziggurat a;
		a.setSeed(j + 1);
		double iterateSpot = o.spot;
		for (int i = 1; i <= numStep; i++) {
			iterateSpot = f(o, iterateSpot, numStep, a.norm());
			if (iterateSpot > o.barrier){
				statistics.Push(0);
				break;
			}
		}
		double thisPayoff;
		if (iterateSpot < o.barrier) {
			thisPayoff = iterateSpot - o.strike;
			thisPayoff = thisPayoff > 0 ? thisPayoff : 0;
			statistics.Push(thisPayoff);
		}
	}
	
	time_t time2 = clock();
	double time_elapsed = double((time2 - time1) / static_cast<double>(CLOCKS_PER_SEC));
	double eulerMC = statistics.Mean() * exp(-o.interest * o.expiry);
	std::cout << eulerMC;
	std::cout << ", std error = " << sqrt(statistics.Variance()/numPath);
	std::cout << ", run time = " << time_elapsed << std::endl;
	return eulerMC;
}

double antitheticMC(Option o, double numPath, double numStep, double(*f)(Option o, double, double, double)) {
	runningStatistics statistics;
	runningStatistics statistics2;
	time_t time1 = clock();
	std::vector<double> price;
	std::vector<double> aPrice;
	for (unsigned long j = 0; j < numPath; j++) {
		bool calcZ = true;
		bool calcNZ = true;
		Ziggurat a;
		a.setSeed(j + 1);
		double iterateSpot = o.spot;
		double antitheticSpot = o.spot;

		for (int i = 1; i <= numStep; i++) {
			double Z = a.norm();
			if (calcZ == true) {
				iterateSpot = f(o, iterateSpot, numStep, Z);
			}
			if (calcNZ == true) {
				antitheticSpot = f(o, antitheticSpot, numStep, -Z);
			}

			if (iterateSpot > o.barrier && calcZ == true) {
				statistics.Push(0);
				calcZ = false;
			}
			if (antitheticSpot > o.barrier && calcNZ == true) {
				statistics2.Push(0);
				calcNZ = false;
			}
			if (calcZ == false && calcNZ == false) {
				break;
			}
		}
		double thisPayoff;
		if (iterateSpot < o.barrier) {
			thisPayoff = iterateSpot - o.strike;
			thisPayoff = thisPayoff > 0 ? thisPayoff : 0;
			statistics.Push(thisPayoff);
		}
		if (antitheticSpot < o.barrier) {
			thisPayoff = antitheticSpot - o.strike;
			thisPayoff = thisPayoff > 0 ? thisPayoff : 0;
			statistics2.Push(thisPayoff);
		}
	}
	time_t time2 = clock();
	double time_elapsed = double((time2 - time1) / static_cast<double>(CLOCKS_PER_SEC));
	double eulerMC = 0.5*(statistics.Mean() + statistics2.Mean()) * exp(-o.interest * o.expiry);
	std::cout << eulerMC;
	std::cout << ", std error = " << sqrt(statistics.Variance() + statistics2.Variance())/sqrt(2*numPath);
	std::cout << ", run time = " << time_elapsed << std::endl;
	return eulerMC;
}

double analytical(Option o) {
	time_t time1 = clock();
	double d1 = (log(o.spot / o.strike) + (o.interest - o.yield + 0.5 * o.volatility * o.volatility)*o.expiry) / (o.volatility*sqrt(o.expiry));
	double d2 = d1 - o.volatility * sqrt(o.expiry);
	double c = o.spot * exp(-o.yield * o.expiry)*CumulativeNormal(d1) - o.strike * exp(-o.interest * o.expiry)*CumulativeNormal(d2);
	double lambda = (o.interest - o.yield + (o.volatility*o.volatility / 2)) / (o.volatility * o.volatility);
	double y = log((o.barrier*o.barrier) / (o.spot*o.strike)) / (o.volatility * sqrt(o.expiry)) + lambda * o.volatility * sqrt(o.expiry);
	double x1 = log(o.spot / o.barrier) / (o.volatility * sqrt(o.expiry)) + lambda * o.volatility * sqrt(o.expiry);
	double y1 = log(o.barrier / o.spot) / (o.volatility * sqrt(o.expiry)) + lambda * o.volatility * sqrt(o.expiry);
	double t1 = o.spot * CumulativeNormal(x1)*exp(-o.yield*o.expiry);
	double t2 = -o.strike * exp(-o.interest * o.expiry) * CumulativeNormal(x1 - o.volatility * sqrt(o.expiry));
	double t3 = -o.spot * exp(-o.yield * o.expiry)*pow(o.barrier / o.spot, 2 * lambda)*(CumulativeNormal(-y) - CumulativeNormal(-y1));
	double t4 = o.strike * exp(-o.interest * o.expiry) * pow(o.barrier / o.spot, 2 * lambda - 2)*(CumulativeNormal(-y + o.volatility*sqrt(o.expiry)) - CumulativeNormal(-y1 + o.volatility*sqrt(o.expiry)));
	time_t time2 = clock();
	double time_elapsed = double((time2 - time1) / static_cast<double>(CLOCKS_PER_SEC));
	double closedPrice = c - (t1 + t2 + t3 + t4);
	std::cout << closedPrice << ", run time = " << time_elapsed << std::endl;
	return closedPrice;
}

int main(){
	const double spot = 100;
	const double strike = 100;
	const double r = 0.05;	
	const double vol = 0.3;
	const double yield = 0.02;
	const double expiry = 1;
	const double barrier = 120;
	const double b = 0.05;
	const int numStep = 252;
	const int numPath = 10000;
	const int numPilotSim = 1000;
	Option callUpAndOut = { spot, strike,r,expiry,vol,yield,barrier};
	//end input
	std::cout << "Closed-form barrier option price = ";
	double closedPrice = analytical(callUpAndOut);
	std::cout << "Monte Carlo barrier option price = ";
	double eulerMC = optionMC(callUpAndOut,numPath,numStep,(*euler));
	std::cout << "Antithetic variates barrier option price = ";
	double aMC = antitheticMC(callUpAndOut,numPath/2,numStep,(*euler));
	std::cout << "Control variates barrier option price = ";
	double cMC = controlMC(callUpAndOut, numPath-numPilotSim, numStep, numPilotSim, (*euler));




	system("pause");
	return 1;
}