// requires Random1.cpp

#include "Random1.h"
#include "Normals.h"
#include <iostream>
#include <cmath>
#include <string>

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

double analytical(double spot, double strike, double r, double expiry, double vol, double yield, double numPath, std::string type) {
	double var = vol*vol*expiry;
	double rootVar = sqrt(var);
	double itoCorrection = 0.5*var;

	double movedSpot = spot*exp(r*expiry - itoCorrection - yield);
	double thisSpot;
	double runningSum = 0;

	for (unsigned long i = 0; i < numPath; i++) {
		double thisGaussian = GetOneGaussianByBoxMuller();
		thisSpot = movedSpot*exp(rootVar*sqrt(expiry)*thisGaussian);
		double thisPayoff;
		if (type == "put") {
			thisPayoff = strike - thisSpot;
		}
		else if (type == "call") {
			thisPayoff = thisSpot - strike;
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

double euler(double spot, double strike, double r, double expiry, double vol, double yield, double numPath, int steps, std::string type) {
	double dt = 1 / static_cast<double>(steps);
	double runningSum = 0;
	for (unsigned long j = 0; j < numPath; j++) {
		double iterateSpot = spot;
		for (int i = 1; i <= steps; i++) {
			double dW = sqrt(dt) * GetOneGaussianByBoxMuller();
			double dS = iterateSpot*((r - yield)*dt + vol*dW);
			iterateSpot = iterateSpot + dS;
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

double eulerLog(double spot, double strike, double r, double expiry, double vol, double yield, double numPath, int steps, std::string type) {
	double dt = 1 / static_cast<double>(steps);
	double runningSum = 0;
	for (unsigned long j = 0; j < numPath; j++) {
		double iterateSpot = log(spot);
		for (int i = 1; i <= steps; i++) {
			double dW = sqrt(dt) * GetOneGaussianByBoxMuller();
			double dS = ((r - yield - pow(vol, 2) / 2)*dt + vol*dW);
			iterateSpot = iterateSpot + dS;
		}
		iterateSpot = exp(iterateSpot);
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

double milstein(double spot, double strike, double r, double expiry, double vol, double yield, double numPath, int steps, std::string type) {
	double dt = 1 / static_cast<double>(steps);
	double runningSum = 0;
	for (unsigned long j = 0; j < numPath; j++) {
		double iterateSpot = spot;
		double a = (r - yield)*iterateSpot;
		double b = vol*iterateSpot;
		for (int i = 1; i <= steps; i++) {
			double dW = sqrt(dt) * GetOneGaussianByBoxMuller();
			double dS = a*dt + b*dW + 0.5*b*vol*(pow(dW, 2) - dt);
			iterateSpot = iterateSpot + dS;
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

int main() {
	//begin input
	const double spot = 100;
	const double strike = 110;
	const double r = 0.05;
	const double vol = 0.3;
	const double yield = 0.02;
	const double expiry = 1;
	const int numSteps = 252;
	const int numPath = 10000;
	//end input
	int choice = 0;
	while (choice != 3) {
		int choice2 = 0;
		std::cout << "European Style Option" << std::endl;
		std::cout << "1. Call" << std::endl;
		std::cout << "2. Put" << std::endl;
		std::cout << "3. Exit" << std::endl;
		std::cin >> choice;
		switch (choice) {
		case 1:
			while (choice2 != 6) {
				std::cout << "European Call Option: " << std::endl;
				std::cout << "1. Closed-form formulae" << std::endl;
				std::cout << "2. Analytical Solution" << std::endl;
				std::cout << "3. Euler numerical solution" << std::endl;
				std::cout << "4. Euler numerical solution for log spot" << std::endl;
				std::cout << "5. Milstein numerical solution" << std::endl;
				std::cout << "6. Back" << std::endl;
				std::cin >> choice2;
				switch (choice2) {
				case 1:
					std::cout << "Option price using closed-form formula: " << closedPut(spot, strike, r, expiry, vol, yield, "call") << std::endl;
					break;
				case 2:
					std::cout << "Option price using single-step exact SDE solution: " << analytical(spot, strike, r, expiry, vol, yield, numPath, "call") << std::endl;
					break;
				case 3:
					std::cout << "Option price using Euler numerical solution of SDE for spot: " << euler(spot, strike, r, expiry, vol, yield, numPath, numSteps, "call") << std::endl;
					break;
				case 4:
					std::cout << "Option price using Euler numerical solution of SDE for log spot:5 " << eulerLog(spot, strike, r, expiry, vol, yield, numPath, numSteps, "call") << std::endl;
					break;
				case 5:
					std::cout << "Option price using Milstein numerical solution of SDE for spot: " << milstein(spot, strike, r, expiry, vol, yield, numPath, numSteps, "call") << std::endl;
					break;
				default:
					break;
				}
			}
			break;
		case 2:
			while (choice2 != 6) {
				std::cout << "European Put Option: " << std::endl;
				std::cout << "1. Closed-form formulae" << std::endl;
				std::cout << "2. Analytical Solution" << std::endl;
				std::cout << "3. Euler numerical solution" << std::endl;
				std::cout << "4. Euler numerical solution for log spot" << std::endl;
				std::cout << "5. Milstein numerical solution" << std::endl;
				std::cout << "6. Back" << std::endl;
				std::cin >> choice2;
				switch (choice2) {
				case 1:
					std::cout << "Option price using closed-form formula: " << closedPut(spot, strike, r, expiry, vol, yield, "put") << std::endl;
					break;
				case 2:
					std::cout << "Option price using single-step exact SDE solution: " << analytical(spot, strike, r, expiry, vol, yield, numPath, "put") << std::endl;
					break;
				case 3:
					std::cout << "Option price using Euler numerical solution of SDE for spot: " << euler(spot, strike, r, expiry, vol, yield, numPath, numSteps, "put") << std::endl;
					break;
				case 4:
					std::cout << "Option price using Euler numerical solution of SDE for log spot: " << eulerLog(spot, strike, r, expiry, vol, yield, numPath, numSteps, "put") << std::endl;
					break;
				case 5:
					std::cout << "Option price using Milstein numerical solution of SDE for spot: " << milstein(spot, strike, r, expiry, vol, yield, numPath, numSteps, "put") << std::endl;
					break;
				default:
					break;
				}
			}
			break;
		default:
			break;
		}
	}
	return 0;
}
