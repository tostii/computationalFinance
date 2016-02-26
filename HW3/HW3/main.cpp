#include <iostream>
#include "Option.h"
#include "LEcuyer.h"
#include "ParkMiller.h"
#include "Random2.h"
#include "Random1.h"
#include <fstream>
int main() {
	//defines parameters of european call option
	Option myOption;
	myOption.spot = 100;
	myOption.volatility = 0.3;
	myOption.interest = 0.05;
	myOption.yield = 0;
	myOption.strike = 110;
	myOption.expiry = 1;
	//define simulation parameters
	double numPaths = 10000;
	//generators
	RandomParkMiller generator(1);
	RandomLEcuyer LEgenerator(1, 1, 2, 3, 5, 7, 12);

	//(a) Uniform RNG 
	std::cout << "Closed-form option price = " << myOption.closedEuropeanCall() << std::endl;
	std::cout << "MC option price with Park-Miller uniform generator = " << myOption.monteCarloGenerator(numPaths, generator) << std::endl;
	std::cout << "MC option price with L'Ecuyer uniform generator = " << myOption.monteCarloGenerator(numPaths, LEgenerator) << std::endl;

	//(b) Fishman Acceptance-Rejection Algorithm
	std::cout << "MC option price with inverse distribution normal generator = " << myOption.monteCarloGaussian(numPaths, "inverse") << std::endl;
	std::cout << "MC option price with Box-Muller normal generator = " << myOption.monteCarloGaussian(numPaths, "park") << std::endl;
	std::cout << "MC option price with Fishman normal generator = " << myOption.monteCarloGaussian(numPaths, "fishman") << std::endl;
	//(b) 2.2 Generating Data
	std::ofstream myFile;
	myFile.open("Box-Muller.txt", std::ios::trunc);
	for (int i = 0; i < 10000; i++) {
		myFile << GetOneGaussianByBoxMullerPark() << std::endl;
	}
	myFile.close();
	myFile.open("Inverse-Transform.txt", std::ios::trunc);
	for (int i = 0; i < 10000; i++) {
		myFile << GetOneGaussianByInverseNormal() << std::endl;
	}
	myFile.close();
	myFile.open("Fishman.txt", std::ios::trunc);
	for (int i = 0; i < 10000; i++) {
		myFile << GetOneGaussianByFishman() << std::endl;
	}
	myFile.close();
	std::cin.get();
	return 1;
}