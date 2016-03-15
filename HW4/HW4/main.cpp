#include "Option.h"
#include "ParkMiller.h"
#include "AntiThetic.h"
#include "SobolSequence.h"
#include <iostream>
#include <fstream>
int main() {
	//european option
	Option vanillaEuropean;
	vanillaEuropean.volatility = 0.3;
	vanillaEuropean.spot = 100;
	vanillaEuropean.interest = 0.05;
	vanillaEuropean.yield = 0;
	vanillaEuropean.strike = 110;
	vanillaEuropean.expiry = 1;
	//paths
	double numPaths = 10000;
	//generators
	RandomParkMiller generatorPM(1);
	AntiThetic antitheticPM(generatorPM);
	SobolSequence generatorSOBOL(2);
	//get data
	/** Data Plotting Generating Code
	std::ofstream myX;
	std::ofstream myY;
	MJArray VariateArray(2);
	myX.open("X.txt", std::ios::trunc);
	myY.open("Y.txt", std::ios::trunc);
	for (int i = 0; i < 1024; i++) {
		generatorSOBOL.GetUniforms(VariateArray);
		myX << VariateArray[0] << std::endl;
		myY << VariateArray[1] << std::endl;

	}
	myX.close();
	myY.close();
	**/
	//print
	std::cout << "Closed-form vanilla call option price = " << vanillaEuropean.closedEuropeanCall() << std::endl;
	std::cout << "MC vanilla call price with Park-Miller uniforms = " << vanillaEuropean.monteCarloGenerator(numPaths, generatorPM) << std::endl;
	std::cout << "MC vanilla call price with Park-Miller uniforms and antithetics = " << vanillaEuropean.monteCarloGenerator(numPaths,antitheticPM) << std::endl;
	std::cout << "QMC vanilla call price with Sobol sequence = " << vanillaEuropean.monteCarloGenerator(numPaths, generatorSOBOL) << std::endl;
	std::cin.get();
	return 1;
}