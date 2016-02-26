#ifndef LECUYER_H
#define LECUYER_H
#include "Random2.h"
#include <math.h>
class LEcuyer
{
public:
	//1 2 3 
	//5 7 11
	LEcuyer(long x10, long x11, long x12, long x20, long x21, long x22);
	long GetOneRandomInteger();
	void SetSeed(long x10, long x11, long x12, long x20, long x21, long x22);
	unsigned long Max();
	unsigned long Min();
private:
	long x10, x11, x12, x13, x20, x21, x22;
};
class RandomLEcuyer: public RandomBase
{
public:
	RandomLEcuyer(unsigned long Dimensionality, long x10_, long x11_, long x12_, long x20_, long x21_, long x22_);
	RandomBase * clone() const;
	virtual void GetUniforms(MJArray& variates);
	virtual void Skip(unsigned long numberOfPaths);
	virtual void Reset();
	virtual void ResetDimensionality(unsigned long NewDimensionality);

private:
	LEcuyer InnerGenerator;
	unsigned long x10, x11, x12, x20, x21, x22;
	double Reciprocal;
};
#endif