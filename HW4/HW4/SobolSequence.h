#ifndef SOBOLSEQUENCE_H
#define SOBOLSEQUENCE_H
#include "Random2.h"

class SobolSequence : public RandomBase
{
public:

	SobolSequence(unsigned long Dimensionality, unsigned long Seed = 1UL);

	virtual RandomBase* clone() const;
	virtual void GetUniforms(MJArray& variates);
	virtual void Skip(unsigned long numberOfPaths);
	virtual void SetSeed(unsigned long Seed);
	virtual void Reset();
	virtual void ResetDimensionality(unsigned long NewDimensionality);

private:

	unsigned long InitialSeed;
	long long int seed;
	double *pdata;

};

#endif /* defined(__CFHW4__SobolSequence__) */