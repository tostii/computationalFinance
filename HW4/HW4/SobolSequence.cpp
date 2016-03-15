#include "SobolSequence.h"
#include "sobol.h"

const long a = 16807;
const long m = 2147483647;
const long q = 127773;
const long r = 2836;

SobolSequence::SobolSequence(unsigned long Dimensionality, unsigned long Seed)
	: RandomBase(Dimensionality),
	InitialSeed(Seed)
{
	seed = Seed;
	pdata = new double[Dimensionality];
}

RandomBase* SobolSequence::clone() const
{
	return new SobolSequence(*this);
}

void SobolSequence::GetUniforms(MJArray& variates)
{
	int d = GetDimensionality();
	i8_sobol(d, &seed, pdata);//I use i8-Sobol
	for (unsigned long j = 0; j < GetDimensionality(); j++)
		variates[j] = pdata[j];
}

void SobolSequence::Skip(unsigned long numberOfPaths)
{
	MJArray tmp(GetDimensionality());
	for (unsigned long j = 0; j < numberOfPaths; j++)
		GetUniforms(tmp);
}

void SobolSequence::SetSeed(unsigned long Seed)
{
	InitialSeed = Seed;
}

void SobolSequence::Reset()
{
	seed = InitialSeed;
}


void SobolSequence::ResetDimensionality(unsigned long NewDimensionality)
{
	RandomBase::ResetDimensionality(NewDimensionality);
	delete[] pdata;
	pdata = new double[NewDimensionality];
	seed = InitialSeed;
}