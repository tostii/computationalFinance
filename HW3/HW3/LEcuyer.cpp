#include "LEcuyer.h"

const long m1 = 2147483647;
const long m2 = 2145483479;
const long a11 = 0;
const long a12 = 63308;
const long a13 = -183326;
const long a21 = 86098;
const long a22 = 0;
const long a23 = -539608;

long q12 = floor(m1 / a12);
long q13 = floor(m1 / a13);
long q21 = floor(m2 / a21);
long q23 = floor(m2 / a23);
long r12 = m1 - a12 * q12;
long r13 = m1 - a13 * q13;
long r21 = m2 - a21 * q21;
long r23 = m2 - a23 * q23;
float k = 1 / (m1 + 1);



LEcuyer::LEcuyer(long x10_, long x11_, long x12_, long x20_, long x21_, long x22_){
	x10 = x10_;
	x11 = x11_;
	x12 = x12_;
	x20 = x20_;
	x21 = x21_;
	x22 = x22_;
}

void LEcuyer::SetSeed(long x10_, long x11_, long x12_, long x20_, long x21_, long x22_) {
	x10 = x10_;
	x11 = x11_;
	x12 = x12_;
	x20 = x20_;
	x21 = x21_;
	x22 = x22_;
}

unsigned long LEcuyer::Max()
{
	return m1 - 1;
}
unsigned long LEcuyer::Min()
{
	return 1;
}

long LEcuyer::GetOneRandomInteger()
{
	//computer a12*x11 mod m1
	long h = floor(x11 / q12);
	long p12 = a12 * (x11 - h * q12) - h * r12;
	if (p12 < 0) {
		p12 = p12 + m1;
	}
	//compute a13 * x10 mod m1
	h = floor(x10 / q13);
	long p13 = a13 * (x10 - h * q13) - h * r13;
	if (p13 < 0) {
		p13 = p13 + m1;
	}
	//compute a21 * x22 mod m2
	h = floor(x22 / q21);
	long p21 = a21 * (x22 - h*q21) - h*r21;
	if (p21 < 0) {
		p21 = p21 + m2;
	}
	//compute a23 * x20 mod m2
	h = floor(x20 / q23);
	long p23 = a23 * (x20 - h * q23) - h * r23;
	if (p23 < 0) {
		p23 = p23 + m2;
	}
	//update x11 and x10
	x10 = x11;
	x11 = x12;
	//compute (p12 + p13) mod m1 and update x12
	p12 = p12 - m1;
	x12 = p12 + p13;
	if (x12 < 0) {
		x12 = x12 + m1;
	}
	//update x21 and x20
	x20 = x21;
	x21 = x22;
	//compute (p21 + p23) mod m2 and update x22
	p21 = p21 - m2;
	x22 = p21 + p23;
	if (x22 < 0) {
		x22 = x22 + m2;
	}
	//compute x1 mod m1
	long x1;
	if (x12 < x22) {
		x1 = x12 - x22 + m1;
	}
	else {
		x1 = x12 - x22;
	}
	//compute u
	if (x1 == 0) {
		x1 = m1;
	}
	return x1;
}

RandomLEcuyer::RandomLEcuyer(unsigned long Dimensionality, long x10_, long x11_, long x12_, long x20_, long x21_, long x22_) : RandomBase(Dimensionality), InnerGenerator(x10_, x11_, x12_, x20_, x21_, x22_), x10(x10_), x11(x11_), x12(x12_), x20(x20_), x21(x21_), x22(x22_){
	Reciprocal = 1 / (1.0 + InnerGenerator.Max());
}

void RandomLEcuyer::GetUniforms(MJArray& variates) {
	for (unsigned long j = 0; j < GetDimensionality(); j++) {
		variates[j] = InnerGenerator.GetOneRandomInteger()*Reciprocal;
	}
}

void RandomLEcuyer::Skip(unsigned long numberOfPaths) {
	MJArray tmp(GetDimensionality());
	for (unsigned long j = 0; j < numberOfPaths; j++) {
		GetUniforms(tmp);
	}
}

RandomBase* RandomLEcuyer::clone() const {
	return new RandomLEcuyer(*this);
}

void RandomLEcuyer::Reset() {
	InnerGenerator.SetSeed(x10, x11, x12, x20, x21, x22);
}

void RandomLEcuyer::ResetDimensionality(unsigned long NewDimensionality) {
	RandomBase::ResetDimensionality(NewDimensionality);
	InnerGenerator.SetSeed(x10, x11, x12, x20, x21, x22);
}
