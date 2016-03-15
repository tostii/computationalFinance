#pragma once
#ifndef RANDOM1_H
#define RANDOM1_H
#include "ParkMiller.h"
#include "Normals.h"
double GetOneGaussianBySummation();
double GetOneGaussianByBoxMullerPark();
double GetOneGaussianByInverseNormal();
double GetOneGaussianByFishman();



#endif