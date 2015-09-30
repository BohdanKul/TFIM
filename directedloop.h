#ifndef dloop_H
#define dloop_H

#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <array>
#include <cmath>

#include "vertex.h"
#include "helper.h"   // flipBit() and sort_indexes()

using namespace std;


array<float, 8> HeatBathSolution(const array<float,8>& Ws, const int& enLeg);

array<float, 8> BounceSolution(const array<float,8>& Ws, const int& enLeg);

array<float, 8> noBounceSolutionB(const array<float,8>& Ws, const int& enLeg);

int getSwitchLegP(const int& enleg, int vtype, array<bool,4>& clamped, array<float, 16>& VWeights, array<float,8>& prob);



#endif
