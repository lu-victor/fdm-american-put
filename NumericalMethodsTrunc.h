#ifndef NUMMETHODSTRUNC_H
#define NUMMETHODSTRUNC_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "DoubleSweep.hpp"
#include "option.h"
#include "CubicSpline.hpp"
#include <iomanip>

class NumericalMethodsTrunc {

public:
	double K;
	double S_0;
	double r;
	double q;
	double sigma;
	double T;
	//int M_prime;
	int M;
	int N;
	double S_max;
	double deltaS;
	double deltaT;

	double BCL;
	double BCR;

	std::vector<double> space_mesh;

	// constructor
	NumericalMethodsTrunc(double strike, double curr_price, double rate, double dividends, double vol, double mat, int space_div, int time_div) {
		K = strike;
		S_0 = curr_price;
		r = rate;
		q = dividends;
		sigma = vol;
		T = mat;
		//M_prime = space_div;
		N = time_div;
		M = space_div;

		S_max = 3 * strike;
		deltaS = S_max / space_div;
		//deltaS = S_0 / static_cast<double>(M_prime);
		
		//M = static_cast<int> ((3 * K) / deltaS);
		//S_max = M * deltaS;

		deltaT = T / static_cast<double>(N);
		BCL = strike;
		BCR = 0;

		space_mesh.resize(space_div + 1);
		for (int i = 0; i <= space_div; i++) {
			space_mesh[i] = i * deltaS;
		}
	}
	
	double explicitEuler();
	std::vector<double> explicitEulerVec();
	double implicitEuler();
	std::vector<double> implicitEulerVec();
	double thetaMethod(double theta);
	std::vector<double> thetaMethodVec(double theta);
	double ADE();
	std::vector<double> ADEVec();

	void printValues();

};




#endif