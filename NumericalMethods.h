#ifndef NUMMETHODS_H
#define NUMMETHODS_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "DoubleSweep.hpp"
#include "option.h"
#include "CubicSpline.hpp"
#include <iomanip>

class NumericalMethods {

public:
	double K;
	double S_0;
	double y_0;
	double r;
	double q;
	double sigma;
	double T;
	
	int M;
	int N;

	double deltaY;
	double deltaT;

	double BCL;
	double BCR;

	std::vector<double> space_mesh_Y;
	
	std::vector<double> S;

	// constructor
	NumericalMethods(double strike, double curr_price, double rate, double dividends, double vol, double mat, int space_div, int time_div) {
		K = strike;
		S_0 = curr_price;
		r = rate;
		q = dividends;
		sigma = vol;
		T = mat;
		//M_prime = space_div;
		N = time_div;
		M = space_div;

		deltaY = 1.0 / static_cast<double>(space_div);

		deltaT = T / static_cast<double>(N);
		BCL = strike;
		BCR = 0;

		space_mesh_Y.resize(space_div + 1);
		for (int i = 0; i <= space_div; i++) {
			space_mesh_Y[i] = i * deltaY;
		}

		// vector of S associated with Y
		S.resize(space_div + 1);
		S[0] = 0;
		for (int i = 1; i <= M - 1; i++) {
			S[i] = (i * deltaY) / (1 - i * deltaY);
		}

		y_0 = S_0 / (S_0 + 1.0);
	}

	double explicitEuler();
	double implicitEuler();
	double thetaMethod(double theta);
	double ADE();

	std::vector<double> explicitEulerVec();
	std::vector<double> implicitEulerVec();
	std::vector<double> thetaMethodVec(double theta);
	std::vector<double> ADEVec();

	void printValues();

};




#endif