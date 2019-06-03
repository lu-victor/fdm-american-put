// Project_FDM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "NumericalMethodsTrunc.h"
#include "NumericalMethods.h"
#include <fstream>

int main() {
	double K = 50;
	double S_0 = 50;
	double r = 0.10;
	double q = 0;
	double sigma = 0.40;
	double T = 1;

	// space subdivisions
	//int M_prime = 10;
	


	//NumericalMethodsTrunc testNumTrunc(K, S_0, r, q, sigma, T, M, N);
	//testNumTrunc.printValues();

	//int stepsPut = 1000;
	
	
	//std::cout << "Binomial Tree: " << std::endl;
	//std::cout << "Time subdivisions: " << stepsPut << std::endl;
	//std::cout << "Put value: " << my_Put.getBinomialTreeValue(S_0, stepsPut) << std::endl;
	//std::cout << "----------------------------------------------------------" << std::endl;
	

	//NumericalMethods testNum(K, S_0, r, q, sigma, T, M, N);
	/*std::cout << testNum.explicitEuler() << std::endl;
	std::cout << testNum.implicitEuler() << std::endl;
	std::cout << testNum.thetaMethod(0.5) << std::endl;*/
	//std::cout << testNum.ADE() << std::endl;
	//testNum.printValues();

	// create text files with values obtained

	std::ofstream outfile;

	/*outfile.open("solutions_binomial.csv");
	outfile << "steps,binom_value" << std::endl;

	std::vector<double> time_steps(12);
	time_steps[0] = 2;
	for (int i = 1; i < time_steps.size(); i++) {
		time_steps[i] = time_steps[i-1] * 2;
	}

	AmericanPut my_Put(K, T, sigma, r);
	for (int i = 0; i < time_steps.size(); i++) {
		outfile << time_steps[i] << "," << my_Put.getBinomialTreeValue(S_0, time_steps[i]) << std::endl;
	}

	outfile.close();
	*/

	// Explicit Euler truncated
	
	// same number of steps
	/*
	std::vector<double> M(9);
	M[0] = 20;
	for (int i = 1; i < M.size(); i++) {
		M[i] = M[i - 1] * 2;
	}
	outfile.open("solutions_expEuler_trunc.csv");
	outfile << "M,N,expEuler_trunc" << std::endl;
	for (int i = 0; i < M.size(); i++) {
		NumericalMethodsTrunc testNumTrunc(K, S_0, r, q, sigma, T, M[i], M[i]);
		outfile << M[i] << "," << M[i] << "," << testNumTrunc.explicitEuler() << std::endl;
	}
	outfile.close();
	*/

	/* Explicit Euler, constant space step, decrease time step
	std::vector<double> N(9);
	N[0] = 20;
	for (int i = 1; i < N.size(); i++) {
		N[i] = N[i - 1] * 2;
	}
	outfile.open("expEuler_trunc_m_n_diff.csv");
	outfile << "M,N,expEuler_trunc" << std::endl;
	for (int i = 0; i < N.size(); i++) {
		NumericalMethodsTrunc testNumTrunc(K, S_0, r, q, sigma, T, 500, N[i]);
		outfile << "500" << "," << N[i] << "," << testNumTrunc.explicitEuler() << std::endl;
	}
	outfile.close();
	*/
	/*
	std::vector<double> M(9);
	M[0] = 20;
	for (int i = 1; i < M.size(); i++) {
		M[i] = M[i - 1] * 2;
	}
	outfile.open("expEuler_m=n.csv");
	outfile << "M,N,expEuler" << std::endl;
	for (int i = 0; i < M.size(); i++) {
		NumericalMethods testNum(K, S_0, r, q, sigma, T, M[i], M[i]);
		outfile << M[i] << "," << M[i] << "," << testNum.explicitEuler() << std::endl;
	}
	outfile.close();
	*/

	/*
	std::vector<double> M(9);
	M[0] = 20;
	for (int i = 1; i < M.size(); i++) {
		M[i] = M[i - 1] * 2;
	}
	outfile.open("data.csv");
	outfile << "M,N,expEuler_trunc,impEuler_trunc,CN_trunc,ADE_trunc,expEuler,impEuler,CN,ADE" << std::endl;
	for (int i = 0; i < M.size(); i++) {
		NumericalMethodsTrunc testNumTrunc(K, S_0, r, q, sigma, T, M[i], M[i]);
		NumericalMethods testNum(K, S_0, r, q, sigma, T, M[i], M[i]);
		outfile << M[i] << "," << M[i] << "," << testNumTrunc.explicitEuler() << "," << testNumTrunc.implicitEuler() 
			<< "," << testNumTrunc.thetaMethod(0.5) << "," << testNumTrunc.ADE() << "," <<
			testNum.explicitEuler() << "," << testNum.implicitEuler()
			<< "," << testNum.thetaMethod(0.5) << "," << testNum.ADE() << std::endl;
	}
	outfile.close();

	*/

	/*
	// returns values of put in grid
	
	// constant steps, we vary the initial stock price
	double M = 100;
	double N = 100;
	AmericanPut my_Put(K, T, sigma, r);

	outfile.open("data_truncM100N100.csv");
	outfile << "S_0,expEuler_trunc,impEuler_trunc,CN_trunc,ADE_trunc,BT" << std::endl;
	NumericalMethodsTrunc testNumTrunc(K, S_0, r, q, sigma, T, M, N);
	std::vector<double> initial_prices = testNumTrunc.space_mesh;
	std::vector<double> expEuler_trunc = testNumTrunc.explicitEulerVec();
	std::vector<double> impEuler_trunc = testNumTrunc.implicitEulerVec();
	std::vector<double> CN_trunc = testNumTrunc.thetaMethodVec(0.5);
	std::vector<double> ADE_trunc = testNumTrunc.ADEVec();
	for (int i = 0; i < initial_prices.size(); i++) {
		outfile << initial_prices[i] << "," << expEuler_trunc[i] << "," << impEuler_trunc[i] << "," <<
			CN_trunc[i] << "," << ADE_trunc[i] << "," << my_Put.getBinomialTreeValue(initial_prices[i],1000) << std::endl;
	}

	outfile.close();

	*/

	// constant steps, we vary the initial stock price
	double M = 500;
	double N = 100;
	AmericanPut my_Put(K, T, sigma, r);

	outfile.open("data_M500N100.csv");
	outfile << "S_0,expEuler,impEuler,CN,ADE,BT" << std::endl;
	NumericalMethods testNum(K, S_0, r, q, sigma, T, M, N);
	std::vector<double> initial_prices = testNum.S;
	std::vector<double> expEuler = testNum.explicitEulerVec();
	std::vector<double> impEuler = testNum.implicitEulerVec();
	std::vector<double> CN = testNum.thetaMethodVec(0.5);
	std::vector<double> ADE = testNum.ADEVec();
	for (int i = 0; i < initial_prices.size(); i++) {
		outfile << initial_prices[i] << "," << expEuler[i] << "," << impEuler[i] << "," <<
			CN[i] << "," << ADE[i] << "," << my_Put.getBinomialTreeValue(initial_prices[i], 250) << std::endl;
	}

	outfile.close();

	return 0;
}