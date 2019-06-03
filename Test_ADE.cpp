

// #include "Utilities.hpp"
// !!!!       #include "ExcelDriverLite.hpp"
// #include "StopWatch.cpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <random>
#include <fstream>

template <typename T>
using Vector = std::vector<T>;

//
//	See chapter 4 of FDM book for Fourier series solution
// 
// Use special code names for specific exact solutions
//
// We have three functions for presentation
//


//	1. Get the solution at one point (x,t)
//									: output double
//	2. Get the solution at a space array and one time (xarr, t)
//									: output Vector
//	3. Get solution for a space array for a number of time levels
//									: output list<Vector>
//
// This idea can be applied to other problems. For example, we could
// use function pointers (even better is std::function<>)
//	
//				double (*f) (double x, double t)
//
// Here are the steps for the IBV problem
//
//		U_t = U_xx for A < x < B, t > 0
//		
//		BC U(A,t) = U(B,t) = 0, t > 0
//		IC
//			U(x,0) = 2x for A <= x <= 0.5
//			U(x,0) = 2(1-x) for 0.5 <= x <= B
//
// In code, take A = 0, B = 1.
//
// Exact solution
//

constexpr double pi = 3.14159265358979;
constexpr int NTrunc = 50;

// First test set

double ExactSolution(double x, double t)
{
	using std::sin; using std::exp;
	double result = 0.0;
	for (int n = 1; n <= NTrunc; ++n)
	{
		result += (sin(n*pi/2)*sin(n*pi*x)*exp(-n*n*pi*pi*t))/(n*n);
	}

	return 8.0*result / (pi*pi);
}

double InitialCondition(double x)
{
	if (x >= 0.0 && x <= 0.5)
	{
		return 2.0 * x;
	}
	else
	{
		return 2.0 * (1.0 - x);
	}

}

double bcl(double t)
{
	return 0.0;
}

double bcr(double t)
{
	return 0.0;
}


/* Second test set
double ExactSolution(double x, double t)
{
	return std::exp(x + t);
}

double InitialCondition(double x)
{
	return std::exp(x);
}

double bcl(double t)
{
	return std::exp(t);
}

double bcr(double t)
{
	return std::exp(1.0 + t);
}
*/

Vector<double> ExactSolutionArray(const Vector<double>& xarr, double t)
{
	Vector<double> result(xarr.size());
	for (std::size_t j = 0; j < result.size(); ++j)
	{
		result[j] = ExactSolution(xarr[j], t);
	}

	return result;
}


int main()
{
	/*std::cout << "*** TIPS ***\n";
	std::cout << "For testing, choose NT = 5 NX\n";
	std::cout << "For not very smooth initial conditon and/or T 'big' choose NT = 1000\n";
	std::cout << "*** END TIPS ***\n\n";
	*/

	long J = 20;

	std::cout << "Number of space subdivisions: ";
	std::cin >> J;

	std::cout << "Number of time subdivisions: ";
	long N = 100; std::cin >> N;

	std::cout << "Expiration: ";
	double T = 1.0; std::cin >> T;

	// Space interval [A,B]
	double A = 0.0; // LHS
	double B = 1.0;	// RHS
	double h = (B - A) / static_cast<double>(J);
	double k = T / static_cast<double>(N);


	// Boundary conditions
	double BCL = 0.0;
	double BCR = 0.0;


	std::vector<double> xarr(J + 1);	xarr[0] = A;
	for (std::size_t j = 1; j < xarr.size(); ++j)
	{
		xarr[j] = xarr[j - 1] + h;
	}

	// Initial condition
	Vector<double> U(xarr.size());	// L-R sweep
	Vector<double> V(xarr.size());	// R-L sweep

	Vector<double> vecNew(xarr.size());	// At time n+1

	// we know the form of the solutions at t=0
	for (std::size_t j = 0; j < U.size(); j++)
	{
		U[j] = V[j] = InitialCondition(xarr[j]);
	}

	// We start at 1st time point
	double current = k;

	// Parameters for FD scheme
	double lambda = k / (h*h);
	double OnePlusLambda = 1.0 / (1.0 + lambda);
	double OneMinusLambda = 1.0 - lambda;

	while (current <= T)
	{
		//	std::cout << current << ",";
		// Update at new time level n+1

		// (Dirichlet) boundary conditions
		vecNew[0] = U[0] = V[0] = bcl(current);
		vecNew[vecNew.size() - 1] = U[vecNew.size() - 1] = V[vecNew.size() - 1] = bcr(current);

		// Up Sweep
		for (std::size_t j = 1; j < U.size() - 1; ++j) // 1..J-1
		{
			U[j] = ((U[j] * OneMinusLambda) + lambda * (U[j + 1] + U[j - 1])) * OnePlusLambda;
		}

		// Down Sweep
		for (std::size_t j = V.size() - 2; j >= 1; --j)  // J-1..1
		{
			V[j] = ((V[j] * OneMinusLambda) + lambda * (V[j + 1] + V[j - 1])) * OnePlusLambda;
		}

		// Barakat and Clark update
		for (std::size_t j = 0; j < vecNew.size(); ++j)
		{
			vecNew[j] = 0.5*(U[j] + V[j]);
		}

		current += k;
	}

	// Examine the error between exact and approximate solutions
	std::cout << "Current time level: " << current-k << '\n';
	std::cout << "Post processing, could take some time\n";

	auto exact = ExactSolutionArray(xarr, T);

	std::cout << "Our solution at T: " << std::endl;
	for (auto elem : vecNew)
	{
		std::cout << elem << ",";
	}
	std::cout << '\n';
	std::cout << "The exact solution at T: " << std::endl;
	for (auto elem : exact)
	{
		std::cout << elem << ",";
	}

	std::ofstream outfile;

	outfile.open("solutions_ADE.csv");
	outfile << "x,exact,approx_solution" << std::endl;

	for (int i = 0; i < xarr.size(); i++) {
		outfile << xarr[i] << "," << exact[i] << "," << vecNew[i] << std::endl;
	}

	outfile.close();

	// Excel display
/*	ExcelDriver xl; xl.MakeVisible(true);
	xl.CreateChart(xarr, exact, "Exact solution");
	xl.CreateChart(xarr, vecNew, "ADE solution");*/
}