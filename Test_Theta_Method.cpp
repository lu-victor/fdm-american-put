#include "DoubleSweep.hpp" // They are templates

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <random>

#include <fstream>

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
//				double (*f) (douible x, double t)
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

#include <cmath>

constexpr double TOL = 1.0e-9; // accuracy desired
constexpr double pi = 3.141592;


double ExactSolution(double x, double t)
{
	using std::sin; using std::exp;
	// Compute the number of terms
	int N = static_cast<int>(std::sqrt(1.0 / TOL));
	double result = 0.0;
	for (int n = 1; n <= N; ++n)
	{
		result += (sin(n*pi / 2)*sin(n*pi*x)*exp(-n * n*pi*pi*t)) / (n*n);
	}

	return 8.0*result / (pi*pi);
}

std::vector<double> ExactSolutionArray(const std::vector<double>& xarr, double t)
{
	std::vector<double> result(xarr.size());
	for (std::size_t j = 0; j < result.size(); ++j)
	{
		result[j] = ExactSolution(xarr[j], t);
	}

	return result;
}
double InitialCondition(double x)
{
	if (x >= 0.0 && x <= 0.5)
	{
		return 2.0 * x;
	}

	return 2.0 * (1.0 - x);
}


int main()
{

	// BTCS scheme for the heat equation

	long J = 20;

	std::cout << "Number of space subdivisions: ";
	std::cin >> J;

	std::cout << "Number of time subdivisions: ";
	long N = 100; std::cin >> N;

	std::cout << "Theta method:  1) Implicit Euler, 2) Crank Nicolson: ";
	long choice = 1; std::cin >> choice;
	double theta = 1.0; if (2 == choice) theta = 0.5;

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

	double lambda = k / (h*h);


	// Dirichlet boundary conditions
	std::vector<double> a(J + 1, -lambda * theta);
	std::vector<double> b(J + 1, (1.0 + 2.0*lambda*theta));
	std::vector<double> c(J + 1, -lambda * theta);
	std::vector<double> r(J + 1, 0);	// Right-hand side NOT CONSTANT ANYMORE
								// Take the boundary conditions into consideration

	// Create mesh in space
	std::vector<double> xarr(J + 1);	xarr[0] = A;
	for (std::size_t j = 1; j < xarr.size(); ++j)
	{
		xarr[j] = xarr[j - 1] + h;
	}


	// Initial condition
	std::vector<double> vecOld(xarr.size());	// At time n
	std::vector<double> vecNew(xarr.size());	// At time n+1

	for (std::size_t j = 0; j < vecOld.size(); j++)
	{
		vecOld[j] = InitialCondition(xarr[j]);
	}

	// We start at 1st time point
	double current = k;


	while (current <= T)
	{
		//	std::cout << current << ",";
			// Update at new time level n+1

			// Compute inhomogeneous term
		for (std::size_t j = 1; j < r.size() - 1; ++j)
		{
			r[j] = (lambda*(1.0 - theta)*vecOld[j + 1]) + (1.0 - (2.0*lambda*(1.0 - theta)))*vecOld[j]
				+ (lambda*(1.0 - theta)*vecOld[j - 1]);
		}


		DoubleSweep<double> mySolver(a, b, c, r, BCL, BCR);
		
		vecNew = mySolver.solve();
		vecOld = vecNew;

		current += k;
	}


	// Examine the error between exact and approximate solutions
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

	outfile.open("solutions_CN.csv");
	outfile << "x,exact,approx_solution" << std::endl;

	for (int i = 0; i < xarr.size(); i++) {
		outfile << xarr[i] << "," << exact[i] << "," << vecNew[i] << std::endl;
	}

	outfile.close();

}