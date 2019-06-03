// TestTridiagonalSolvers.cpp
//
// Testing matrix solvers (LU decomposition) and Godunov
// Double Sweep method
//
// 2006-1-10 DD Kick-off
// 2017-2-13 DD new code
//
// (C) Datasim Education BV 2003-2017

#include "LUSolver.hpp"
#include "DoubleSweep.hpp" // They are templates
#include "Utilities.hpp"

#include <complex>
#include <iostream>
#include <cassert>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include "StopWatch.cpp"

/* Solve BVP u'' = 1 in (0, 1) with u(0) = u(1) = 0
   Solution u(x) = x(1-x)
   Solve by FDM */

int main()
{
	using value_type = double;
	using Vector = std::vector<value_type, boost::pool_allocator<double> >;
	
	std::size_t J = 20;
	std::cout << "Number of subdivisions J ";
	std::cin >> J;
	
	double h = 1.0 / static_cast<double>(J);

	// Boundary conditions
	double BCL = 0.0;	// LHS
	double BCR = 0.0;	// RHS

	// Double Sweep
	Vector a(J+1,1.0);
	Vector b(J+1,-2.0);
	Vector c(J+1,1.0);
	Vector r(J+1, -2.0*h*h);	// Right-hand side
	/*
	// Thomas algorithm
	Vector a2(J-1, 1.0);
	Vector b2(J-1, -2.0);
	Vector c2(J-1, 1.0);
	Vector r2(J-1, -2.0*h*h);	// Right-hand side

	// Take the boundary conditions into consideration
	r2[0] -= BCL;
	r2[r2.size()-1] -=  BCR;

	LUTridiagonalSolver<double> mySolver2(a2, b2, c2, r2);
	std::cout << "Matrix has a solution? " << mySolver2.diagonallyDominant() << '\n';
	DoubleSweep<value_type> mySolver(a, b, c, r, BCL, BCR);

	StopWatch<> sw;
	sw.Start();
	Vector result = std::move(mySolver());
	Vector result2 = std::move(mySolver2());
	sw.Stop();
	std::cout << "Elapsed time: " << sw.GetTime() << '\n';

	auto exact = [](double x) { return x*(1.0 - x); };
	double val = h;		// Double Sweep

	// Comapare output from Double Sweep and Thomas
	for (std::size_t j = 1; j < result.size()-1; ++j)
	{ // The values shoule be zero

		std::cout << j << ", " << result[j]-result2[j-1] << ", " << exact(val) << '\n';
		val += h;
	}
	*/
	{	// Boundary conditions
		double BCL = 0.0;	// LHS
		double BCR = 0.0;	// RHS

							// Double Sweep
		Vector a(J + 1, 1.0);
		Vector b(J + 1, -2.0);
		Vector c(J + 1, 1.0);
		Vector r(J + 1, -2.0*h*h);	// Right-hand side

		// Double Sweep with a memory allocator
		/*#include <boost/pool/pool.hpp>
		#include <boost/pool/pool_alloc.hpp>*/
		DoubleSweep<double, std::vector, boost::pool_allocator<double>> mySolver3(a, b, c, r, BCL, BCR);
		Vector result = mySolver3();
		auto exact = [](double x) { return x*(1.0 - x); };
		double val = 0;		// Double Sweep

		// Comapare output from Double Sweep and Thomas
		for (std::size_t j = 0; j < result.size(); ++j)
		{ // The values shoule be zero

			std::cout << j << ", " << result[j] << ", " << exact(val) << '\n';
			val += h;
		}
	}

	return 0;
}
