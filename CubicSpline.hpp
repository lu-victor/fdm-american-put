// CubicSpline.hpp
//
// Class to represent cubic spline interpolation
//
// Code is default inline and we include some C utility functions
//
// Last Modification Dates:
//
// DD 2006-7-31 Kick-off code
// DD 2006-8-1 Tested: 1) BC 2nd order terms
// DD 2017-2-17 C++11 update + core functionality
// DD 2017-7-1 cubic spline as a function object
//
// (C) Datasim Education BV 2006-2018
//

// All equations based on my C++ book 1st edition, especially chapter 18

// Description: Given a set of mesh points Xarr (0 to N), a set of function values Yarr
// (0 to N) at these mesh points find the cubic spline function
// that agrees with Yarr at the mesh points and having either of the Boundary conditions.

#ifndef CUBICSPLINE_HPP
#define CUBICSPLINE_HPP

#include "LUSolver.hpp" // Process involves solution of a tridiagonal matrix
#include "Utilities.hpp"
#include <vector> 
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

/* STEPS

	1. Give Xarr and Yarr (Same dimensions!)
	2. Give B.C. 
	3. Calculate the sub, main and superdiagonal arrays of LU matrix
	4. Solve AM = b 
	5. Choose the abscissa value 'x' where you want the interpolant.
	6. Use this value in the formula for the cubic spline S
*/

// 2018-3-5 DD method for 2nd derivative

enum CubicSplineBC {SecondDeriv, FirstDeriv};

class CubicSplineInterpolator
{
private:
	std::vector<double> x;		// Abscissa x-values x_0,...,x_N
	std::vector<double> y;		// Function values
	std::vector<double> h;		// Array of mesh sizes in x direction

	CubicSplineBC type;			// Type of BC

	std::vector<double> M;		// Special calculated coefficients of spline
	std::vector<double> A, B, C, r;	// Input arrays for LU decomposition

	// For first order derivatives
	double a, b;

private:
		// Private member functions
		void CalculateVectors()
		{ // A, B, C and r

			std::size_t N = x.size()-1;

			if (type == SecondDeriv)
			{ 
				C[0] = 0.0;
				r[0] = 0.0;

				A[A.size()-1] = 0.0;
				r[r.size()-1] = 0.0;
			}
			else
			{
				C[0] = 1.0;
				r[0] = 6.0 * ((y[1] - y[0])/h[1] - a)/h[1];

				A[A.size()-1] = 1.0;
				r[r.size()-1] = 6.0 * (b - ((y[N] - y[N-1])/h[N]))/h[N];
			}
		
			double tmp;

			for (std::size_t j = 1; j < x.size()-1; ++j)
			{ // Optimise later

				double fac = 1.0 / (h[j] + h[j + 1]);
				C[j] = h[j+1] * fac;
				A[j] = h[j] * fac;

				tmp = ((y[j+1] - y[j])/h[j+1]) - ((y[j] - y[j-1])/h[j]);

				r[j] = (6.0 * tmp) * fac;
			}

		}
	
	std::size_t findSAbscissa(double xvar) const
	{ // Will give index of LHS value <= x. 
		if (xvar < x[0] || xvar >  x[x.size() - 1])
		{
			std::string s = "\nValue " + boost::lexical_cast<std::string>(xvar) + " not in range "
				+ "(" + boost::lexical_cast<std::string>(x[0]) + ","
				+ boost::lexical_cast<std::string>(x[x.size() - 1]) + ")";
			throw std::out_of_range(s);
		}

		auto posA = std::lower_bound(std::begin(x), std::end(x), xvar); // Log complexity
			
		std::size_t index = std::distance(std::begin(x), posA);

		return index;
	}

public:
	CubicSplineInterpolator(const std::vector<double>& xarr,
							const std::vector<double>& yarr,
							CubicSplineBC BCType, 
							double alpha = 0.0,	double beta = 0.0)
	{ // Discrete function value case

		// Arrays must have the same size
		x = xarr;
		y = yarr;
		type = BCType;

	
		a = alpha;	// LHS
		b = beta;	// RHS
		std::size_t N = xarr.size();
	
		// Calculate array of offset
		h = std::vector<double>(N, 0.0);
		for (std::size_t j = 1; j < h.size(); ++j)
		{
			h[j] = x[j] - x[j-1];
		}

		// All arrays have start index 1
		// Compared to the equations in the book, M(j) --> M(j+1)
		M = std::vector<double>(N, 0.0); // Solution

		// LU Coefficients
		A = std::vector<double>(N, 0.0);
		B = std::vector<double>(N, 2.0);	// Diagonal vector, constant == 2
		C = std::vector<double>(N, 0.0);
		r = std::vector<double>(N, 0.0);
		
		// Calculate the elements 
		CalculateVectors();

		LUTridiagonalSolver<double> mySolver(A, B, C, r);
		
		M = mySolver.solve();
	}

	CubicSplineInterpolator(const std::vector<double>& xarr,
							const std::function<double (double)>& fun,
							CubicSplineBC BCType,
							double alpha = 0.0, double beta = 0.0)
	{ // Continuous function value case

		std::size_t N = xarr.size();

		// Arrays must have the same size
		x = xarr;
//		print(x);
		y = CreateDiscreteFunction(xarr, fun);
	//	print(y);
		type = BCType;


		a = alpha;	// LHS
		b = beta;	// RHS
		
		// Calculate array of offset
		h = std::vector<double>(N, 0.0);
		for (std::size_t j = 1; j < h.size(); ++j)
		{
			h[j] = x[j] - x[j - 1];
		}

		// All arrays have start index 1
		// Compared to the equations in the book, M(j) --> M(j+1)
		M = std::vector<double>(N, 0.0); // Solution

		// LU Coefficients
		A = std::vector<double>(N, 0.0);
		B = std::vector<double>(N, 2.0);	// Diagonal vector, constant == 2
		C = std::vector<double>(N, 0.0);
		r = std::vector<double>(N, 0.0);

		// Calculate the elements 
		CalculateVectors();

		LUTridiagonalSolver<double> mySolver(A, B, C, r);

		M = mySolver.solve();
	}

	double Solve(double xvar) const
	{ // Find the interpolated valued at a value x)

		std::size_t j = findAbscissa(x, xvar);	// will give index of LHS value <= x
		// Now use the formula
		double tmp = xvar - x[j];
		double tmpA = x[j+1] - xvar;

		double tmp3 = tmp * tmp * tmp;
		double tmp4 = tmpA * tmpA * tmpA;

		double A = (y[j+1] - y[j])/h[j+1] - (h[j+1] * (M[j+1] - M[j]))/6.0;
		double B = y[j] - (M[j] * h[j+1] * h[j+1])/6.0; 

		double result = (M[j] * tmp4)/(6.0 * h[j+1])
							+ (M[j+1] * tmp3)/(6.0 * h[j+1])
								+ (A * tmp)
									+ B;
		return result;

	}
	
	double operator () (double xvar) const
	{
		return Solve(xvar);
	}
	
	std::vector<double> Curve(const std::vector<double>& xarr) const
	{ // Create the interpolated curve

//		print(xarr);
		std::vector<double> result(xarr.size());

		for (std::size_t j = 0; j < xarr.size(); ++j)
		{

			result[j] = Solve(xarr[j]);
		}

	
		return result;
	}

	std::vector<double> Curve() const
	{ // Create the interpolated curve, MEMBER DATA AS ABSCISSAE

		return Curve(x);
	}
	
	double Integral() const
	{ // ANW page 45

		double r1 = 0.0; double r2 = 0.0;

		for (std::size_t j = 1; j < x.size(); ++j)
		{
			double hj = h[j];
			r1 += (y[j - 1] + y[j]) * hj;
			r2 -= (M[j - 1] + M[j]) * hj * hj * hj;
		}

		r1 *= 0.5;
		r2 /= 24.0;

		return r1 + r2;
	}

	std::tuple<double, double, double> ExtendedSolve(double xvar)
	{ // Solve for S and derivatives S', S"

		auto j = findAbscissa(x, xvar);

		double Mj = M[j]; double Mjp1 = M[j + 1];
		
		double hjp1 = 1.0 / h[j + 1];
		double xj = x[j]; double xjp1 = x[j + 1];
		double tmp = xvar - x[j];
		double tmpA = x[j + 1] - xvar;

		double tmp3 = tmp * tmp * tmp;
		double tmp4 = tmpA * tmpA * tmpA;

		// S"
		double s2 = hjp1*(Mj*(xjp1 - xvar) + Mjp1*(xvar - xj));

		double A = (y[j + 1] - y[j]) / h[j + 1] - (h[j + 1] * (M[j + 1] - M[j])) / 6.0;
		double B = y[j] - (M[j] * h[j + 1] * h[j + 1]) / 6.0;

		// S'
		double s1 = -0.5*hjp1*Mj*tmpA*tmpA + 0.5*hjp1*Mjp1*tmp*tmp + A;

		// S
		double s0 = (M[j] * tmp4) / (6.0 * h[j + 1])
			+ (M[j + 1] * tmp3) / (6.0 * h[j + 1])
			+ (A * tmp)
			+ B;

		return std::make_tuple(s0, s1, s2);
	}

	double Derivative(double xvar)
	{ // Solve for derivative S'

		auto j = findAbscissa(x, xvar);

		double Mj = M[j]; double Mjp1 = M[j + 1];

		double hjp1 = 1.0 / h[j + 1];
		double xj = x[j]; double xjp1 = x[j + 1];
		double tmp = xvar - x[j];
		double tmpA = x[j + 1] - xvar;

		double tmp3 = tmp * tmp * tmp;
		double tmp4 = tmpA * tmpA * tmpA;

		double A = (y[j + 1] - y[j]) / h[j + 1] - (h[j + 1] * (M[j + 1] - M[j])) / 6.0;

		// S'
		double s1 = -0.5*hjp1*Mj*tmpA*tmpA + 0.5*hjp1*Mjp1*tmp*tmp + A;

		return s1;
	}

	double SecondDerivative(double xvar)
	{ // Solve for 2nd derivative S''

		auto j = findAbscissa(x, xvar);

		double Mj = M[j]; double Mjp1 = M[j + 1];

		double hjp1 = 1.0 / h[j + 1];
		double xj = x[j]; double xjp1 = x[j + 1];
	
		// S"
		double s2 = hjp1*(Mj*(xjp1 - xvar) + Mjp1*(xvar - xj));

		return s2;
	}

};

#endif