// LUSolver.hpp
//
// Classes and functions for solving linear systems of equations 
// (numerical linear algebra).
//
// Modification Dates:
//
//	DD 2003-1-14 First code (tridiagonal)
///	DD 2012-9-25 bug fixed
//	2017-2-13 DD std::Vector<T>.
//
// (C) Datasim Education BV 2003-2017

#ifndef LUSolver_HPP
#define LUSolver_HPP

#include <vector>

template <typename T>
	using Vector = std::vector<T>;

template <class T> class LUTridiagonalSolver
{ // Solve tridiagonal matrix equation
private:

	// Defining arrays (input)
	// V2 optimise so to work with pointers
	Vector<T> a;	// The lower-diagonal array [0..J]
	Vector<T> b;	// The diagonal array [0..J] "baseline array"
	Vector<T> c;	// The upper-diagonal array [0..J]
	Vector<T> r;	// The right-hand side of the equation Au = r [0..J]

					// Work arrays

					// Coefficients of Lower and Upper matrices: A = LU
					// V2 use of Templated static Vector<T>s, but we must be careful
	Vector<T> beta;	// Range [0..J]
	Vector<T> gamma;// Range [0..J]
					// Solutions of temporary and final problems
	Vector<T> z;	// Range [0..J]

	std::size_t Size;

	void calculateBetaGamma_ZU(Vector<T>& r)
	{

		beta[0] = b[0];
		gamma[0] = c[0] / beta[0];

		for (std::size_t j = 1; j < Size - 1; ++j)
		{
			beta[j] = b[j] - (a[j] * gamma[j - 1]);
			gamma[j] = c[j] / beta[j];

		}

		beta[Size - 1] = b[Size - 1] - (a[Size - 1] * gamma[Size - 2]);

		// Calculate z and u

		// Forward direction
		z[0] = r[0] / beta[0];


		for (std::size_t j = 1; j < Size; ++j)
		{
			z[j] = (r[j] - (a[j] * z[j - 1])) / beta[j];

		}

		// Backward direction
		r[Size - 1] = z[Size - 1];

		for (long i = Size - 2; i >= 0; --i)
		{
			r[i] = z[i] - (gamma[i] * r[i + 1]);
		}

	}

public:
	LUTridiagonalSolver() = delete;
	LUTridiagonalSolver(const LUTridiagonalSolver<T>& source) = delete;
	virtual ~LUTridiagonalSolver() = default;
	LUTridiagonalSolver<T>& operator = (const LUTridiagonalSolver<T>& source) = delete;

LUTridiagonalSolver(Vector<T>& lower, Vector<T>& diagonal, Vector<T>& upper, Vector<T>& RHS)
{

	a = lower;
	b = diagonal;
	c = upper;
	r = RHS;

	Size = diagonal.size();

	beta = Vector<T>(Size);
	gamma = Vector<T>(Size);

	z = Vector<T>(Size);
}

	
Vector<T> solve()
{
		calculateBetaGamma_ZU(r);		// Calculate beta and gamma

		return r;
}

Vector<T> operator () ()
{
	return solve();
}

bool diagonallyDominant() const
{
		if (std::abs(b[0]) < std::abs(c[0]))
			return false;

		if (std::abs(b[Size - 1]) < std::abs(a[Size - 1]))
			return false;

		for (std::size_t j = 1; j < Size - 1; ++j)
		{
			if (std::abs(b[j]) < std::abs(a[j]) + std::abs(c[j]))
				return false;
		}

		return true;
}
};

#endif