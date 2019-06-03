// DoubleSweep.hpp
// 
// The Godunov Balayage method for solving tridiagonal systems. We
// should be careful with start and end indices. The solution vector
// is in the range [0, j] as well as the vectors representing the 
// the main and of-diagonals of the tridiagonal system.
//
// The solution has size J+1 and contains the boundary conditions.
//
//	Note that this kind of problem originates typically from problems 
//  with Dirichlet boundary conditions. It is a competitor of the LU
//	decomposition for tridiagonal matrices.
//
// (C) Datasim Component Technology BV 2000-2017

#ifndef DOUBLESWEEP_HPP
#define DOUBLESWEEP_HPP

#include <vector>
#include <iostream>

// Using template template parameters to model vectors
template <typename T, 
			template <typename T, typename Alloc> class Container = std::vector,
			typename Alloc = std::allocator<T>>
class DoubleSweep
{ // The Balayage method from Godounov

private:
	
		// The vectors of length J and start index 0
		Container<T, Alloc> a, b, c, f;

		// Dirichlet boundary conditions
		T left;		// Left boundary condition
		T right;	// Right boundary condition

		// Work arrays
		Container<T, Alloc> L;
		Container<T, Alloc> K;

public:
		// Constructors and destructor
		DoubleSweep() = delete;
		DoubleSweep(const DoubleSweep<T,Container, Alloc>& s2) = delete;

		// Create members to initialise input for AU = F, A = (a,b,c)
		DoubleSweep(const Container<T, Alloc>& lowerDiagonal, const Container<T, Alloc>& diagonal,
			const Container<T, Alloc>& upperDiagonal, const Container<T, Alloc>& F,
			const T& bc_left, const T& bc_right)
		{

			// Vectors are copied
			a = lowerDiagonal;
			b = diagonal; 
			c = upperDiagonal;
			f = F;

			left = bc_left;
			right = bc_right;

			std::size_t N = a.size();

			// Work arrays
			L = Container<T, Alloc>(N, 0);
			K = Container<T, Alloc>(N, 0);
		}
		virtual ~DoubleSweep() = default;


		// Operator overloading
		DoubleSweep<T,Container, Alloc>& operator = (const DoubleSweep<T, Container>& i2) = delete;

		Container<T, Alloc> solve() 
		{ // Result; this is a vector in closed range [0, J], a vector of size J+1.
			std::size_t N = a.size();
		
			// equation 13.7	
			L[0] = 0.0;
			K[0] = left;

			// Equation 13.6
			std::size_t SZ = L.size();
			//quick fix not SZ but SZ-1 -> J-1
			for (std::size_t j = 1; j < SZ-1; ++j)
			{ // L

				double tmp = b[j] + (a[j] * L[j - 1]);

				L[j] = -c[j] / tmp;
				K[j] = (f[j] - (a[j] * K[j - 1])) / tmp;
			}

				
			// Equation 13.5. Recycle array f for u
			f[0] = left;
			f[N - 1] = right;
			for (std::size_t j = f.size() - 2; j >= 1; --j)
			{ // U

				f[j] = (L[j] * f[j + 1]) + K[j];
			}
		
			return f;
		}

		Container<T, Alloc> operator () ()
		{
			return solve();
		}

};



#endif

//VERIFIED