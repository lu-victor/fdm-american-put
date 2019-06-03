// Utilities.hpp
//
// Handy functions to help with numerical linear algebra.
//
// (C) Datasim Education BV 2017
//

#ifndef Utilities_HPP
#define Utilities_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <boost/lexical_cast.hpp>

inline void print(const std::vector<double>& v)
{
	for (auto e : v)
	{
		std::cout << e << ",";
	}
	std::cout << '\n';
}

inline double LInfinityNorm(const std::vector<double>& v1, const std::vector<double>& v2)
{ // Max of absolute value of elements of v1 - v2

	double result = std::abs(v1[0] - v2[0]);

	for (std::size_t j = 1; j < v1.size(); ++j)
	{
		result = std::max<double>(result, std::abs(v1[j] - v2[j]));
	}

	return result;
}

inline std::tuple<double, std::size_t> HotSpotError(const std::vector<double>& v1, const std::vector<double>& v2)
{ // Max of absolute value of elements of v1 - v2 and identify *where* it occurs

	double result = std::abs(v1[0] - v2[0]);
	std::size_t index = 0;

	for (std::size_t j = 1; j < v1.size(); ++j)
	{
		double tmp = std::max<double>(result, std::abs(v1[j] - v2[j]));

		if (result < tmp)
		{ // Find the max difference

			result = tmp;
			index = j;
		}
	}

	return std::make_tuple(result, index);
}

inline std::size_t findAbscissa(const std::vector<double>& x, double xvar) 
{ // Will give index of LHS value <= xvar. 
  
//	std::cout << xvar; int yy; std::cin >> yy;
	if (xvar < x[0] || xvar >  x[x.size() - 1])
	{
		std::string s = "\nValue " + boost::lexical_cast<std::string>(xvar) + " not in range "
			+ "(" + boost::lexical_cast<std::string>(x[0]) + ","
			+ boost::lexical_cast<std::string>(x[x.size() - 1]) + ")";
		std::cout << "Outside range\n";
		throw std::out_of_range(s);
	}

	auto posA = std::lower_bound(std::begin(x), std::end(x), xvar); // Log complexity
												
	std::size_t index = std::distance(std::begin(x), posA)-1;

	return index;
}


inline std::vector<double> CreateMesh(const std::vector<double>& x)
{ // Create a refined mesh from a coarser mesh

	std::size_t N = x.size();
	std::vector<double> result(2 * N - 1);

	for (std::size_t j = 0; j < result.size(); j += 2)
	{ // Even 

		result[j] = x[j/2];
	}

	for (std::size_t j = 1; j < result.size(); j += 2)
	{ // Odd

		result[j] = 0.5*(result[j - 1] + result[j + 1]);
	}

	return result;
}

inline std::vector<double> CreateMesh(std::size_t n, double a, double b)
{ 

	std::vector<double> result(n+1);
	result[0] = a;
	result[result.size()-1] = b;

	double h = (b - a) / static_cast<double>(n);
	for (std::size_t j = 1; j < result.size()-1; ++j)
	{ 
		result[j] = result[j - 1] + h;
	}

	return result;
}




/*
std::vector<double> CreateDiscreteFunction(std::size_t n, double a, double b, 
											const std::function<double (double)>& f)
{ // Create a discrete function from a continuous function y = f(x)

	std::vector<double> y(n + 1);

	double h = (b - a) / static_cast<double>(n);
	double x = a;
	for (std::size_t j = 0; j < y.size(); ++j)
	{
		y[j] = f(x);
		x += h;
	}

	return y;
}
*/
inline std::vector<double> CreateDiscreteFunction(const std::vector<double>& x, const std::function<double(double)>& f)
{ // Create a discrete function from a continuous function y = f(x)

	std::vector<double> y(x.size());
	
	for (std::size_t j = 0; j < y.size(); ++j)
	{
		y[j] = f(x[j]);
	}

	return y;
}

#endif
