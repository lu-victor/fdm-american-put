#include <iostream>
#include <vector>

#include "DoubleSweep.hpp"

int main()
{
	using value_type = double;
	using Vector = std::vector<value_type>;

	std::size_t J = 20;
	std::cout << "Number of subdivisions J ";
	std::cin >> J;

	double h = 1.0 / static_cast<double>(J);

	double BCL = 0.0;
	double BCR = 0.0;

	Vector a(J + 1, 1.0);
	Vector b(J + 1, -2.0);
	Vector c(J + 1, 1.0);
	// right hand side for u'' = -2
	Vector r(J + 1, -2 * h*h);

	DoubleSweep<value_type> mySolver(a, b, c, r, BCL, BCR);
	Vector result = std::move(mySolver());

	auto exact = [](double x) { return x * (1.0 - x); };
	double val = 0;

	for (std::size_t j = 0; j < result.size() ; ++j)
	{ // The values should be zero
		std::cout << j << ", " << result[j] << ", " << exact(val)<< '\n';
			val += h;
	}
	return 0;

}