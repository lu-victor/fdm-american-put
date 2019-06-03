#include "CubicSpline.hpp"
#include <iomanip>

int main() {
	std::vector<double> y = { 1,2,3 };
	std::vector<double> x = { 0.0, 1,2 };
	CubicSplineInterpolator csi(x, y, SecondDeriv);
	
	double xvar = 1.5;

	try
	{
		double result = csi.Solve(xvar);
		std::cout << "Interpolated value at " << xvar << ": " <<
			std::setprecision(16) << result << std::endl;
		
	}
	catch (std::exception& e)
	{ // Catch not in range values
		std::cout << e.what() << '\n';
	}
	return 0;
}