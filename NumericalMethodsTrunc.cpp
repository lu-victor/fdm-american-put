#include "NumericalMethodsTrunc.h"

double NumericalMethodsTrunc::explicitEuler() {


	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i < vecOld.size(); i++) {
		vecOld[i] = std::max(K - i * deltaS, 0.0);
	}
	std::vector<double> vecNew(M + 1);
	vecNew[0] = K;
	vecNew[M] = 0;

	std::vector<double> a(M);
	std::vector<double> b(M);
	std::vector<double> c(M);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = (-(r - q)*0.5*j + sigma * sigma*j*j*0.5)*deltaT;
		b[j] = -r * deltaT + 1 - sigma * sigma*j*j*deltaT;
		c[j] = (0.5*(r - q)*j + 0.5*sigma*sigma*j*j)*deltaT;
	}
	// we iterate through time until arriving at t=T
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(a[j] * vecOld[j - 1] + b[j] * vecOld[j] + c[j] * vecOld[j + 1], K - j * deltaS);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];
	CubicSplineInterpolator csi(space_mesh, vecOld, SecondDeriv);

	return csi.Solve(S_0);
	
}

std::vector<double> NumericalMethodsTrunc::explicitEulerVec() {


	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i < vecOld.size(); i++) {
		vecOld[i] = std::max(K - i * deltaS, 0.0);
	}
	std::vector<double> vecNew(M + 1);
	vecNew[0] = K;
	vecNew[M] = 0;

	std::vector<double> a(M);
	std::vector<double> b(M);
	std::vector<double> c(M);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = (-(r - q)*0.5*j + sigma * sigma*j*j*0.5)*deltaT;
		b[j] = -r * deltaT + 1 - sigma * sigma*j*j*deltaT;
		c[j] = (0.5*(r - q)*j + 0.5*sigma*sigma*j*j)*deltaT;
	}
	// we iterate through time until arriving at t=T
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(a[j] * vecOld[j - 1] + b[j] * vecOld[j] + c[j] * vecOld[j + 1], K - j * deltaS);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];
	//CubicSplineInterpolator csi(space_mesh, vecOld, SecondDeriv);

	//return csi.Solve(S_0);
	return vecOld;

}




double NumericalMethodsTrunc::implicitEuler() {

	std::vector<double> vecOld(M + 1);
	for (int i = 0; i < vecOld.size(); i++) {
		vecOld[i] = std::max(K - i * deltaS, 0.0);
	}
	std::vector<double> vecNew(M + 1);


	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = deltaT / (1 - r * deltaT) * ((r - q)*0.5 * j - 0.5*sigma*sigma*j*j);
		b[j] = 1 / (1 - r * deltaT) * (1 + sigma * sigma*j*j*deltaT);
		c[j] = deltaT / (1 - r * deltaT) * (-0.5*(r - q)*j - 0.5*sigma*sigma*j*j);
	}
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		DoubleSweep<double> mySolver(a, b, c, vecOld, BCL, BCR);
		std::vector<double> vecNew = std::move(mySolver());

		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(vecNew[j], K - j * deltaS);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	CubicSplineInterpolator csi(space_mesh, vecOld, SecondDeriv);

	return csi.Solve(S_0);
}

std::vector<double> NumericalMethodsTrunc::implicitEulerVec() {

	std::vector<double> vecOld(M + 1);
	for (int i = 0; i < vecOld.size(); i++) {
		vecOld[i] = std::max(K - i * deltaS, 0.0);
	}
	std::vector<double> vecNew(M + 1);


	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = deltaT / (1 - r * deltaT) * ((r - q)*0.5 * j - 0.5*sigma*sigma*j*j);
		b[j] = 1 / (1 - r * deltaT) * (1 + sigma * sigma*j*j*deltaT);
		c[j] = deltaT / (1 - r * deltaT) * (-0.5*(r - q)*j - 0.5*sigma*sigma*j*j);
	}
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		DoubleSweep<double> mySolver(a, b, c, vecOld, BCL, BCR);
		std::vector<double> vecNew = std::move(mySolver());

		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(vecNew[j], K - j * deltaS);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	return vecOld;
}
//Crank Nicolson is theta method with theta = 1/2
double NumericalMethodsTrunc::thetaMethod(double theta) {
	// number of space subdivisions
	// int M_prime = 10;
	//double deltaS = S_0 / static_cast<double>(M_prime);

	//int M = static_cast<int> ((3 * K) / deltaS);
	//double S_max = M * deltaS;

	// number of time subdivisions
	// int N = 100;
	//double deltaT = T / static_cast<double>(N);

	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i < vecOld.size(); i++) {
		vecOld[i] = std::max(K - i * deltaS, 0.0);
	}
	std::vector<double> vecNew(M + 1);
	//BC
	//double BCL = K;
	//double BCR = 0;

	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);
	std::vector<double> alpha(M + 1);
	std::vector<double> beta(M + 1);
	std::vector<double> gamma(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = theta * deltaT*(0.5*(r - q)*j - 0.5*sigma*sigma*j*j);
		b[j] = 1 + theta * deltaT*sigma*sigma*j*j;
		c[j] = theta * deltaT*(-0.5*(r - q)*j - 0.5*sigma*sigma*j*j);
		alpha[j] = (1 - theta)*deltaT*0.5*j*(-(r - q) + sigma * sigma*j);
		beta[j] = (1 - r * theta*deltaT - (1 - theta)*deltaT*(r + sigma * sigma*j*j));
		gamma[j] = (1 - theta)*deltaT*0.5*j*((r - q) + sigma * sigma*j);
	}

	std::vector<double> rightHandSide(M + 1);
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		for (int k = 1; k <= M - 1; k++) {
			rightHandSide[k] = alpha[k] * vecOld[k - 1] + beta[k] * vecOld[k] + gamma[k] * vecOld[k + 1];
		}
		DoubleSweep<double> mySolver(a, b, c, rightHandSide, BCL, BCR);
		std::vector<double> vecNew = std::move(mySolver());

		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(vecNew[j], K - j * deltaS);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	CubicSplineInterpolator csi(space_mesh, vecOld, SecondDeriv);

	return csi.Solve(S_0);
}

std::vector<double> NumericalMethodsTrunc::thetaMethodVec(double theta) {
	// number of space subdivisions
	// int M_prime = 10;
	//double deltaS = S_0 / static_cast<double>(M_prime);

	//int M = static_cast<int> ((3 * K) / deltaS);
	//double S_max = M * deltaS;

	// number of time subdivisions
	// int N = 100;
	//double deltaT = T / static_cast<double>(N);

	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i < vecOld.size(); i++) {
		vecOld[i] = std::max(K - i * deltaS, 0.0);
	}
	std::vector<double> vecNew(M + 1);
	//BC
	//double BCL = K;
	//double BCR = 0;

	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);
	std::vector<double> alpha(M + 1);
	std::vector<double> beta(M + 1);
	std::vector<double> gamma(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = theta * deltaT*(0.5*(r - q)*j - 0.5*sigma*sigma*j*j);
		b[j] = 1 + theta * deltaT*sigma*sigma*j*j;
		c[j] = theta * deltaT*(-0.5*(r - q)*j - 0.5*sigma*sigma*j*j);
		alpha[j] = (1 - theta)*deltaT*0.5*j*(-(r - q) + sigma * sigma*j);
		beta[j] = (1 - r * theta*deltaT - (1 - theta)*deltaT*(r + sigma * sigma*j*j));
		gamma[j] = (1 - theta)*deltaT*0.5*j*((r - q) + sigma * sigma*j);
	}

	std::vector<double> rightHandSide(M + 1);
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		for (int k = 1; k <= M - 1; k++) {
			rightHandSide[k] = alpha[k] * vecOld[k - 1] + beta[k] * vecOld[k] + gamma[k] * vecOld[k + 1];
		}
		DoubleSweep<double> mySolver(a, b, c, rightHandSide, BCL, BCR);
		std::vector<double> vecNew = std::move(mySolver());

		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(vecNew[j], K - j * deltaS);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	return vecOld;
}

double NumericalMethodsTrunc::ADE() {
	double current = deltaT;

	std::vector<double> U_old(M + 1);
	std::vector<double> U(M + 1); // L-R sweep
	std::vector<double> V_old(M + 1);
	std::vector<double> V(M + 1); // R-L sweep

	//std::vector<double> vecNew(M + 1); // vector at time i+1

	// vectors at t=0 are known
	for (int i = 0; i <= M; i++) {
		U_old[i] = V_old[i] = std::max(K - i * deltaS, 0.0);
	}
	double aj, bj, cj, dj, ej;
	while (current <= T) {
		// update at new time level i+1
		U[0] = BCL;
		U[M] = BCR;

		V[0] = BCL;
		V[M] = BCR;
		// up sweep
		for (int j = 1; j <= M - 1; j++) {
			aj = 0.5*(r - q)*j*deltaT;
			bj = r * deltaT - 1 + 0.5*sigma*sigma*j*j*deltaT;
			cj = -0.5*(r - q)*j*deltaT - 0.5*sigma*sigma*j*j*deltaT;
			dj = -0.5*sigma*sigma*j*j*deltaT;
			ej = -1.0 / (1.0 + 0.5*sigma*sigma*j*j*deltaT);

			U[j] = std::max(ej * (aj*U_old[j - 1] + bj * U_old[j] + cj * U_old[j + 1] + dj * U[j - 1]), K-j*deltaS);
		}
	
		U_old = U;

		//down sweep
		for (int j = M - 1; j >= 1; j--) {
			aj = 0.5*(r - q)*j*deltaT - 0.5*sigma*sigma*j*j*deltaT;
			bj = r * deltaT - 1 + 0.5*sigma*sigma*j*j*deltaT;
			cj = -0.5*(r - q)*j*deltaT;
			dj = -0.5*sigma*sigma*j*j*deltaT;
			ej = -1.0 / (1.0 + 0.5*sigma*sigma*j*j*deltaT);

			V[j] = std::max(ej * (aj*V_old[j - 1] + bj * V_old[j] + cj * V_old[j + 1] + dj * V[j + 1]), K-j*deltaS);
		}
		V_old = V;

		current += deltaT;
	}
	std::vector<double> vecNew(M + 1);
	vecNew[0] = K;
	vecNew[M] = 0;
	for (int j = 1; j <= M - 1; j++) {
		vecNew[j] = 0.5*(U[j] + V[j]);
	}

	CubicSplineInterpolator csi(space_mesh, vecNew, SecondDeriv);

	return csi.Solve(S_0);
	//return 0.5*(U[M_prime] + V[M_prime]);
}

std::vector<double> NumericalMethodsTrunc::ADEVec() {
	double current = deltaT;

	std::vector<double> U_old(M + 1);
	std::vector<double> U(M + 1); // L-R sweep
	std::vector<double> V_old(M + 1);
	std::vector<double> V(M + 1); // R-L sweep

	//std::vector<double> vecNew(M + 1); // vector at time i+1

	// vectors at t=0 are known
	for (int i = 0; i <= M; i++) {
		U_old[i] = V_old[i] = std::max(K - i * deltaS, 0.0);
	}
	double aj, bj, cj, dj, ej;
	while (current <= T) {
		// update at new time level i+1
		U[0] = BCL;
		U[M] = BCR;

		V[0] = BCL;
		V[M] = BCR;
		// up sweep
		for (int j = 1; j <= M - 1; j++) {
			aj = 0.5*(r - q)*j*deltaT;
			bj = r * deltaT - 1 + 0.5*sigma*sigma*j*j*deltaT;
			cj = -0.5*(r - q)*j*deltaT - 0.5*sigma*sigma*j*j*deltaT;
			dj = -0.5*sigma*sigma*j*j*deltaT;
			ej = -1.0 / (1.0 + 0.5*sigma*sigma*j*j*deltaT);

			U[j] = std::max(ej * (aj*U_old[j - 1] + bj * U_old[j] + cj * U_old[j + 1] + dj * U[j - 1]), K - j * deltaS);
		}

		U_old = U;

		//down sweep
		for (int j = M - 1; j >= 1; j--) {
			aj = 0.5*(r - q)*j*deltaT - 0.5*sigma*sigma*j*j*deltaT;
			bj = r * deltaT - 1 + 0.5*sigma*sigma*j*j*deltaT;
			cj = -0.5*(r - q)*j*deltaT;
			dj = -0.5*sigma*sigma*j*j*deltaT;
			ej = -1.0 / (1.0 + 0.5*sigma*sigma*j*j*deltaT);

			V[j] = std::max(ej * (aj*V_old[j - 1] + bj * V_old[j] + cj * V_old[j + 1] + dj * V[j + 1]), K - j * deltaS);
		}
		V_old = V;

		current += deltaT;
	}
	std::vector<double> vecNew(M + 1);
	vecNew[0] = K;
	vecNew[M] = 0;
	for (int j = 1; j <= M - 1; j++) {
		vecNew[j] = 0.5*(U[j] + V[j]);
	}

	return vecNew;
}

void NumericalMethodsTrunc::printValues() {
	std::cout << "Current stock price: " << S_0 << ", Strike: " << K << ", Maturity (years): " << T << std::endl;
	std::cout << "Risk-free rate: " << r << ", Dividends: " << q << ", Volatility: " << sigma << std::endl;

	std::cout << "Space subdivisions: " << M << ", time subdivisions: " << N << std::endl;
	std::cout << "Explicit Euler Trunc: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->explicitEuler() << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
	std::cout << "Implicit Euler Trunc: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->implicitEuler() << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
	std::cout << "Crank-Nicolson Trunc: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->thetaMethod(0.5) << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
	std::cout << "ADE Trunc: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->ADE() << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
}

