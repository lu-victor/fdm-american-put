#include "NumericalMethods.h"

double NumericalMethods::explicitEuler() {

	// init 
	std::vector<double> vecOld(M + 1);

	// payoff known at t=0
	for (int i = 1; i <= M-1; i++) {
		vecOld[i] = std::max(K - S[i], 0.0);
	}
	vecOld[0] = K;
	vecOld[M] = 0;
	// value f at the next time step
	std::vector<double> vecNew(M + 1);
	vecNew[0] = K;
	vecNew[M] = 0;

	std::vector<double> a(M);
	std::vector<double> b(M);
	std::vector<double> c(M);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = deltaT * (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		b[j] = deltaT * (1.0 / deltaT - r - sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		c[j] = deltaT * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
	}
	// we iterate through time until arriving at t=T
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(a[j] * vecOld[j - 1] + b[j] * vecOld[j] + c[j] * vecOld[j + 1], K - S[j]);
		}
		vecOld = vecNew;
	}

	//for (int i = 0; i <= M; i++) {
	//	std::cout << i << ", space = " << space_mesh_Y[i] << ", vecOld = " << vecOld[i] << std::endl;
	//}

	CubicSplineInterpolator csi(space_mesh_Y, vecOld, SecondDeriv);
	
	return csi.Solve(y_0);
}

double NumericalMethods::implicitEuler() {

	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i <= M-1; i++) {
		vecOld[i] = std::max(K - S[i], 0.0);
	}
	vecOld[M] = 0;

	std::vector<double> vecNew(M + 1);
	//BC
	//double BCL = K;
	//double BCR = 0;

	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = 1.0/(r-1.0/deltaT)* (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		b[j] = 1.0 / (r - 1.0 / deltaT) * (-1.0 / deltaT  - sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		c[j] = 1.0 / (r - 1.0 / deltaT) * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
	}
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		DoubleSweep<double> mySolver(a, b, c, vecOld, BCL, BCR);
		std::vector<double> vecNew = std::move(mySolver());

		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(vecNew[j], K - S[j]);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	CubicSplineInterpolator csi(space_mesh_Y, vecOld, SecondDeriv);

	return csi.Solve(y_0);
}

//Crank Nicolson is theta method with theta = 1/2
double NumericalMethods::thetaMethod(double theta) {

	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i <= M-1; i++) {
		vecOld[i] = std::max(K - S[i], 0.0);
	}
	vecOld[M] = 0;
	std::vector<double> vecNew(M + 1);


	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);
	std::vector<double> alpha(M + 1);
	std::vector<double> beta(M + 1);
	std::vector<double> gamma(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = theta*(0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) - 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		b[j] = (1.0 / deltaT + theta*sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		c[j] = theta*(-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) - 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		alpha[j] = (1-theta) * (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		beta[j] = (1.0 / deltaT -r - (1-theta) * sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		gamma[j] = (1-theta) * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
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
			vecNew[j] = std::max(vecNew[j], K-S[j]);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	CubicSplineInterpolator csi(space_mesh_Y, vecOld, SecondDeriv);

	return csi.Solve(y_0);
}

double NumericalMethods::ADE() {
	double current = deltaT;

	std::vector<double> U_old(M + 1);
	std::vector<double> U(M + 1); // L-R sweep
	std::vector<double> V_old(M + 1);
	std::vector<double> V(M + 1); // R-L sweep

	//std::vector<double> vecNew(M + 1); // vector at time i+1

	// vectors at t=0 are known
	for (int i = 0; i <= M-1; i++) {
		U_old[i] = V_old[i] = std::max(K - S[i], 0.0);
	}
	U_old[M] = V_old[M] = 0;

	std::vector<double> a(M), b(M), c(M), d(M), e(M);
	for (int j = 1; j <= M - 1; j++) {
		a[j] = 0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY);
		b[j] = -1.0 / deltaT + r + 0.5* sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY);
		c[j] = -0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) - 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY);
		d[j] = -0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY);
		e[j] = -1.0 / (1.0 / deltaT + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
	}
	while (current <= T) {
		// update at new time level i+1
		U[0] = BCL;
		U[M] = BCR;

		V[0] = BCL;
		V[M] = BCR;
		// up sweep
		for (int j = 1; j <= M - 1; j++) {
			U[j] = std::max(e[j] * (a[j]*U_old[j - 1] + b[j] * U_old[j] + c[j] * U_old[j + 1] + d[j] * U[j - 1]), K - S[j]);
		}

		U_old = U;

		//down sweep
		// the new coeffs can be deduced from the up sweep
		for (int j = M - 1; j >= 1; j--) {

			V[j] = std::max(e[j] * ((a[j]+d[j])*U_old[j - 1] + b[j] * U_old[j] -a[j] * U_old[j + 1] + d[j] * U[j - 1]), K - S[j]);
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

	CubicSplineInterpolator csi(space_mesh_Y, vecNew, SecondDeriv);

	return csi.Solve(y_0);

}

std::vector<double> NumericalMethods::explicitEulerVec() {

	// init 
	std::vector<double> vecOld(M + 1);

	// payoff known at t=0
	for (int i = 1; i <= M - 1; i++) {
		vecOld[i] = std::max(K - S[i], 0.0);
	}
	vecOld[0] = K;
	vecOld[M] = 0;
	// value f at the next time step
	std::vector<double> vecNew(M + 1);
	vecNew[0] = K;
	vecNew[M] = 0;

	std::vector<double> a(M);
	std::vector<double> b(M);
	std::vector<double> c(M);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = deltaT * (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		b[j] = deltaT * (1.0 / deltaT - r - sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		c[j] = deltaT * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
	}
	// we iterate through time until arriving at t=T
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(a[j] * vecOld[j - 1] + b[j] * vecOld[j] + c[j] * vecOld[j + 1], K - S[j]);
		}
		vecOld = vecNew;
	}

	//for (int i = 0; i <= M; i++) {
	//	std::cout << i << ", space = " << space_mesh_Y[i] << ", vecOld = " << vecOld[i] << std::endl;
	//}

	return vecOld;
}

std::vector<double> NumericalMethods::implicitEulerVec() {

	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i <= M - 1; i++) {
		vecOld[i] = std::max(K - S[i], 0.0);
	}
	vecOld[M] = 0;

	std::vector<double> vecNew(M + 1);
	//BC
	//double BCL = K;
	//double BCR = 0;

	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = 1.0 / (r - 1.0 / deltaT)* (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		b[j] = 1.0 / (r - 1.0 / deltaT) * (-1.0 / deltaT - sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		c[j] = 1.0 / (r - 1.0 / deltaT) * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
	}
	for (int i = 1; i <= N; i++) {
		// computes the value of option at time t* = i*deltaT
		DoubleSweep<double> mySolver(a, b, c, vecOld, BCL, BCR);
		std::vector<double> vecNew = std::move(mySolver());

		for (int j = 1; j <= M - 1; j++) {
			vecNew[j] = std::max(vecNew[j], K - S[j]);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	return vecOld;
}

//Crank Nicolson is theta method with theta = 1/2
std::vector<double> NumericalMethods::thetaMethodVec(double theta) {

	// init 
	std::vector<double> vecOld(M + 1);
	for (int i = 0; i <= M - 1; i++) {
		vecOld[i] = std::max(K - S[i], 0.0);
	}
	vecOld[M] = 0;
	std::vector<double> vecNew(M + 1);


	std::vector<double> a(M + 1);
	std::vector<double> b(M + 1);
	std::vector<double> c(M + 1);
	std::vector<double> alpha(M + 1);
	std::vector<double> beta(M + 1);
	std::vector<double> gamma(M + 1);

	for (int j = 1; j <= M - 1; j++) {
		a[j] = theta * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) - 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		b[j] = (1.0 / deltaT + theta * sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		c[j] = theta * (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) - 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		alpha[j] = (1 - theta) * (-0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		beta[j] = (1.0 / deltaT - r - (1 - theta) * sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
		gamma[j] = (1 - theta) * (0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
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
			vecNew[j] = std::max(vecNew[j], K - S[j]);
		}
		vecOld = vecNew;
	}
	// actually with this procedure we get all the values at time t=0 for a range of spot prices
	// we return only the value of the put for S_0
	// return vecOld[M_prime];

	return vecOld;
}

std::vector<double> NumericalMethods::ADEVec() {
	double current = deltaT;

	std::vector<double> U_old(M + 1);
	std::vector<double> U(M + 1); // L-R sweep
	std::vector<double> V_old(M + 1);
	std::vector<double> V(M + 1); // R-L sweep

	//std::vector<double> vecNew(M + 1); // vector at time i+1

	// vectors at t=0 are known
	for (int i = 0; i <= M - 1; i++) {
		U_old[i] = V_old[i] = std::max(K - S[i], 0.0);
	}
	U_old[M] = V_old[M] = 0;

	std::vector<double> a(M), b(M), c(M), d(M), e(M);
	for (int j = 1; j <= M - 1; j++) {
		a[j] = 0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY);
		b[j] = -1.0 / deltaT + r + 0.5* sigma * sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY);
		c[j] = -0.5*j*(1 - j * deltaY)*(r - sigma * sigma*j*deltaY) - 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY);
		d[j] = -0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY);
		e[j] = -1.0 / (1.0 / deltaT + 0.5*sigma*sigma*j*j*(1 - j * deltaY)*(1 - j * deltaY));
	}
	while (current <= T) {
		// update at new time level i+1
		U[0] = BCL;
		U[M] = BCR;

		V[0] = BCL;
		V[M] = BCR;
		// up sweep
		for (int j = 1; j <= M - 1; j++) {
			U[j] = std::max(e[j] * (a[j] * U_old[j - 1] + b[j] * U_old[j] + c[j] * U_old[j + 1] + d[j] * U[j - 1]), K - S[j]);
		}

		U_old = U;

		//down sweep
		// the new coeffs can be deduced from the up sweep
		for (int j = M - 1; j >= 1; j--) {

			V[j] = std::max(e[j] * ((a[j] + d[j])*U_old[j - 1] + b[j] * U_old[j] - a[j] * U_old[j + 1] + d[j] * U[j - 1]), K - S[j]);
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

void NumericalMethods::printValues() {
	std::cout << "Current stock price: " << S_0 << ", Strike: " << K << ", Maturity (years): " << T << std::endl;
	std::cout << "Risk-free rate: " << r << ", Dividends: " << q << ", Volatility: " << sigma << std::endl;

	std::cout << "Space subdivisions: " << M << ", time subdivisions: " << N << std::endl;
	std::cout << "Explicit Euler: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->explicitEuler() << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
	std::cout << "Implicit Euler: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->implicitEuler() << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
	std::cout << "Crank-Nicolson: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->thetaMethod(0.5) << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;

	std::cout << "ADE: " << std::endl;
	std::cout << "Put value: " << std::setprecision(16) << this->ADE() << std::endl;
	std::cout << "----------------------------------------------------------" << std::endl;
}

