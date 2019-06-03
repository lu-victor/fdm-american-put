#include "option.h"
#include <cmath>
#include <vector>

Option::Option(double K, double T, double sigma, double rate){
	Strike = K;
	TimeToMaturity = T;
	Sigma = sigma;
	r = rate;
}

double Option::getExerciseValue(double s, double t)
{
	return 0.0;
}

double Option::getBlackScholesValue(double s)
{
	return 0.0;
}

double Option::getValue(double s)
{
	return 0.0;
}

double EuropeanCall::getExerciseValue(double s, double t) {
	if (t != TimeToMaturity) {
		return 0;
	}
	else {
		return (s - Strike > 0) ? (s - Strike) : 0;
	}
}
double EuropeanPut::getExerciseValue(double s, double t) {
	if (t != TimeToMaturity) {
		return 0;
	}
	else {
		return (Strike-s > 0) ? (Strike-s) : 0;
	}
}
double AmericanCall::getExerciseValue(double s, double t) {
	if (t <= TimeToMaturity) {
		return (s - Strike > 0) ? (s - Strike) : 0;
	}
	else {
		return 0;
	}
}
double AmericanPut::getExerciseValue(double s, double t) {
	if (t <= TimeToMaturity) {
		return (Strike - s > 0) ? (Strike - s) : 0;
	}
	else {
		return 0;
	}
}

double normalCDF(double x) 
{
	return std::erfc(-x / std::sqrt(2)) / 2;
}

double EuropeanCall::getBlackScholesValue(double s) {
	double d1 = 1 / (Sigma*sqrt(TimeToMaturity)) * (log(s / Strike) + (r + Sigma * Sigma / 2)*TimeToMaturity);
	double d2 = d1 - Sigma * sqrt(TimeToMaturity);

	return normalCDF(d1)*s - normalCDF(d2)*Strike*exp(-r * TimeToMaturity);
}

double EuropeanPut::getBlackScholesValue(double s) {
	double d1 = 1 / (Sigma*sqrt(TimeToMaturity)) * (log(s / Strike) + (r + Sigma * Sigma / 2)*TimeToMaturity);
	double d2 = d1 - Sigma * sqrt(TimeToMaturity);

	return normalCDF(-d2)*Strike*exp(-r * TimeToMaturity)-normalCDF(-d1)*s;
}

double Option::getBinomialTreeValue(double s, int N) {
	double deltaT = TimeToMaturity / N;
	double up = exp(Sigma*sqrt(deltaT));

	double p0 = (up - exp(-r * deltaT)) / (up*up - 1);
	double p1 = exp(-r * deltaT) - p0;

	std::vector<double> p(N+1);

	double spotAtTime;

	for (int i = 0; i <= N; i++) {
		spotAtTime = s * pow(up, 2 * i - N);
		p[i] = getExerciseValue(spotAtTime, TimeToMaturity);
	}

	double time = TimeToMaturity;
	double exercise;
	double presentValue;
	for (int j = N - 1; j >= 0; j--) {
		time = time - deltaT;
		for (int i = 0; i <= j; i++) {
			spotAtTime = s * pow(up, 2 * i - j);
			exercise = getExerciseValue(spotAtTime, time);
			presentValue = p0 * p[i + 1] + p1 * p[i];

			p[i] = (exercise > presentValue) ? exercise : presentValue;
		}
	}
	return p[0];
}

double EuropeanOption::getValue(double s) {
	return getBlackScholesValue(s);
}

double AmericanOption::getValue(double s) {
	return getBinomialTreeValue(s,1000);
}