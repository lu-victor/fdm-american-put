#pragma once

#ifndef OPTION_HEADER
#define OPTION_HEADER

class Option {
protected:
	double Strike;
	double TimeToMaturity;
	double Sigma;
	double r;
public:
	Option(double K, double T, double sigma, double r);

	virtual double getExerciseValue(double s, double t);
	virtual double getBlackScholesValue(double s);

	double getBinomialTreeValue(double s, int N);

	virtual double getValue(double s);
};

class EuropeanOption : public Option {
public:
	EuropeanOption(double K, double T, double sigma, double rate) : Option(K, T, sigma, rate) {};
	virtual double getValue(double s);

};

class AmericanOption : public Option {
public:
	AmericanOption(double K, double T, double sigma, double rate) : Option(K, T, sigma, rate) {};
	virtual double getValue(double s);

};

class EuropeanCall : public EuropeanOption {
public:
	EuropeanCall(double K, double T, double sigma, double rate) : EuropeanOption(K, T, sigma, rate) {};
	virtual double getExerciseValue(double s, double t);
	virtual double getBlackScholesValue(double s);

};

class EuropeanPut : public EuropeanOption {
public:
	EuropeanPut(double K, double T, double sigma, double rate) : EuropeanOption(K, T, sigma, rate) {};
	virtual double getExerciseValue(double s, double t);
	virtual double getBlackScholesValue(double s);

};

class AmericanCall : public AmericanOption {
public:
	AmericanCall(double K, double T, double sigma, double rate) : AmericanOption(K, T, sigma, rate) {};
	virtual double getExerciseValue(double s, double t);

};

class AmericanPut : public AmericanOption {
public:
	AmericanPut(double K, double T, double sigma, double rate) : AmericanOption(K, T, sigma, rate) {};
	virtual double getExerciseValue(double s, double t);

};

#endif