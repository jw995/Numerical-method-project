/*
 * 24-703 Group Project
 * Group Member: Jiayi Wang, Di Wu, Kai Ge
 * Andrew id   : jiayiw2, dwu1, kge
 *
 * Description:
 * Together with Integration1.cpp, this file implements the integration of
 * certain function needed in matrix K and P.
 */

#ifndef INTEGRATION1_H_IS_INCLUDED
#define INTEGRATION1_H_IS_INCLUDED
#include <functional>

double K(int i, int j, int n, double h);
double P(int i, int n, double h);

class Integrate
{
public:
	double integrate(double x0, double x1,
		double step, std::function <double(double)> func);
};

// k(nn) class
class k11
{
public:
	double h, xc;

	k11(double i, double HH);
	double Evaluate(double x) const;
	double run(k11 item);
};

class k12
{
public:
	double h, xc;

	k12(double i, double HH);
	double Evaluate(double x) const;
	double run(k12 item);
};

class k13
{
public:
	double h, xc;

	k13(double i, double HH);
	double Evaluate(double x) const;
	double run(k13 item);
};

class k22
{
public:
	double h, xc;

	k22(double i, double HH);
	double Evaluate(double x) const;
	double run(k22 item);
};

class k23
{
public:
	double h, xc;

	k23(double i, double HH);
	double Evaluate(double x) const;
	double run(k23 item);
};

class k33
{
public:
	double h, xc;

	k33(double i, double HH);
	double Evaluate(double x) const;
	double run(k33 item);
};

//p(n) class
class p1
{
public:
	double h, xc;

	p1(double i, double HH);
	double Evaluate(double x) const;
	double run(p1 item);
};

class p2
{
public:
	double h, xc;

	p2(double i, double HH);
	double Evaluate(double x) const;
	double run(p2 item);
};

class p3
{
public:
	double h, xc;

	p3(double i, double HH);
	double Evaluate(double x) const;
	double run(p3 item);
};
#endif
