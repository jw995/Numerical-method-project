/*
 * 24-703 Group Project
 * Group Member: Jiayi Wang, Di Wu, Kai Ge
 * Andrew id   : jiayiw2, dwu1, kge
 * 
 * Description:
 * Together with Integration1.h, this file implements the integration of 
 * certain function needed in matrix K and P. 
 */


#include <stdio.h>
#include <functional>
#include <math.h>
#include "Integration1.h"

double K(int i, int j, int n, double h) {
	double result = 0.;
	if (i == 1) {
		if (j == 1) {
			k11 item11(n, h);
			result = h*item11.run(item11);
			return result;
		}
		if (j == 2) {
			k12 item12(n, h);
			result = h*item12.run(item12);
			return result;
		}
		if (j == 3) {
			k13 item13(n, h);
			result = h*item13.run(item13);
			return result;
		}
	}
	if (i == 2) {
		if (j == 1) {
			k12 item12(n, h);
			result = h*item12.run(item12);
			return result;
		}
		if (j == 2) {
			k22 item22(n, h);
			result = h*item22.run(item22);
			return result;
		}
		if (j == 3) {
			k23 item23(n, h);
			result = h*item23.run(item23);
			return result;
		}
	}
	if (i == 3) {
		if (j == 1) {
			k13 item13(n, h);
			result = h*item13.run(item13);
			return result;
		}
		if (j == 2) {
			k23 item23(n, h);
			result = h*item23.run(item23);
			return result;
		}
		if (j == 3) {
			k33 item33(n, h);
			result = h*item33.run(item33);
			return result;
		}
	}
}

double P(int i, int n, double h) {
	double result = 0.;
	if (i == 1) {
		p1 item1(n, h);
		result = h*item1.run(item1);
		return result;
	}
	if (i == 2) {
		p2 item2(n, h);
		result = h*item2.run(item2);
		return result;
	}
	if (i == 3) {
		p3 item3(n, h);
		result = h*item3.run(item3);
		return result;
	}
	return result;
}




/////////////////////////////////////////////////////////////

double Integrate::integrate(double x0, double x1,
	double step, std::function <double(double)> func)
{
	double sum = 0.0;
	for (double x = x0; x <= x1; x += step)
	{
		sum += func(x)*step;
	}
	return sum;
}

/////////////////////////////////////////////////////////////
// k(nn) functions
k11::k11(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */double k11::Evaluate(double x) const
{
	double a = xc + h*x;
	double b = 1 / h*(x - 0.5);
	double fun = pow(a, 0.5)*pow(b, 2);
	
	return fun;
}

double k11::run(k11 item)
{
	Integrate itg;
	std::function <double(double)> func11 =
		std::bind(&k11::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func11);
}


k12::k12(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */ double k12::Evaluate(double x) const
{
	double a = xc + h*x;
	double b = 1 / h*(x - 0.5);
	double fun = pow(a, 0.5)*b*(-2 * x / h);

	return fun;
}

double k12::run(k12 item)
{
	Integrate itg;
	std::function <double(double)> func12 =
		std::bind(&k12::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func12);
}



k13::k13(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */ double k13::Evaluate(double x) const
{
	double a = xc + h*x;
	double b = 1 / h*(x - 0.5);
	double c = 1 / h*(x + 0.5);
	double fun = pow(a, 0.5)*b*c;

	return fun;
}

double k13::run(k13 item)
{
	Integrate itg;
	std::function <double(double)> func13 =
		std::bind(&k13::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func13);
}



k22::k22(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */ double k22::Evaluate(double x) const
{
	double a = xc + h*x;
	double fun = pow(a, 0.5)*pow(-2*x/h,2);

	return fun;
}

double k22::run(k22 item)
{
	Integrate itg;
	std::function <double(double)> func22 =
		std::bind(&k22::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func22);
}




k23::k23(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */ double k23::Evaluate(double x) const
{
	double a = xc + h*x;
	double c = 1 / h*(x + 0.5);
	double fun = pow(a, 0.5)*c*(-2 * x / h);

	return fun;
}

double k23::run(k23 item)
{
	Integrate itg;
	std::function <double(double)> func23 =
		std::bind(&k23::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func23);
}


k33::k33(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */ double k33::Evaluate(double x) const
{
	double a = xc + h*x;
	double c = 1 / h*(x + 0.5);
	double fun = pow(a, 0.5)*pow(c, 2);

	return fun;
}

double k33::run(k33 item)
{
	Integrate itg;
	std::function <double(double)> func33 =
		std::bind(&k33::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func33);
}

//p(n) functions
p1::p1(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */double p1::Evaluate(double x) const
{
	double a = xc + h*x;
	double b = 0.5 * x *(x - 1);
	double fun = pow(a, 2) * b;

	return fun;
}
double p1::run(p1 item)
{
	Integrate itg;
	std::function <double(double)> func1 =
		std::bind(&p1::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func1);
}

p2::p2(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */double p2::Evaluate(double x) const
{
	double a = xc + h*x;
	double b = (1 - x)*(1 + x);
	double fun = pow(a, 2) * b;

	return fun;
}
double p2::run(p2 item)
{
	Integrate itg;
	std::function <double(double)> func2 =
		std::bind(&p2::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func2);
}

p3::p3(double i, double HH)
{
	h = HH;
	xc = 2 * i*h - h;
}
/* */double p3::Evaluate(double x) const
{
	double a = xc + h*x;
	double b = 0.5 * x *(x + 1);
	double fun = pow(a, 2) * b;

	return fun;
}
double p3::run(p3 item)
{
	Integrate itg;
	std::function <double(double)> func3 =
		std::bind(&p3::Evaluate, &item, std::placeholders::_1);
	return itg.integrate(-1.0, 1.0, 0.0001, func3);
}