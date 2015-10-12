#include <iostream>
#include <cmath>
#include "cavlib/constants.hh"
#include "cavlib/gslint_integ.hh"

using namespace std;

// This program evaluates the C(u) and S(u) integrals
// in the form they are given in the manual. It uses the
// gslint_integ.hh class-based interface to the 
// gsl_integration_qag routine. The generate_spiral()
// procedure then loops through different values of u to
// generate data points from which the Cornu spiral can
// be restored.

double const pi = C::pi; // define pi

// Define the integrand of the C(u) integral
double integrand_function1(double x, void* param)
{
	return cos(pi*x*x/2);
}


// Define the integrand of the S(u) integral
double integrand_function2(double x, void* param)
{
	return sin(pi*x*x/2);
}


// Evaluates the C(u) integral
double integral_C(double u, double *error)
{
	cav::Gsl_Integrator Integrator(&integrand_function1);
	
	return Integrator.integrate(0, u, error);
}


// Evaluates the S(u) integral
double integral_S(double u, double *error)
{
	cav::Gsl_Integrator Integrator(&integrand_function2);
	
	return Integrator.integrate(0, u, error);
}


// Prints data points, which are then used to plot the spiral
void generate_spiral()
{
	double error1 = 0, error2 = 0;
	
	for (int i=0; i<60000; i++)
	{
		double result1 = integral_C(i*(1e-4), &error1);
		double result2 = integral_S(i*(1e-4), &error2);
		cout<<result1<<" "<<result2<<endl;
	}
}

int main()
{
	generate_spiral();
	return 0;
}
