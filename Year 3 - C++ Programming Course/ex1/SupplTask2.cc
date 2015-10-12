#include <iostream>
#include <cmath>
#include "cavlib/constants.hh"
#include "cavlib/gslint_integ.hh"

using namespace std;

// This program evaluates the C(u) and S(u) integrals
// in the form they are given in the manual. It uses the
// gslint_integ.hh class-based interface to the 
// gsl_integration_qag routine. The generate_amplitude()
// and generate_phase() procedures then use the integrals
// to estimate the relative amplitude and phase against the
// distance from the centre of the screen for a Fresnel
// diffraction from a slit.

double const pi = C::pi; // define pi

double lambda = 1, d = 10, D = 30;
double unit = d/sqrt(lambda*D); // define scaling unit

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


// Evaluates the realtive amplitude between x0 and x1
double rel_amplitude(double x0, double x1)
{
	double Cu_0, Cu_1, Su_0, Su_1, error;
	
	Cu_0 = integral_C(x0, &error);
	Cu_1 = integral_C(x1, &error);
	Su_0 = integral_S(x0, &error);
	Su_1 = integral_S(x1, &error);
	
	return sqrt((Cu_0 - Cu_1)*(Cu_0 - Cu_1) + (Su_0 - Su_1)*(Su_0 - Su_1));
}


// Evaluates the raltive phase between x0 and x1 [0; pi/2]
double phase(double x0, double x1)
{
	double Cu_0, Cu_1, Su_0, Su_1, error;
	
	Cu_0 = integral_C(x0, &error);
	Cu_1 = integral_C(x1, &error);
	Su_0 = integral_S(x0, &error);
	Su_1 = integral_S(x1, &error);
	
	return atan(fabs((Su_0 - Su_1)/(Cu_0 - Cu_1)));
}


// Prints data points for the amplitude against distance from centre
// x0 and x1 represent the two ends of the amplitude vector on the
// Cornu spiral. The integral length along the spiral between x0 and x1
// remains unchanged, while they are translated along the spiral.
// (x0+x1)/2 corresponds to the distance from the centre of the screen
// to the point of evaluation.
void generate_amplitude()
{
	double x0 = -unit/2, x1 = unit/2;
	
	for (int i=0; i<10000; i++)
	{
		cout<<(x0+x1)/2<<" "<<rel_amplitude(x0, x1)<<endl;
		cout<<-(x0+x1)/2<<" "<<rel_amplitude(x0, x1)<<endl; //symmetry
		x0 += unit*(1e-3);
		x1 += unit*(1e-3);
	}
}

// Prints data points for the phase against distance from centre
// x0 and x1 represent the two ends of the amplitude vector on the
// Cornu spiral. The integral length along the spiral between x0 and x1
// remains unchanged, while they are translated along the spiral.
// (x0+x1)/2 corresponds to the distance from the centre of the screen
// to the point of evaluation.
void generate_phase()
{
	double x0 = -unit/2, x1 = unit/2;
	
	for (int i=0; i<20000; i++)
	{
		cout<<(x0+x1)/2<<" "<<phase(x0, x1)<<endl;
		cout<<-(x0+x1)/2<<" "<<phase(x0, x1)<<endl; //symmetry
		x0 += unit*(1e-4);
		x1 += unit*(1e-4);
	}
}

int main()
{
	generate_amplitude();
	//generate_phase();
	return 0;
}
