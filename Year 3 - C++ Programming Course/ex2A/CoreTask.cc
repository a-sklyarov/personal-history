#include <iostream>
#include <cmath>
#include <vector>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv.h"
#include "cavlib/string_util.hh"
#include "cavlib/constants.hh"

using namespace std;

double const pi = C::pi;


// Function to evaluate the derivatives.
// We work with y[0] = theta, y[1] = d(theta)/dt

int calc_derivs(double t, const double y[], double f[], void *params)
{
	// Extract the parameters from *param. In this case we have two
	// parameters - q and F, which are passed by an array
	
	double q=0, F=0; 
	q = ((double *)params)[0]; // extract q
	F = ((double *)params)[1]; // extract F
	
	// Describe the two coupled first order ODEs
	
	f[0] = y[1];
	f[1] = -sin(y[0]) - q*y[1] + F*sin(2*t/3);
	
	return GSL_SUCCESS;
}


// Function that evaluates the energy for given
// deflection theta and angular velocity omega.
// A factor of 1/2 is left in front of omega^2
// after g/l is taken to be unity.

double calc_energy(double theta, double omega)
{
	double energy = omega*omega/2 - cos(theta);
	return energy;
}


// Function that prints the archive, i.e. the values
// of a two dimensional vector (or a vector of vectors).
// The archive represents the output from the integration,
// having 4 columns (t, theta, omega, energy).

void print_archive(vector< vector<double> > archive)
{
	for (unsigned int i=0; i<archive.size(); i++)
	{
		for (unsigned int j=0; j<(archive.at(i)).size(); j++)
			cout<<(archive.at(i)).at(j)<<" ";
		
		cout<<endl;
	}
}


// Function that calculates the mean of a set of values,
// which are passed as a vector. It is used to calculate
// the mean period of the oscilations in the function 
// estimate_period().

double calc_mean(vector<double> values)
{
	double sum = 0;
	
	for (unsigned int i = 0; i < values.size(); i++) sum += values.at(i);
	
	return sum/(values.size());
}


// Function which estimates the mean period of the oscillations.
// It takes as an argument the archive (i.e. the output from the 
// integration), which has the matrix representation of many rows
// and 4 columns. The function then builds up many estimates of
// the period by measuring the time between the points where the
// deflection of theta changes sign from plus to minus. The returned 
// value is the mean value of all the individual estimtes.

double estimate_period(vector< vector<double> > archive)
{
	// container of the individual period estimates
	vector<double> period_values;
	
	// t1 and t2 are the markers whose difference gives the 
	// period. t is used as a help variable, theta and theta_next
	// are used to take the values of the deflections theta of two 
	// consecutive iterations. 
		
	double t1 = 0, t2 = 0, t, theta, theta_next;
	
	// Loop through the output of the integration by taking two
	// neighbouring values of theta at each step. If the first value 
	// is positive and the next one is negative, then we put one
	// of the markers (t1 or t2) at the time when that happens. When 
	// the next such instance occurs, we put our second marker at
	// that time. Every time we put a marker (except the very first
	// time), we take the difference between them and collect the result
	// as an estimate of the period.
	
	for (unsigned int i=0; i<(archive.size()-1); i++)
	{
		theta = (archive.at(i)).at(1); // deflection at row i
		theta_next = (archive.at(i+1)).at(1); // deflection at row i+1
		t = (archive.at(i)).at(0); // the time at row i
		
		if ((theta > 0) && (theta_next < 0))
			{					//decide which marker to move
				if (t1 > t2) 
				{
					t2 = t;
					period_values.push_back(fabs(t2 - t1));
				}
				else if (t1 < t2)
					{
						t1 = t;
						period_values.push_back(fabs(t2 - t1));
					}
					else t1 = t; //if t1 = t2, i.e. at the start
			}
	}
	
	return calc_mean(period_values); // return the mean of all values
}
	

// Function to do the actual integration using RK4 method. 

void integrate(double q, double F, double h, double err, double tmax, 
			   double theta0, double omega0)
{
	const int n_eqs = 2; // number of equations
	double params[2] = {q, F}; // parameters to be passed to calc_derivs
	double y[2] = {theta0, omega0}; // initial values of theta and omega
	double t = 0; // time starts from 0
	int n_steps = 0; // counter for the integration steps used to 
	                 // limit the output to only once in a given number
	                 // of steps
	
	vector< vector<double> > archive; // for the output
	
	// RK4 method to be used
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
	
	// Create a stepping function
	gsl_odeiv_step *step = gsl_odeiv_step_alloc(T, n_eqs);
	
	// Use adaptive step control with an absolute error "err"
	gsl_odeiv_control *control = gsl_odeiv_control_y_new(err, 0.0);
	
	// Create evolution function
	gsl_odeiv_evolve *evolve = gsl_odeiv_evolve_alloc(n_eqs);
	
	// Set up the system needed by GSL. No need for jacobian (NULL)
	gsl_odeiv_system sys = {calc_derivs, NULL, n_eqs, params};
	
	// Main loop of the integration. Advance until tmax is reached.
	while (t < tmax)
	{
		n_steps++; // incease the counter of steps at each step
		
		vector<double> row; // create a row of the output
		row.push_back(t); // column 1 - time
		row.push_back(y[0]); // column 2 - theta
		row.push_back(y[1]); // column 3 - omega
		row.push_back(calc_energy(y[0],y[1])); // column 4 - energy
		
		//if (n_steps % 1000 == 0) // uncomment to control the output frequency
		archive.push_back(row); // add entry to the archive
		
		int status = gsl_odeiv_evolve_apply (evolve, control, step,
											 &sys, &t, tmax, &h, y);
		if (status != GSL_SUCCESS) break;
	}
	
	// Tidy up politely
	
	gsl_odeiv_evolve_free(evolve);
	gsl_odeiv_control_free(control);
	gsl_odeiv_step_free(step);
	
	// Print the output archive.
	
	print_archive(archive); 
	
	// Uncomment this to output amplitude and period.
	// This was used to plot the Period vs Amplitude plot. 
	
	// cout<<theta0<<" "<<estimate_period(archive)<<endl;
	
}


// Main

int main(int argc, char* argv[])
{
	// Define the needed variables for the problem.
	// h is the step, err is the desired absolute error
	
	double q, F, theta0, omega0;
	double h, err, tmax;
	
	// Gather up the parameters in an array, so that they can
	// be easily defined as input from the command line.
	// Note: the order of the parameters is not the same as
	// the order in which they are fed into integrate(). It is:
	// q, F, theta0, omega0, h, err, tmax
	// i.e. as defined in the beginning of main().
	
	double params[7] = {0, 0, 0.01, 0, 1e-4, 1e-9, 20};
	
	// You can input only part of the paramenters from the command
	// line and the ones that are left will stay with their default 
	// values.
	
	for (int i=1; i<argc; i++) 
		params[i-1] = cav::string_to_double(argv[i]);
		
	// Assign the values to the corresponding variables.
	
	q = params[0]; F = params[1]; theta0 = params[2]; 
	omega0 = params[3]; h = params[4]; err = params[5]; 
	tmax = params[6]; 
	
	// Call the integrate(...) function, which in this version of
	// the program will output 4 columns: t, theta, omega, energy.
		
	integrate(q, F, h, err, tmax, theta0, omega0);
	
	// Uncomment this and the corresponding expression in integrate(...)
	// if you want to loop through values of theta0 in the range [0.01; 3.13)
	// at steps of 0.01 and output two columns: theta0 and period.
	// This was used to produce the Period vs Amplitude plot.
	
	// for (int i=1; i<313; i++) integrate(q, F, h, err, tmax, ((double)i)/100, omega0);
	
	// Uncomment this and the respective section of integrate(...)
	// if you want to get the value for the period exactly at theta0 = pi/2.
	// The value for pi is defined at the top and is taken from the 
	// cavlib library of constants.
	
	// integrate(q, F, h, err, 20000, pi/2, omega0);
	
	return 0;
}
