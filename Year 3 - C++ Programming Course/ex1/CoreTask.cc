#include <iostream>
#include <cmath>
#include <vector>
#include "gsl/gsl_rng.h"
#include "cavlib/constants.hh"

using namespace std;

int const n = 8; //number of dimensions
double const s = (C::pi)/n; // Defining the upper boundary s of the integrals.
double const V = pow(s, n); // define the 8 dimensional volume

double x[n] = {0}; // arguments x0 to x7 for the integrand
double mean_f = 0, mean_squared_f = 0; // mean and mean squared of f

// Using the function random_uniform() from /ex/prng5.cc,
// as suggested in the manual, to produce random numbers
// in the interval [0,1] with uniform distrubution, using 
// <gsl/gsl_rng.h> from the GSL Library. The function is 
// used with its default seed in this program.
double random_uniform( unsigned long int seed = 0 )
{
  static gsl_rng* rng = 0;
  
  if ( rng == 0 ) rng = gsl_rng_alloc( gsl_rng_default ); 
  if ( seed != 0 ) gsl_rng_set( rng, seed );

  return gsl_rng_uniform( rng );
}


// This function fills the array x[] with n=8 random
// values in the interval [0,s]. The array x[] represents
// the arguments for the integrand function.
void generate_arguments()
{
	for (int i=0; i<n; i++) x[i] = s*random_uniform();
}


// This function represents the integrand function.
double f(double x[])
{
	double sum = 0; // the sum of the 8 arguments
	
	for (int i=0; i<n; i++) sum += x[i]; // calculate the sum
	
	return sin(sum);
}


// This function calculated the mean and mean squared values
// of the integrand function f, i.e. <f> and <f^2>. It takes
// as an argument the number of Monte Carlo sample points, N.
void calculate_means(int N) 
{
	// initialize the sum of the function values
	// and the sum of their squares
	double sum_f = 0, sum_f_sq = 0; 
	
	//variable to hold the value of f(x), so that
	//it is calculated only once per cycle
	double f_x;
	
	for (int i=0; i < N; i++)
	{
		generate_arguments();
		f_x = f(x);
		sum_f += f_x;
		sum_f_sq += f_x * f_x;
	}
	
	mean_f = sum_f / N;
	mean_squared_f = sum_f_sq / N;
}


// This function returns an estimate of the integral's value
// using N sample points
double estimate_integral(int N) 
{	
	calculate_means(N);
	return 1000000*V*mean_f; // estimate of the integral
}


// This function returns an estimate of the error in the
// integral using N sample points
double estimate_error(int N)
{
	return 1000000*V*sqrt((mean_squared_f - mean_f*mean_f)/N);
}


// This function calculates the mean of an array of integral values
double best_value(vector<double> integrals)
{
	double sum = 0;
	int nt = integrals.size();
	
	for (int i = 0; i < nt; i++) sum += integrals.at(i);
	
	return sum/nt;
}


// This function calculates the RMS of an array of error estimates
double rms_error(vector<double> errors)
{
	double sum = 0;
	int nt = errors.size();
	
	for (int i = 0; i < nt; i++) sum += errors.at(i) * errors.at(i);
	
	return sqrt(sum/nt);
}


int main()
{	
	vector<double> integrals; // array of integral values
	vector<double> errors; // array of error estimates
	int nt = 25; // number of integral estimates to be made
	int N; // number of sample points
	const int N_max = 1e8; // could be a long wait with this value :)
	
	cout.precision(4); // set an appropriate precision after experiments
	
	// Loop through exponentially increasing values of N and
	// find a best esimate and error for each.
	for (N = 10; N < N_max; N = N*2)
	{
		// Collect nt estimates of the integral
		for (int i = 0; i < nt; i++) 
		{
			integrals.push_back(estimate_integral(N));
			errors.push_back(estimate_error(N));
		}
		
		// Print the best estimate and error for the collection
		// of integral values for the current N
		cout<<N<<" "<<fixed<<best_value(integrals)<<" "<<rms_error(errors)<<endl;
		
		// Empty the arrays so that they can be used again in 
		// the next cycle of the loop
		integrals.clear();
		errors.clear();
	} 
	
	return 0;
}
