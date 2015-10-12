// simulation.hh

#ifndef SIMUL_HH
#define SIMUL_HH

#include "lattice.hh"

// Function to extract the quilibrium region from a vector array
// of E or M values. This means to get only the second half of the
// values.
vector<long> extract_equilibrium ( vector<long> data);

// Using the function random_uniform() from /ex/prng5.cc
// to produce random numbers in the interval [0,1] with
// uniform distrubution, using <gsl/gsl_rng.h> from the
// GSL Library. The function is used with its default 
// seed in this program.
double random_uniform( unsigned long int seed);

// Function to choose n random values from a vector array.
// This function is used in the implementation of the bootstrap
// method when calculating the error in C = simga^2/kT^2.
vector<long> choose_random ( long n, vector<long> data );
			
// Function to calculate the mean of an array of integral values
double mean_value(vector<long> data);

// Function to calculate the standard deviation squared
double sigma_sqr_value_int( vector<long> data );

// Function to calculate the standard deviation squared
double sigma_sqr_value_double( vector<double> data );
			
// Function to create T=inf initial conditions
vector<long> random_init (long N);

// Function to print a picture of the current spin state 
// which to be used then for visualization
void take_picture ( vector<long> state, long at_time, double at_T );
				
// Function to carry out a whole time step of N^2 flips
void metropolis ( Lattice& latt );

class Simulation
{
	private:
	
		vector<long> M_values; // values of the magnetization
		vector<long> E_values; // values of the energy
		
		vector<double> T_values; // values in the range [T_min, T_max]
		
		vector<double> H_values; // values in the range [H_min, H_max]
		
		vector<double> M_mean_values; // values of the mean magnetization
		vector<double> M_mean_errors; // errors of the mean magnetization
		
		vector<double> E_mean_values; // values of the mean energy
		vector<double> E_mean_errors; // errors of the mean energy
		
		vector<double> C_values_I; // heat capacity values C = dE/dT
		vector<double> C_errors_I; // errors in C = dE/dT
		vector<double> T_values_I; // T values for C = dE/dT
		
		vector<double> C_values_II; // values of C = sigma^2/kT^2
		vector<double> C_errors_II; // errors of C = sigma^2/kT^2
		vector<double> T_values_II; // T values for C = sigma^2/kT^2
						
		// Function to simulate the system lattice for max_time steps
		void simulate ( Lattice& lattice, long max_time );

	public:

		// Methods
		
		// Function to print the data for the current simulation.
		void print_data ( long time_step );
		
		// Function to print the data for the current simulation
		// per site of the lattice.
		void print_data_per_site ( long time_step, long N );
		
		// Function to print the the mean values of M and E over some
		// temperature range.
		void print_T_data ();
		
		// Function to print the the mean values of M and E over some
		// temperature range.
		void print_T_data_per_site ( long N );
		
		// Function to carry out a simulation of a system with 
		// dimension N and temperature T for a number of time steps
		// max_time. The bool init indicates whether the initial 
		// configuration of the spins should correspond to T = 0
		// (i.e. all spins in the same direction) or to T -> inf
		// (i.e. all spins are randomly oriented). The two cases
		// correspond to bool = 0 and 1 respectively.
		void evolve_system ( long N, double T, long max_time, bool init );
		
		// Function to calculate the mean values of the equilibrium
		// energy and magnetization for the current simulation.
		// It assummes that the second half of the interval max_time
		// corresponds to equilibrium.
		void calc_means ();
		
		// Function to calculate the error in the value of C, which
		// has been calculated by C = dE/dT. The absolute error is 
		// simply twice the absolute error of E divided by dT.
		void calc_error_C_I ( double E_mean, double dE, double dT );
		
		// Function to calculate the error in the value of C, which
		// has been calculated by C = sigma^2/kT^2. The bootstrap
		// method is used to achieve this. (see more in report)
		void calc_error_C_II ( double T );
		
		// Function to calculate the heat capacity using the first
		// formula from the instruction handout - C = dE/dT. It uses
		// the stored values of the mean energies and their
		// corresponding temperatures.
		void calc_C_I ();
		
		// Function to calculate the heat capacity using the second
		// formula from the instruction handout - C = sigma^2/kT^2. 
		// It uses the stored values of the mean energies and their
		// corresponding temperatures.
		void calc_C_II ( double T );
		
		// Function to carry out a series of simulations of a system
		// for different temperatures in a range [T_min, T_max] with 
		// a step of T_step. For every individual simulation the mean
		// energy and magnetization of the equilibrium state are 
		// measured. 
		void evolve_in_T ( long N, double T_min, double T_max, double T_step, double H );		
		
};
		

#endif
