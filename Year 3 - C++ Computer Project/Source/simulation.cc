// simulation.cc

#include <iostream>
#include <cmath>
#include "gsl/gsl_rng.h"
#include "simulation.hh"

using namespace std;

// Function to extract the quilibrium region from a vector array
// of E or M values. This means to get only the second half of the
// values.
vector<long> extract_equilibrium ( vector<long> data)
{
	long size = data.size();
	vector<long> result ( data );
	
	result.erase ( result.begin(), result.begin() + 2*size/3 );
	
	return result;
}

// Using the function random_uniform() from /ex/prng5.cc
// to produce random numbers in the interval [0,1] with
// uniform distrubution, using <gsl/gsl_rng.h> from the
// GSL Library. The function is used with its default 
// seed in this program.
double random_uniform( unsigned long int seed = 0 )
{
  static gsl_rng* rng = 0;
  
  if ( rng == 0 ) rng = gsl_rng_alloc( gsl_rng_default ); 
  if ( seed != 0 ) gsl_rng_set( rng, seed );

  return gsl_rng_uniform( rng );
}

// Function to choose n random values from a vector array.
// This function is used in the implementation of the bootstrap
// method when calculating the error in C = simga^2/kT^2.
vector<long> choose_random ( long n, vector<long> data )
{
	long size = data.size();
	long el = 0;
	vector<long> result;
	
	for ( long i = 0; i < n; i++ )
	{
		el = (long) ( size * random_uniform() );
		result.push_back( data.at( el ) );
	}
	
	return result;
}

// This function calculates the mean of an array of integral values
double mean_value(vector<long> data)
{
	long sum = 0;
	long nt = data.size();
	
	for (long i = 0; i < nt; i++) sum += data.at(i);
	
	return ((double)sum)/nt;
}

// Function to calculate the standard deviation squared
double sigma_sqr_value_int( vector<long> data )
{
	/*
	long sum = 0, ssum = 0;
	long nt = data.size();
	
	for (long i = 0; i< nt; i++) 
	{
		long el = data.at(i);
		
		sum += el;
		ssum += el*el;
	}
	
	double sigma = ((double)ssum)/nt - (((double)sum)/nt)*(((double)sum)/nt);
	
	return fabs( sigma );
	* */
	
	double mean = mean_value ( data );
	long size = data.size();
	
	double sigma_sqr = 0;
	
	for ( int i = 0; i < size; i++ )
	{
		long el = data.at(i);
		sigma_sqr += ( el - mean ) * ( el - mean );
	}
	
	sigma_sqr /= size;
	
	return sigma_sqr;
	
}

// Function to calculate the standard deviation
double sigma_sqr_value_double( vector<double> data )
{
	double sum = 0, ssum = 0;
	long nt = data.size();
	
	for (long i = 0; i< nt; i++) 
	{
		double el = (double)data.at(i);
		
		sum += el;
		ssum += el*el;
	}
	
	double sigma = ssum/nt - (sum/nt)*(sum/nt);
	
	return sigma;
}

// Function to create T=inf initial conditions
vector<long> random_init (long N)
{
	vector<long> rand;
	double r = 0.0;
	
	for (long i = 0; i < N*N; i++)
	{
		r = random_uniform();
		
		if ( r < 0.5 ) rand.push_back( -1 );
		else rand.push_back( 1 );
	}
	
	return rand;
}

// Function to print a picture of the current spin state 
// which to be used then for visualization
void take_picture ( vector<long> state, long at_time, double at_T )
{
	long N = (long) ( sqrt( (double)state.size() ) );
	
	cout << "System size: " << N << " x " << N << ", at time: " 
		 << at_time << ", at temperature: " << at_T << endl;
	
	for ( long i = 0; i < N; i++ )
	{
		for ( long j = 0; j < N; j++ )
		{
			cout<< state.at( i*N + j ) <<" ";
		}
		
		cout << endl;
	}
	
	cout << endl;
}		
		

// a whole time step of N^2 flips
void metropolis ( Lattice& latt )
{
	
	long ran_row = 0; // to carry a random row number
	long ran_col = 0; // to carry a ranodm col number
	double p = 0.0; // to carry a random number to decide flipping
	long N = latt.GetN(); // dimension of lattice
	
	for (long i = 0; i < (N*N); i++)
	{
		ran_row = (long) ( random_uniform() * (double)N );
		ran_col = (long) ( random_uniform() * (double)N );
		p = random_uniform();
		
		latt.flip(ran_row, ran_col, p);
	}
	
}

// Function to simulate the system lattice for max_time steps
void Simulation::simulate ( Lattice& lattice, long max_time )
{
	long M = 0;
	long E = 0;
	
	M_values.clear();
	E_values.clear();
	
	for (long t = 0; t < max_time; t++)
	{
		M = lattice.GetM();
		E = lattice.GetE();
		
		/* This snippet was used to print lattice configurations
		 * for visualization. Not needed for ant other activities.
		
		if ( ( t == 0) || ( t == max_time / 1000 ) || 
			( t == max_time / 100 ) || ( t == max_time / 10 )
			|| ( t == (max_time - 1) ) )
		{
			double T = lattice.GetT();
			take_picture ( lattice.GetState(), t, T );
		} 
		* Note: This must be used only if no other print options
		* 		 are specified in the main program.
		* */
		
		M_values.push_back( M );
		E_values.push_back( E );
		
		metropolis( lattice );
	}
}

// Function to print the data for the current simulation
void Simulation::print_data ( long time_step )
{
	long size = E_values.size();
	
	for ( long t = 0; t < size; t += time_step )
	{
		cout<<t<<" "<<E_values.at(t)<<" "<<M_values.at(t)<<endl;
	}
}

// Function to print the data for the current simulation
// per site of the lattice.
void Simulation::print_data_per_site ( long time_step, long N )
{
	long size = E_values.size();
	
	for ( long t = 0; t < size; t += time_step )
	{
		cout<<t<<" "<<((double)E_values.at(t))/(N*N)<<
				 " "<<((double)M_values.at(t))/(N*N)<<endl;
	}
}

// Function to print the the mean values of M and E over some
// temperature range.
void Simulation::print_T_data ()
{
	long size = T_values.size();
	
	// T_values_I, C_values_I, C_errors_I have one element less, so
	// we can just add the last one once more to make the
	// printing easier
	T_values_I.push_back( T_values_I.at(size-2) );
	C_values_I.push_back( C_values_I.at(size-2) );
	C_errors_I.push_back( C_errors_I.at(size-2) );
	
	for ( long i = 0; i < size; i++ )
	{
		
			cout<< T_values.at(i) <<" " 				// 1)
				<< E_mean_values.at(i) <<" "			// 2)
				<< E_mean_errors.at(i) <<" "			// 3)
				<< fabs(M_mean_values.at(i)) <<" "		// 4)
				<< M_mean_errors.at(i) <<" "			// 5)
				<< T_values_I.at(i) <<" "				// 6)
				<< C_values_I.at(i) <<" "				// 7)
				<< C_errors_I.at(i) <<" "				// 8)
				<< T_values_II.at(i) <<" "				// 9)
				<< C_values_II.at(i) <<" "				// 10)
				<< C_errors_II.at(i) <<endl;			// 11)
		
	}
}

// Function to print the the mean values of M and E over some
// temperature range.
void Simulation::print_T_data_per_site ( long N )
{
	long size = T_values.size();
	
	// T_values_I, C_values_I, C_errors_I have one element less, so
	// we can just add the last one once more to make the
	// printing easier
	T_values_I.push_back( T_values_I.at(size-2) );
	C_values_I.push_back( C_values_I.at(size-2) );
	C_errors_I.push_back( C_errors_I.at(size-2) );
	
	for ( long i = 0; i < size; i++ )
	{
		
			cout<< T_values.at(i) <<" " 					// 1)
				<< E_mean_values.at(i) / (N*N)<<" "			// 2)
				<< E_mean_errors.at(i) / (N*N)<<" "			// 3)
				<< fabs(M_mean_values.at(i)) / (N*N)<<" "	// 4)
				<< M_mean_errors.at(i) / (N*N)<<" "			// 5)
				<< T_values_I.at(i) <<" "					// 6)
				<< C_values_I.at(i) / (N*N)<<" "			// 7)
				<< C_errors_I.at(i) / (N*N)<<" "			// 8)
				<< T_values_II.at(i) <<" "					// 9)
				<< C_values_II.at(i) / (N*N) <<" "			// 10)
				<< C_errors_II.at(i) / (N*N)<<endl;			// 11)
		
	}
}

// Function to calculate the mean values of the equilibrium
// energy and magnetization for the current simulation.
// It assummes that the second half of the interval max_time
// corresponds to equilibrium.
void Simulation::calc_means ()
{
	// equilibrium regions of E and M
	vector<long> eq_M_values ( extract_equilibrium( M_values ) );
	vector<long> eq_E_values ( extract_equilibrium( E_values ) );
	
	// calculate the means and add them to the corresponding vectors
	double mean_M = mean_value ( eq_M_values );
	double mean_E = mean_value ( eq_E_values );
	
	// calculate the errors in the means
	double error_M = sqrt ( sigma_sqr_value_int( eq_M_values ) );
	double error_E = sqrt ( sigma_sqr_value_int( eq_E_values ) );
	
	// record the mean values
	M_mean_values.push_back ( mean_M );
	E_mean_values.push_back ( mean_E );
	
	// record the errors in the means
	M_mean_errors.push_back ( error_M );
	E_mean_errors.push_back ( error_E );
}

// Function to calculate the error in the value of C, which
// has been calculated by C = dE/dT. The relative error is 
// simply the relative error of E divided by dT.
void Simulation::calc_error_C_I ( double E_mean, double dE, double dT )
{
	// equilibrium region of E values
	vector<long> eq_E_values ( extract_equilibrium( E_values ) );
	
	// absolute error in the equilibrium value of E
	double sigma = sqrt( sigma_sqr_value_int( eq_E_values ) );
	
	// relative error in E
	double R_E = sigma / E_mean;
	
	// relative error in C
	double R_C = ( R_E / dT );
	
	// absolute error in C
	double error = R_C * ( dE / dT );
	
	C_errors_I.push_back( error );
}

// Function to calculate the error in the value of C, which
// has been calculated by C = sigma^2/kT^2. The bootstrap
// method is used to achieve this. (see more in report)
void Simulation::calc_error_C_II ( double T )
{
	// equilibrium region of E values
	vector<long> eq_E_values ( extract_equilibrium( E_values ) );
	
	// number of resamplings to be done
	const long N_resample = 20;
	
	// 1/fraction of the whole array to be used at each resampling
	const long fraq = 4; // i.e. 1/4 of the size of the array
	
	// size of the array to be resampled
	long size = eq_E_values.size();
	
	// number of entries corresponding to the fraction specified above
	long n = size / fraq;
	
	// vector array to contain the calculated C values after each 
	// resampling cycle
	vector<double> c_values;
	
	// do the resamplings
	for ( long i = 0; i < N_resample; i++ )
	{
		vector<long> sample ( choose_random( n, eq_E_values ) );
		double sigma_sqr = sigma_sqr_value_int(sample);
		
		c_values.push_back( sigma_sqr / (T*T) );
	}
	
	double sigma = sqrt ( sigma_sqr_value_double( c_values ) );
	
	C_errors_II.push_back( sigma );
}
	

// Function to calculate the heat capacity using the first
// formula from the instruction handout - C = dE/dT. It uses
// the stored values of the mean energies and their
// corresponding temperatures. The method consists in taking
// the difference of two neighbouring energies and dividing
// it by the difference of the two corresponding neighbouring
// temperatures. The resulting value is the heat capacity at the
// average of the two temperatures.
void Simulation::calc_C_I ()
{
	long size = T_values.size();
	
	double dE = 0;
	double dT = 0;
	
	double E_mean = 0;
	
	for ( long i = 0; i < (size - 1); i++ )
	{
		dE = fabs(E_mean_values.at(i+1) - E_mean_values.at(i));
		dT = T_values.at(i+1) - T_values.at(i);
		
		C_values_I.push_back( dE/dT );
		T_values_I.push_back( (T_values.at(i+1) + T_values.at(i)) / 2 );
		
		E_mean = ( E_mean_values.at(i+1) + E_mean_values.at(i) ) / 2;
		
		calc_error_C_I( E_mean, dE, dT );
	}
}
		
// Function to calculate the heat capacity using the second
// formula from the instruction handout - C = sigma^2/kT^2. 
// It uses the stored values of the mean energies and their
// corresponding temperatures.
void Simulation::calc_C_II ( double T )
{
	vector<long> eq_E_values ( extract_equilibrium ( E_values ) );
	
	double sigma_sqr = sigma_sqr_value_int( eq_E_values );
	
	C_values_II.push_back ( sigma_sqr / ( T*T ) );
	T_values_II.push_back ( T );
	
	calc_error_C_II ( T );
}

// Function to carry out a simulation of a system with 
// dimension N and temperature T for a number of time steps
// max_time. The bool init indicates whether the initial 
// configuration of the spins should correspond to T = 0
// (i.e. all spins in the same direction) or to T -> inf
// (i.e. all spins are randomly oriented). The two cases
// correspond to bool = 0 and 1 respectively.
void Simulation::evolve_system ( long N, double T, long max_time, bool init )
{
	Lattice lattice (N, T, 0);
	
	if ( init ) lattice = Lattice( N, T, random_init(N), 0 );
	
	simulate ( lattice, max_time );
}

// Function to carry out a series of simulations of a system
// for different temperatures in a range [T_min, T_max] with 
// a step of T_step. For every individual simulation the mean
// energy and magnetization of the equilibrium state are 
// measured. Also, the heat capacity is estimated in the two
// ways suggested in the instruction handout. 
//NOTE: At every temperature step, the simulated system is
// initialized with the equilibrium lattice state of the system
// from the previous step. This helps for a much faster achievement
// of equilibrium compared to the case where the initial state
// is T=0 or T->inf at every step.
void Simulation::evolve_in_T ( long N, double T_min, double T_max, double T_step, double H )
{
	// The time to reach equilibrium should not be more than N*N
	// because by that time all spins will have flipped at least once
	// on average. Experimentation shows that this time is usually
	// less than N*N/2 and if you add the fact that at each step
	// the starting configuration is the equilibrium state of the
	// previous step, max_time = N*N should be more than enough to 
	// have equilibrium at least during the second half of it.
	
	long max_time = 2.5 * N * N; // equilibrium during the second half
	
	long steps = (long)( (T_max - T_min) / T_step ); // number of steps
	
	Lattice lattice ( N, T_min, H ); // initialize ( try random_init(N) )
	
	T_values.clear(); // to be safe
	
	for ( long i = 0; i < steps; i++ )
	{		
		double T = T_min + i*T_step; // step forward
		
		// use the lattice configuration reached by the system
		// from the previous iteration to achieve equilibrium faster
		
		lattice = Lattice( N, T, lattice.GetState(), H );
		
		simulate ( lattice, max_time );
		
		// take note of the current T
		
		T_values.push_back( T );
		
		// Calculate the heat capacity for the current T
		
		calc_C_II( T );
		
		// Measure mean energy and magnetization
		
		calc_means();
	}
	
	// Calculate heat capacity
	
	calc_C_I();
}









