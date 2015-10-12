// main.cc

#include <ctime>
#include "simulation.hh"
#include "cavlib/string_util.hh"

//-------------------------------------------------------------------//

void performance_test1 ( long N_min, long N_max, long N_step, double T, long max_time, bool init )
{
	using namespace std;
	
	// to mark the begin and end times
	clock_t begin, end;
	// number of steps 	
	long steps = ( N_max - N_min ) / N_step;
	
	for ( int i = 0; i < steps; i++ )
	{
		long N = N_min + i*N_step; // step forward
		
		cout << N << " ";
		
		begin = clock(); // take initial time
		
		Simulation sim;		
		sim.evolve_system( N, T, max_time, init );

		end = clock(); // take final time		
		
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		cout << elapsed_secs << endl;
	}
}

void performance_test2 ( long N_min, long N_max, long N_step, 
						 double T_min, double T_max, double T_step)
{
	using namespace std;
	
	// to mark the begin and end times
	clock_t begin, end;
	// number of steps 	
	long steps = ( N_max - N_min ) / N_step;
	
	for ( int i = 0; i < steps; i++ )
	{
		long N = N_min + i*N_step; // step forward
		
		cout << N << " ";
		
		begin = clock(); // take initial time
		
		Simulation sim;		
		sim.evolve_in_T ( N, T_min, T_max, T_step, 0 );

		end = clock(); // take final time		
		
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		
		cout << elapsed_secs << endl;
	}
}

//------------------------------------------------------------------//

int main( int argc, char* argv[])
{	
	if ( cav::string_to_int(argv[1]) == 0 )
	{
	// This is help screen showing the available options.
	
		cout<<"This program offers five possible regimes of work: "<<endl
			<<" 1) If you want to perform a time evolution simulation of a specified system (Task 1), you should start your input with 1, followed by four values, corresponding to N (dimension of the lattice), T (temperature), max_time (number of time steps for which the system will be simulated), init (0 or 1 determining the initial conditions to be T=0 or T=inf respectively)."<<endl
			<<" 2) If you want to perform a temperature evolution of a specified system (Task 2, 3, 4), you should start your input with 2, followed by four values, corresponding to N (dimension of the lattice), T_min (starting temperature), T_max (final temperature), T_step (temperature step with which the evolution will proceed)."<<endl
			<<" 3) If you want to perform a temperature evolution of a specified system in non-zero magnetic field (Task 5), you should start your input with 3, followed by five values, corresponding to N (dimension of the lattice), T_min (starting temperature), T_max (final temperature), T_step (temperature step with which the evolution will proceed), H (the value of the magnetic field)."<<endl
			<<" 4) If you want to perform the first performance test, you should start your input with 4, followed by six values, corresponding to N_min (starting dimension of the lattice), N_max (final dimenstion of the lattice), N_step (step in the dimension with which the test will proceed), T (temperature), init (0 or 1 determining the initial conditions to be T=0 or T=inf respectively)."<<endl
			<<" 5) If you want to perform the first performance test, you should start your input with 5, followed by six values, corresponding to N_min (starting dimension of the lattice), N_max (final dimenstion of the lattice), N_step (step in the dimension with which the test will proceed), T_min (starting temperature), T_max (final temperature), T_step (temperature step with which the evolution will proceed)."<<endl<<endl;
					
	} else if ( cav::string_to_int(argv[1]) == 1 )
	{
	// To perform time evolution of a specified system. Related to
	// Task 1. from the instructions handout.
		long N = cav::string_to_int( argv[2] );
		double T = cav::string_to_double( argv[3] );
		long max_time = cav::string_to_int( argv[4] );
		bool init = (bool) cav::string_to_int( argv[5] );
	
		Simulation sim;
	
		sim.evolve_system( N, T, max_time, init );
		sim.print_data( 1 + max_time / 500 );
		
		// uncomment this if you want to print the same things
		// but with their values per spin
		//sim.print_data_per_site( 1 + max_time / 500, N );
		
	} else if ( cav::string_to_int(argv[1]) == 2 )
	{
	// To perform temperature evolution of a specified system.
	// Related to Tasks 2., 3. and 4. from the instructions 
	// handout. Produce temperature dependence of various quantities.
		long N = cav::string_to_int( argv[2] );
		double T_min = cav::string_to_double( argv[3] );
		double T_max = cav::string_to_double( argv[4] );
		double T_step = cav::string_to_double( argv[5] );
	
		Simulation sim;
		
		sim.evolve_in_T ( N, T_min, T_max, T_step, 0 );
		sim.print_T_data ();
		
		// uncomment this if you want to print the same things
		// but with their values per spin		
		//sim.print_T_data_per_site( N );
		
	} else if ( cav::string_to_int(argv[1]) == 3 )
	{
	// Same as the previous option, just this time you can also
	// specify a magnetic field H in the last parameter. Related to 
	// Task 5. from the instructions handout.
		long N = cav::string_to_int( argv[2] );
		double T_min = cav::string_to_double( argv[3] );
		double T_max = cav::string_to_double( argv[4] );
		double T_step = cav::string_to_double( argv[5] );
		double H = cav::string_to_double( argv[6] );
	
		Simulation sim;
		
		sim.evolve_in_T ( N, T_min, T_max, T_step, H );
		sim.print_T_data ();
		
		// uncomment this if you want to print the same things
		// but with their values per spin	
		//sim.print_T_data_per_site( N );
		
	} else if ( cav::string_to_int(argv[1]) == 4 )
	{
	// Performance test 1
		long N_min = cav::string_to_int( argv[2] );
		long N_max = cav::string_to_int( argv[3] );
		long N_step = cav::string_to_int( argv[4] );
		double T = cav::string_to_double( argv[5] );
		long max_time = cav::string_to_int( argv[6] );
		bool init = (bool) cav::string_to_int( argv[7] );
		
		performance_test1 ( N_min, N_max, N_step, T, max_time, init );
		
	} else if ( cav::string_to_int(argv[1]) == 5 )
	{
	// Performance test 2
		long N_min = cav::string_to_int( argv[2] );
		long N_max = cav::string_to_int( argv[3] );
		long N_step = cav::string_to_int( argv[4] );
		double T_min = cav::string_to_double( argv[5] );
		double T_max = cav::string_to_double( argv[6] );
		double T_step = cav::string_to_double( argv[7] );
		
		performance_test2 ( N_min, N_max, N_step, T_min, T_max, T_step );
		
	} else cout<<" Error: invalid input "<<endl<<"Input 0 for help,"<<endl;
	
	return 0;
}
