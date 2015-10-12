// SupplTask1a.cc - the main file for the first part of Supplementary Task 1

// Same as Core Task 1 code with minor modifications to fulfil the
// needs of the task

#include <iostream>
#include <cmath>
#include "cavlib/constants.hh"
#include "aperture.hh"
#include "screen.hh"

using namespace std;

double // constants for Supplementary Task 1a
	lambda = 500e-9, // lambda = 500nm - wavelength
	d = 100e-6, // d = 100 microns - slit width
	D = 5e-3, // D = 5mm - distance to screen
	L = 5e-3; // L = 5mm - maximum extent of aperture
	
Aperture aperture;
Screen screen (0); // the actual D value is set in SetUp_Screen()

complex_number aperture_function_1 ( double x ) // A'(x) in Suppl Task 1a
{
	complex_number result = {0};
	double const pi = (C::pi);
	
	if ( fabs(x) <= d/2 ) 
	{
		result.real_part = cos( 2*pi*x*x/(2*lambda*D) );
		result.imag_part = sin( 2*pi*x*x/(2*lambda*D) );
	}
	else if ( fabs(x) <= L/2 ) result.real_part = 0;
		else if ( fabs(x) > L/2 ) result.real_part = 1;
		
	return result;
}

void SetUp_Aperture( int N_discrete_points )
{
	aperture.SetLambda( lambda );
	aperture.SetN( N_discrete_points );
	aperture.SetMaximumExtent( L );
	aperture.SetDelta();
	
	aperture.SetAnalyticalFunction( aperture_function_1 );
	
	aperture.discretize();
}

void SetUp_Screen()
{
	screen.SetD( D );
	
	screen.GeneratePattern( aperture );
}

void produce_plot_data()
{
		cav::Complex_Array1 
		discretized_aperture_function = aperture.GetDiscretizedRepresentation(),
		screen_pattern = screen.GetPattern();
		
	unsigned int N = screen_pattern.n();
	
	vector<double> 
		x_values = aperture.GetDiscretized_X_Values(), 
		y_values = screen.GetY_Values(), 
		intensity = screen.GetIntensityPattern();
		

	
	for (unsigned int i=0; i < N; i++)
	{	
		// plot only the first ~20 maxima for clearness
		
		if (fabs(y_values.at(i)) < 20*lambda*D / d) {
				
			cout<< x_values.at(i) << " " 						// 1) x
				<< discretized_aperture_function.amp(i) << " "	// 2) Re(A(x))
				<< discretized_aperture_function.phi(i) << " "	// 3) Im(A(x))
				<< y_values.at(i) << " "						// 4) y
				<< intensity.at(i) << " "						// 5) Intensity
				<< screen_pattern.real(i) << " "				// 6) Re(Psi(y))
				<< screen_pattern.imag(i) << " "				// 7) Im(Psi(y))
				<< screen_pattern.amp(i) << " "					// 8) Amp(Psi(y))
				<< screen_pattern.phi(i) << endl;				// 9) Phase(Psi(y))
		}
	}
}

int main()
{
	SetUp_Aperture( two_to_the_20 ); // const is defined in aperture.hh
	SetUp_Screen();
	
	produce_plot_data();

	return 0;
}
