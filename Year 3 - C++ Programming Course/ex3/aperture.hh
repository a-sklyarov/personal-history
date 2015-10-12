// aperture.hh - defines the Aperture class

/* This is the definition of the Aperture class. An Aperture has
 * private variables corresponding to the wavelnegth lambda, the
 * number of discrete fourier sample regions N, its maximum extent 
 * L, the size of a sample region delta, a boolean saying if the 
 * aperture has been discretized or not, an analytical function 
 * representing its complex transmittance, a complex array containing
 * its discretized representation (if it has been discretized) and a 
 * vector array containing the x values corresponding to the discretized
 * complex transmittance. 
 */

#ifndef APERTURE_HH
#define APERTURE_HH

#include <iostream>
#include <vector>
#include "complex_array1.hh"

using namespace std;

const int two_to_the_30 = 1073741824;
const int two_to_the_20 = 1048576;

// Define the type complex_number which is expected to be returned
// by the function which represents the analytical complex_transmittance
// of the aperture. The analytical function will be discretized later.

typedef struct { double real_part, imag_part; } complex_number;
typedef complex_number (*complex_transmittance)(double x);

class Aperture
{
	private:
		double m_lambda;
		int m_N;
		double m_maximum_extent; // L
		double m_delta;
		bool m_discretized;
		
		complex_transmittance m_analytical_function; // A(x)
				
		cav::Complex_Array1 m_discretized_representation;	
		vector<double> m_discretized_x_values; // size N	
		
	public:
		// constructors
		Aperture();
		Aperture(double lambda_in, int N, double maximum_extent_in);
		
		// Setters and Getters
		void SetLambda ( double lambda_in );
		double GetLambda ();
		
		void SetN ( double N_in);
		int GetN ();
		
		void SetMaximumExtent ( double max_ext_in );
		double GetMaximumExtent ();
		
		void SetDelta ();
		double GetDelta ();
		
		void SetAnalyticalFunction ( complex_transmittance a );
		
		cav::Complex_Array1 GetDiscretizedRepresentation ();
		
		vector<double> GetDiscretized_X_Values ();
		
		// Methods
		
		// This method will discretize the analytical complex
		// transmittance.
		
		void discretize ();
		
};

#endif
		
		
