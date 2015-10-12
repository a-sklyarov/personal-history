// screen.cc - implements the Screen class

/* This is the implementation of the Screen class. A screen has private
 * variables corresponding to the distance to the screen D, a boolean 
 * saying if the pattern has been generated or not, a comlex array
 * containing the pattern in the form of complex numbers, a vector array
 * containing the intensity pattern and a vector array containing the
 * y values corresponding to the pattern. */

#include "screen.hh"

using namespace std;

// constructor
Screen::Screen (double D) : m_D(D), pattern_generated(0), m_pattern(0) {}

// Setters and Getters
void Screen::SetD(double D) { m_D = D; }
double Screen::GetD() { return m_D; }

cav::Complex_Array1 Screen::GetPattern() { return m_pattern; }

vector<double> Screen::GetIntensityPattern() { return m_intensity_pattern; }
		
vector<double> Screen::GetY_Values() { return m_y_values; }

// Methods

// This method takes as an input an aperture which it then
// uses to generate the diffraction pattern from it.

void Screen::GeneratePattern ( Aperture aperture )
{
	if (pattern_generated) return; // generate only once
	
	aperture.discretize();
	
	cav::Complex_Array1 discretized_aperture = aperture.GetDiscretizedRepresentation();
	cav::Complex_Array1 pattern ( aperture.GetN() );
	
	cav::fft( discretized_aperture, pattern ); // call the FFT method 

	m_pattern = pattern; // the raw pattern in the form of a comlex array     
	int N = pattern.n();
	double freq_to_y = ( aperture.GetLambda() * m_D ) ; // use to convert to y values
	
	for (int i = 0; i < N; i++)
	{
		m_intensity_pattern.push_back((pattern.amp(i))*(pattern.amp(i))); // fill in the intensity pattern
		
		// fill in the y values do that the pattern gets symmetrical about 0.
		
		if (i <= (N/2) ) 
		{
			m_y_values.push_back( i * (freq_to_y) / ( N * aperture.GetDelta() ) );
		} 
		else 
		{
			m_y_values.push_back( (((double)i)/N - 1)*(freq_to_y) / aperture.GetDelta() );
		}
	}
	
	pattern_generated = 1;
} 
	


