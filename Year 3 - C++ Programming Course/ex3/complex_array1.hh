// $Id: complex_array.hh,v 1.1 2009/01/25 11:57:44 jsr Exp jsr $
// Classes for one dimensional arrays of complex numbers

#ifndef cavlib_complex_array1_hh
#define cavlib_complex_array1_hh

#include <vector>
#include <iostream>
#include <cmath>

#include <fftw3.h>

namespace cav
{

// Define a class to represent an array of complex numbers, storing
// the numbers as real, complex, real, complex etc...

class Complex_Array1
{
private:
  
  int                 m_n;    // number of complex points
  std::vector<double> m_data; // storage array length 2 * m_n
  
public:
  
  Complex_Array1( int n ) 
    : m_n( n ), m_data( 2* n, 0.0 )
  {}
  
  int n() const { return m_data.size() / 2; }

  // Pointer to start of data:

  double* data_ptr() { return &m_data[0]; }
  
  // The data array itself:

  std::vector<double>& data() { return m_data; }

  // Access the data as real and imaginary parts:

  double& real( int i ) { return m_data[2*i]; }
  double& imag( int i ) { return m_data[2*i+1]; }

  // Access the data as amp/phase:
  
  double amp( int i ) { return sqrt( pow( real(i), 2 ) + pow( imag(i), 2 ) ); }
  double phi( int i ) { return atan2( imag(i), real(i) ); }
};

// Function to do the FFT:

inline void fft( Complex_Array1& in, Complex_Array1& out, bool forward = true )
{
  fftw_plan plan = fftw_plan_dft_1d ( in.n(), 
				      (fftw_complex*) in.data_ptr(),
				      (fftw_complex*) out.data_ptr(),
				      forward ? FFTW_FORWARD : FFTW_BACKWARD, 
				      FFTW_ESTIMATE );
  fftw_execute ( plan );
  fftw_destroy_plan( plan );

  // Normalise inverse transform:

  if ( !forward )
    for ( unsigned int i = 0; i < out.data().size(); i++ )
      out.data()[i] /= out.n();
}



} // end namespace

#endif
