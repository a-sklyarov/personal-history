#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>

using namespace std;

#define ranf() \
  ((double)random()/(1.0+(double)RAND_MAX)) // Uniform from interval [0,1) */

struct particle {
  double x    ; // x coordinate
  double im   ; // inverse mass
  double v    ; // velocity
  double p    ; // momentum
  double T    ; // energy
} ; // Note the definition of a structure ends with a semicolon

void v2p( particle &a )
  // propagate a velocity into a momentum
  // p = m * v
{
  if( a.im > 0.0 ) 
    a.p =  a.v / a.im ;
  else // infinite mass particle can't be handled correctly
    a.p = 0.0 ;

  // update the energy too
  a.T = 0.5*a.v * a.p ; 
}

void collide ( particle & a , particle & b  )
  // two particles collide elastically
{
  // find the relative velocity 
  double velocity_along_line =  a.v - b.v;
  // find a's inverse-mass fraction
  double af = a.im / ( a.im + b.im ) ; 
  double bf = b.im / ( a.im + b.im ) ; 
  // reverse the c.o.m. velocity of each along the line of collision
  a.v -= 2.0 * af * velocity_along_line ;
  b.v += 2.0 * bf * velocity_along_line ;
  // The collision has now happened.
  // Finally update the other state information
  v2p(a); 
  v2p(b);
}

// Test the collide function using random initial conditions
int main(int argc, char* argv[])
{
  particle   a , b ;
  int Tests = 10 ; 
  // set defaults
  // read in any command-line arguments
  if(argc>1)   {
    sscanf( argv[1], "%d", &Tests ) ; // put the first command-line argument 
  }

  for ( int t = 1 ; t <= Tests ; t ++ ) {
    a.im    = ranf() + 0.1 ; // random inverse-mass
    a.v     = ranf()-0.5 ;   // random velocity
    b.im    = ranf() + 0.1 ; // random inverse-mass
    b.v     = ranf()-0.5 ;   // random velocity
    v2p(a);
    v2p(b);
    double MomentumBefore = a.p + b.p ;
    double EnergyBefore = a.T + b.T ;
    cout << "a.m:" << 1.0/a.im << "," << "b.m:" << 1.0/b.im << endl ;
    cout << a.v << "," << b.v << " -> " ;
    collide(a,b) ; 
    cout << a.v << "," << b.v << endl;
    double MomentumAfter = a.p + b.p ;
    double EnergyAfter = a.T + b.T ;
    cout << "(E,P) before: "
	 << EnergyBefore << ", "
	 << MomentumBefore << ";" << endl ; 
    cout << "(E,P) after:  "
	 << EnergyAfter << ", "
	 << MomentumAfter << "." << endl ;
    cout << "Percentage change in energy: "
	 << 100.0*(EnergyAfter-EnergyBefore)/
      (0.5*(EnergyAfter+EnergyBefore)) << endl; 
    cout << "Percentage change in momentum: "
	 << 100.0*(MomentumAfter-MomentumBefore)/
      (0.5*(MomentumAfter+MomentumBefore)) << endl;
    cout << endl; 
  }
  return 0;
}
