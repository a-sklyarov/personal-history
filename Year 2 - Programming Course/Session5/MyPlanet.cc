// Here begins MyPlanet.cc

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;


#define ranf() \
  ((double)random()/(1.0+(double)RAND_MAX)) // Uniform from interval [0,1) */

#define D  2  // number of dimensions
struct particle {
  double x[D] ; // (x,y) coordinates
  double p[D] ; // momentum
  double F[D] ; // force
  double im   ; // inverse mass
  double GMm  ; // gravitational parameter of this particle
  double v[D] ; // velocity
  double T    ; // kinetic energy
  double V    ; // potential energy
  double r    ; // absolute distance from origin
} ; // Note the definition of a structure ends with a semicolon

struct control {
  int verbose ; // program verbosity
  int printing ; // period with which to print
};

void Force( particle &a )
  // sets the force vector and the potential energy
{
  double R = 0.0 ; 
  for ( int i = 0 ; i < D ; i++ ) {
    R += a.x[i] * a.x[i] ;
  }
  double r = sqrt(R) ;
  a.V = - a.GMm / r ;
  a.r = r ; 
  for ( int i = 0 ; i < D ; i++ ) {
    a.F[i] = - a.GMm * a.x[i] / (r*R) ;  // inverse sq
  }
}

void PositionStep ( particle &a , double dt )
{ // the "&" indicates 'this is a pass by reference'
  for ( int i = 0 ; i < D ; i++ ) { 
    a.x[i] += dt * a.p[i] * a.im ; 
  } 
}

void MomentumStep ( particle &a , double dt )
{ // the "&" indicates 'this is a pass by reference. 
  for ( int i = 0 ; i < D ; i++ ) { 
    a.p[i] += dt * a.F[i] ;
  }
}

void v2p( particle &a )
  // propagate a change in velocity into the momentum vector
{
  for ( int i = 0 ; i < D ; i++ ) {
    a.p[i] =  a.v[i] / a.im ;
  }
}

void p2v ( particle &a ) {
  for ( int i = 0 ; i < D ; i++ ) { 
    a.v[i]  = a.p[i] * a.im ;
  }
}

void pv2T ( particle &a ) {
  a.T=0.0;
  for ( int i = 0 ; i < D ; i++ ) { 
    a.T    += 0.5*a.v[i] * a.p[i] ; 
  }
}

void showState ( particle &a )
{ 
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\t"<<a.x[i];
  }
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\t"<<a.v[i];
  }
  cout << "\t" << a.T << "\t" << a.V << "\t" << a.T+a.V << endl;
  cout << "#" ;
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\tx["<<i<<"]\t";
  }
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\tv["<<i<<"]\t";
  }
  cout << "\tT\t\tV\t\tT+V" ; 
  cout << endl;
}

void LeapfrogDynamics( particle &a, double dt, double &t, int N, control &c)
{
  for ( int i=0 ; i < N ; i ++ ) {
    if( c.printing && !(i%c.printing) ) {
      // this iteration we will print out the state
      p2v(a);   pv2T(a) ;     Force( a ) ;
      cout << t ; showState( a ) ;
    }
    // each iteration we move the position a half-step
    PositionStep( a , dt*0.5 ) ;
    Force( a ) ; // (computes the force at this position)
    // then we move the momentum full-step
    MomentumStep( a , dt )     ;
    // then we move the position a half-step
    PositionStep( a , dt*0.5 ) ;
    t += dt ;
  }
}

void EulerDynamics( particle &a, double dt, double &t, int N, control &c)
{
  for ( int i=0 ; i < N ; i ++ ) {
    if( c.printing && !(i%c.printing) ) {
      // this iteration we will print out the state
      p2v(a);   pv2T(a) ;     Force( a ) ;
      cout << t ; showState( a ) ;
    }
    Force( a ) ; // (computes the force at this position)
    PositionStep( a , dt ) ;
    MomentumStep( a , dt ) ;
    t += dt ;
  }
}


void makeCircular(particle &a) // changes the velocity of a particle s that its orbit becomes circular
{
	double R = sqrt((a.x[0])*(a.x[0])+(a.x[1])*(a.x[1])); // distance from the sun
	
	a.v[0] = sqrt((a.GMm*a.im)/(R*(1+(a.x[0]/a.x[1])*(a.x[0]/a.x[1])))); // follows from newton's law
	a.v[1] = -sqrt((a.GMm*a.im)/(R*(1+(a.x[1]/a.x[0])*(a.x[1]/a.x[0])))); // follows from newton's law
}

// Writes the state into a file
void fileState ( particle &a, char fileName[], double t )
{ 
  ofstream fout;
  fout.open(fileName, ios::app);

  fout << t ;

  for ( int i = 0 ; i < D ; i++ ) { 
    fout << "\t"<<a.x[i];
  }
  for ( int i = 0 ; i < D ; i++ ) { 
    fout << "\t"<<a.v[i];
  }
  fout << "\t" << a.T << "\t" << a.V << "\t" << a.T+a.V << endl;
  fout << "#" ;
  for ( int i = 0 ; i < D ; i++ ) { 
    fout << "\tx["<<i<<"]\t";
  }
  for ( int i = 0 ; i < D ; i++ ) { 
    fout << "\tv["<<i<<"]\t";
  }
  fout << "\tT\t\tV\t\tT+V" ; 
  fout << endl;

  fout.close();
}

// This is a leapfrog method which writes in a file
void LeapfrogToFile ( particle &a, double dt, double t, int N, control &c, char fileName[] )
{
  ofstream fout;
  fout.open(fileName);
  fout.close();
  // first we delete the previous content of the file

  for ( int i=0 ; i < N ; i ++ ) {
    if( c.printing && !(i%c.printing) ) {
      // this iteration we will print out the state
      p2v(a);   pv2T(a) ;     Force( a ) ;
      fileState( a, fileName, t ) ;
    }
    // each iteration we move the position a half-step
    PositionStep( a , dt*0.5 ) ;
    Force( a ) ; // (computes the force at this position)
    // then we move the momentum full-step
    MomentumStep( a , dt )     ;
    // then we move the position a half-step
    PositionStep( a , dt*0.5 ) ;
    t += dt ;
  }
}

// Gives random direction of the velocity 'v' of a particle
void RandomVelocityDirection(double v, particle &a, int seed) //v is the absolute value of the velocity we need
{
	srandom(time(0)+seed); //initialise random generator
	a.v[0] = ranf()*v; //uniform in the range [0, v)
	a.v[1] = sqrt(v*v-(a.v[0])*(a.v[0])); // pithagoras
}

 // Five particles with similar initial conditions differing only in the direction of the initial velocity
void CollectionParticles(double x0, double x1, double G, double invM, double velocity)
{
	particle a1, a2, a3, a4, a5;
	control c;
	c.verbose = 1;
	c.printing = 5;
	double v = velocity; // absolute value
	double dt = 0.005 ; // time per step
	double t = 0.0 ; // start time
	int N = 8000 ;  // number of steps to make
	
	a1.im = a2.im = a3.im = a4.im = a5.im = invM; // initialising the particles
	a1.x[0] = a2.x[0] = a3.x[0] = a4.x[0] = a5.x[0] = x0;
	a1.x[1] = a2.x[1] = a3.x[1] = a4.x[1] = a5.x[1] = x1;
	a1.GMm = a2.GMm = a3.GMm = a4.GMm = a5.GMm = G;
	RandomVelocityDirection(v, a1, 123); // give a random direction of the velocities
	RandomVelocityDirection(v, a2, 321);
	RandomVelocityDirection(v, a3, 132);
	RandomVelocityDirection(v, a4, 231);
	RandomVelocityDirection(v, a5, 312);
	v2p(a1); v2p(a2); v2p(a3); v2p(a4); v2p(a5); // set up the momentum of each

	char f1[50]="output/file1.dat";
	char f2[50]="output/file2.dat";
	char f3[50]="output/file3.dat";
	char f4[50]="output/file4.dat";
	char f5[50]="output/file5.dat";


	LeapfrogToFile(a1, dt, t, N, c, f1);
	LeapfrogToFile(a2, dt, t, N, c, f2);
	LeapfrogToFile(a3, dt, t, N, c, f3);
	LeapfrogToFile(a4, dt, t, N, c, f4);
	LeapfrogToFile(a5, dt, t, N, c, f5);

}

 // Five particles with similar circular orbits with small perturberances in the momentum
void CollectionParticles2(double x0, double x1, double G, double invM, double dp)
{
	particle a1, a2, a3, a4, a5;
	control c;
	c.verbose = 1;
	c.printing = 5;
	double dt = 0.005 ; // time per step
	double t = 0.0 ; // start time
	int N = 8000 ;  // number of steps to make
	
	a1.im = a2.im = a3.im = a4.im = a5.im = invM; // initialising the particles
	a1.x[0] = a2.x[0] = a3.x[0] = a4.x[0] = a5.x[0] = x0;
	a1.x[1] = a2.x[1] = a3.x[1] = a4.x[1] = a5.x[1] = x1;

	a1.GMm = a2.GMm = a3.GMm = a4.GMm = a5.GMm = G;
	makeCircular(a1); // give each particle such speed that its orbit becomes circular
	makeCircular(a2);
	makeCircular(a3);
	makeCircular(a4);
	makeCircular(a5);
	v2p(a1); v2p(a2); v2p(a3); v2p(a4); v2p(a5); // set up the momentum of each

	a1.p[0]+=dp; //perturb the orbits
	a2.p[1]+=dp;
	a3.p[0]-=dp;
	a4.p[1]-=dp;
	a5.p[0]+=dp;
	a5.p[1]+=dp;

	char f1[50]="output/file1.dat";
	char f2[50]="output/file2.dat";
	char f3[50]="output/file3.dat";
	char f4[50]="output/file4.dat";
	char f5[50]="output/file5.dat";

	LeapfrogToFile(a1, dt, t, N, c, f1);
	LeapfrogToFile(a2, dt, t, N, c, f2);
	LeapfrogToFile(a3, dt, t, N, c, f3);
	LeapfrogToFile(a4, dt, t, N, c, f4);
	LeapfrogToFile(a5, dt, t, N, c, f5);

}


int main()
{
//  CollectionParticles(-1.0, -2.0, 1.0, 1.0, 0.8);
//  CollectionParticles2(-1.0, -2.0, 1.0, 1.0, 0.1);
  return 0;
}

