#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

#define ranf() \
  ((double)random()/(1.0+(double)RAND_MAX)) // Uniform from interval [0,1)

int main(int argc, char* argv[])
{
  int    outcome, N=100000, count_in=0 , seed=1234 ; //N=100000 to get 2 decimal places
  double fraction_in ;

  if(argc>1)   {
    sscanf( argv[1], "%d", &N ) ; // put the first command-line argument in N
  }
  if(argc>2) {
    sscanf( argv[2], "%d", &seed ) ; // put the 2nd argument in seed
  }
  // Write out a summary of parameters
  cout << "# " << argv[0] << " N=" << N << " seed=" << seed << endl ;
  
  // Initialise random number generator 
  srandom(seed);
  
  // Perform N experiments. 
  for(int n=1; n<=N; n++) {
      double x = 1 + ranf(); // x is in [1,2)
      double y = ranf(); // y is in [0,1)
      outcome = ( x*y < 1.0 ) ? 1 : 0 ; 
      count_in += outcome ; 
  }

  fraction_in = static_cast<double>(count_in)/N;

  cout<<count_in<<"\t"<<N<<"\t"<<fraction_in<<endl;
// gives the number of points below the graph, the total number of points and the estimation of ln(2)

  return 0;
}

