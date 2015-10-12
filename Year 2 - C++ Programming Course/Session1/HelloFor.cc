// HelloFor.cc
//   usage:
//     HelloFor
//   features:
//     uses cout to write text
//     for loop
//     if 

#include <iostream>
#include <cstdlib>

using namespace std;

int main()
{
  for (int i=1; i<=10; i++) {
    cout << i << " hello world" << endl ;
    if ( i == 7 ) {
      cout << "that was lucky!" << endl ;
    } else {
      cout << endl ;
    }
  }
}
