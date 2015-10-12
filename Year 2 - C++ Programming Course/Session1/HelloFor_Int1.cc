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
  int p;
  cout<<"Choose a lucky number:"<<endl;
  cin>>p;
  for (int i=1; i<=10; i++) {
    cout << i << " hello world" << endl ;
    if ( i == p ) {
      cout << "that was lucky!" << endl ;
    } else {
      cout << endl ;
    }
  }
}
