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
  int p, q;
  cout<<"Choose a lucky number:"<<endl;
  cin>>p;
  cout<<"Now choose the number of cycles"<<endl;
  cin>>q;
  for (int i=1; i<=q; i++) {
    cout << i << " hello world" << endl ;
    if ( i == p ) {
      cout << "that was lucky!" << endl ;
    } else {
      cout << endl ;
    }
  }
}
