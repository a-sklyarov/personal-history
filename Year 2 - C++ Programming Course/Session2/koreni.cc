#include <iostream>
#include <cmath>
using namespace std;

int main()
{
 float total = 1.0;

 for (int i=0; i<20; i++) {
  cout.precision(10);
  cout<<total<<", ";
  total = sqrt(1.0+total);
 }

 cout<<endl;

 return 0;
}
