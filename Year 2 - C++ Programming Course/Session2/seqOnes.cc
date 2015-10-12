#include <iostream>
using namespace std;

int main()
{
 int x;

 cin>>x;

 float total = static_cast<float>(x);

 for (int i=0; i<20; i++) {
  cout.precision(10);
  cout<<total<<", ";
  total = static_cast<float>(x)+1/total;
 }

 cout<<endl;

 return 0;
}
