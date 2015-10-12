// Program that produces all the numbers up to a given number (set by the user),
// raised to a power, again chosen by the user
#include <iostream>
using namespace std;

int main()
{
 int N, p; // N is the final number, p is the power

 cout<<"Choose a number:"<<endl;
 cin>>N;
 cout<<"Now choose a power:"<<endl;
 cin>>p;
 cout<<"Results:"<<endl;
 
 for (int i=1; i<=N; i++) {
  int result=1; //accumulation variable
  for (int j=0; j<p; j++) {
   result = result*i;
  }
  if (i==N) {cout<<result;} else {cout<<result<<", ";}
 }
 cout<<endl;

return 0;
}
