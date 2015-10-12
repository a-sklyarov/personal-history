#include <iostream>
#include <cmath>
using namespace std;

int main()
{
 int N;
 bool b;
 cin>>N;
 cout<<2<<endl;
 for (int i=3; i<=N; i++) {
  b=1;
  for (int j=2; j<(sqrt(i)+1); j++) {
   if (i%j == 0) {b=0; break;}
  }
  if (b==1) cout<<i<<endl;
 }

 return 0;
}
