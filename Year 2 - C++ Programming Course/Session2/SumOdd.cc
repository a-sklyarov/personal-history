#include <iostream>
using namespace std;

int main()
{
 int N, sum=0;
 cin>>N;
 
 for (int i=1; i<=N; i++) sum += 2*i-1;
 cout<<sum<<", "<<N*N<<", "<<(sum==N*N)<<endl;

 return 0;
}
