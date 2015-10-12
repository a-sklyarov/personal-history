#include <iostream>
using namespace std;

int main()
{
 // i is the variable which will be used to find consecutive prime numbers
 // S is the sum we are aiming for
 // N is the counter for how many prime numbers we have used
 // sum is the sum of the prime numbers we have used
 // prime is the boolean variable which will be used to check if a number is prime or not

 int i=2, S, N=0, sum=0; 
 bool prime;

 cout<<"Choose a number S:"<<endl; // the user chooses the aim S 
 cin>>S;

 cout<<"Numbers are: "; // start to write all the needed prime numbers

 while (sum<S) {
  prime = 1;
  for (int j=2; j*j <= i; j++) { // here I have used a mathematical property which makes the algorithm more efficient
   if (i % j == 0) {prime=0; break;} // we stop the loop once we have found that the number is not prime
  }
  if (prime) {cout<<i<<" "; sum+=i;N++;}
  i++;
 }

 cout<<endl<<"We need the first "<<N<<" prime numbers"<<endl;
 cout<<"Their sum is "<<sum<<endl;

 return 0;
}
