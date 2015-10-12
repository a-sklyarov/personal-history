// Here begind file NRoot.cc

#include <iostream>
using namespace std;

float NthRoot(float c, int n); // finds nth root of c
float power(float x, int n); // raises x to power of n

int main()
{
 float number;
 int n;
 cout<<"Choose a number of which you want to find the Nth root, followed by N: "<<endl;
 cin>>number>>n; // user chooses the number and n
 cout<<NthRoot(number, n)<<endl;

 return 0; 
}

float NthRoot(float c, int n)
{
 float lower=0.1, upper=10.1, root, sign;

 for (int i=0; i<25; i++) { // 25 times to get accuracy of about 10^-5
  root = (lower + upper)/2;
  sign = (c - power(root, n))*(c - power(lower, n));
  if (sign < 0) upper = root; else lower = root;
 }
 
 return root;
}

float power(float x, int n)
{
 float result = 1.0;
 for (int i=0; i<n; i++) result = result*x;
 return result;
}
