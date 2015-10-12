#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

float MySqrt(float x);

int main()
{
 float number;
 cin>>number;
 cout<<MySqrt(number)<<endl;

 return 0; 
}

float MySqrt(float c)
{
 float lower=0.1, upper=10.1, root, sign;

 for (int i=0; i<21; i++) {
  root = (lower + upper)/2;
  sign = (c - root*root)*(c - lower*lower);
  if (sign < 0) upper = root; else lower = root;
 }
 
 return root;
}
