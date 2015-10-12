// Program to add two integers typed by user at keyboard
#include <iostream>
using namespace std;

int main ()
{
	int a, b, total;

	cout << "Enter integers to be added:" << endl;
	cin >> a >> b;
	total = a + 2*b;
	cout << "The sum is " << total << endl;

	return 0;
}
