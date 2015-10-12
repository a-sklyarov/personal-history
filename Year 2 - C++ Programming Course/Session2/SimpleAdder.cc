// Program to add two integers typed by user at keyboard
#include <iostream>
using namespace std;

int main ()
{
	int a, total;

	cout << "Enter number of years:" << endl;
	cin >> a;
	total = static_cast<float>(365.25*24*3600*a);
	cout << "This is equal to " << total <<" seconds :P"<< endl;

	return 0;
}
