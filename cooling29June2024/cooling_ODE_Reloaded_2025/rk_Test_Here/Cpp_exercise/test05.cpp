#include <iostream>

using namespace std;

void swap_values(int &x, int &y)
{

  int tmp = x;
  x = y;
  y = tmp;

}


int main()
{

  int x = 5;
  int y = 10;
  
  swap_values(x, y);
  
  cout << x << ", " << y << endl;

}
