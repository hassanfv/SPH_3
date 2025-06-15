#include <iostream>


using namespace std;



void add(int &x, int &y, int &z)
{

  z = x + y;

}


typedef void addx(int&, int&, int&);

void do_op2(int &x, int&y, int &z, addx *func)
{

  func(x, y, z);

}


int main()
{

  int a = 3;
  int b = 5;
  int c = 4;
  
  do_op2(a, b, c, add);
  
  cout << "c = " << c << endl;

}



