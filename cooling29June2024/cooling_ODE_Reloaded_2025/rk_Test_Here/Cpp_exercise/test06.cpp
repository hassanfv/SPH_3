#include <iostream>

using namespace std;

struct Star
{
  int x;
  int y;
  int z;
  
  string print_name()
  {
    return name;
  }
  
  private:
    string name = "Vegas";
}; 


int main()
{

  Star obj1;
  
  obj1.x = 110;
  obj1.y = 123;
  
  cout << obj1.x << '\n';
  cout << obj1.print_name() << '\n';  
  return 0;
}
