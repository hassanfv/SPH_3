#include <iostream>
#include <string>

using namespace std;

class Student
{
  public:
    //Constructor - This function is called every time new object is created.
    Student()
    {
      cout << "An object created !" << endl;
    }
    
    //Deconstructor
    ~Student()
    {
      cout << "The object destructed !" << endl;
    }
    
};

int main()
{

  Student student1;
  
  int x = 4;
  
  cout << "x = " << x << endl;

}
