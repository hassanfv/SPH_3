#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int main()
{
  ifstream file ("data.txt");
  
  if (file.is_open())
  {
    cout << "File is open to read !";
  }
  else
  {
    cout << "file does not exist !";
  }
  cout << '\n';
  
  string line;
  
  getline(file, line);
  
  cout << line << '\n';
  
}
