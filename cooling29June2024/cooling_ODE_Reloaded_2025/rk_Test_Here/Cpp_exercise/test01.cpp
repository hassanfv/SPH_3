#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main()
{

  vector<string> names;
  names.push_back("Hassan");
  names.push_back("John");
  names.push_back("Peter");
  
  ofstream file ("output.txt");
  
  if (file.is_open())
  {
    cout << "output.txt successfully opened !" << '\n';
  }
  
  for (int i = 0; i < names.size(); i++)
  {
    file << names[i] << '\n';
  }
  
  file.close();
}



