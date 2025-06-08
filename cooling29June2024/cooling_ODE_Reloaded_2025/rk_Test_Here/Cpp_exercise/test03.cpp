#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;


int main()
{

  ifstream file ("data.txt");
  
  string input;
  
  vector<string> names;
  
  while (file >> input)
  {
    names.push_back(input);
  }

  for (string x : names)
  {
    cout << x << '\n';
  }

}
