#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main()
{

  ofstream file ("output2.txt");
  
  if (file.is_open())
  {
    cout << "file is open to write !";
  }
  else
  {
    cout << "File is not allowed for writting !";
  }
  cout << '\n';
  
  vector<string> vect = {"C", "E", "G", "B", "D"};
  
  for (int i = 0; i < vect.size(); i++)
  {
    file << i << "," << vect[i] << "\n";
  }
  file.close();

}
