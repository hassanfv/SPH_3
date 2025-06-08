#include <iostream>
#include <string>
#include "customer.h"

using namespace std;


int main()
{

  Customer customer1("Hassan", "Fathivavsari", 44, 3655.0);
  customer1.print_account();
  cout << "\n\n" << endl;
  customer1.add_deposit(500.0);
  customer1.print_account();

}




