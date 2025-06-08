#ifndef CUSTOMER_H
#define CUSTOMER_H

#include <string>

using namespace std;

class Customer
{

  public:
  
    string m_first_name;
    string m_last_name;
    int m_age;
    float m_balance;
  
    // Constructor - this member function is called every time we create a new object.
    Customer(string first_name, string last_name, int age, float balance);
    
    void add_deposit(float deposit);
    
    void print_account();
  
};


#endif
