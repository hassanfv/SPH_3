#include "customer.h"
#include <iostream>

using namespace std;

// Constructor - this member function is called every time we create a new object.
Customer::Customer(string first_name, string last_name, int age, float balance)
{
  cout << "A new customer created ! " << endl;
  m_first_name = first_name;
  m_last_name = last_name;
  m_age = age;
  m_balance = balance;
}

void Customer::add_deposit(float deposit)
{
  m_balance += deposit;
}

void Customer::print_account()
{
  cout << "Customer: " << m_first_name << " " << m_last_name << endl;
  cout << "balance: " << m_balance << " Euro" << endl;
}
