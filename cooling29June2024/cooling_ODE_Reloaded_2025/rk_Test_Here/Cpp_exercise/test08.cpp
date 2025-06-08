#include <iostream>
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
    Customer(string first_name, string last_name, int age, float balance)
    {
      cout << "A new customer created ! " << endl;
      m_first_name = first_name;
      m_last_name = last_name;
      m_age = age;
      m_balance = balance;
    }
    
    void add_deposit(float deposit)
    {
      m_balance += deposit;
    }
    
    void print_account()
    {
      cout << "Customer: " << m_first_name << " " << m_last_name << endl;
      cout << "balance: " << m_balance << " Euro" << endl;
    }
  
};



int main()
{

  Customer customer1("Hassan", "Fathivavsari", 44, 3655.0);
  customer1.print_account();
  cout << "\n\n" << endl;
  customer1.add_deposit(500.0);
  customer1.print_account();

}




