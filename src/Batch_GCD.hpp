#include <iostream>
#include <stdio.h>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include<vector>
#include<fstream>
#include<algorithm>
using namespace std;

typedef struct pkey pkey;

struct pkey{
	string n;
	int pos;
	bool rep;
	mpz_class result;
};

class Batch_GCD {
  //// public methods&attributes can be accessed directly from calling methods
 public:
  // constructor
  Batch_GCD(vector<mpz_class>set_of_numbers);
  // constructor from a vector of strings
  Batch_GCD(vector<string>set_of_strings);
  // contructor from a file name
  Batch_GCD(const char* txtname);
  // empty constructor
  Batch_GCD();
  // function aimed to be used only when combined to the empty constructor
  void choose_set_of_numbers(vector<string>set_of_strings);
  //removes from the strig vector repeated keys
  vector<string> nettoyer (vector<string>& vs);
  // reads the input and put it in a vector of strings
  vector<mpz_class> commencer (vector<string>& vs2);
  //saves the results in the vector propre
  void saveResults (vector<mpz_class>& results);
  //prints the results into a file named results
  void printResults();
  // gets the content in my_numbers
  vector<mpz_class> get_my_numbers();
  // gets the content in my_strings
  vector<string> get_my_strings();
  // computes product of its numbers
  mpz_class product(void);
  // returns the product tree of its numbers
  vector<vector<mpz_class> >product_tree(void);
  // prints a given product tree
  void print_product_tree(vector<vector<mpz_class> >p_T);
  // computes the remainders of an integer M when divided by its squared numbers  
  vector<mpz_class> remainders_square(mpz_class M);
  // print the remainders stored at rM and obtained by the method remainders_square
  void print_remainders(vector<mpz_class> rM);
  // computes the gcd between two objects of class mpz_class
  // apparently the machines from our classroom don't recognize the method "gcd()" for this class 
  mpz_class gcd2 (mpz_class c1, mpz_class c2); 
  // computes the gcds G_i = gcd (N_i, M_i/N_i) for each number in my_numbers
  vector<mpz_class> compute_gcds(void);
  // computes the gcds G_i = gcd (N_i, M/N_i) for each number in my_numbers and a given number M
  // will be very important in the parallel implementation
  vector<mpz_class> compute_gcds_from_given_M(mpz_class M);
  // treats the gcds obained from compute_gcds, returning the smallest prime factor that we found for each element
  // if it couldn't break the key, it returns 0 for the element
  vector<mpz_class> gcds_treatment(vector<mpz_class>&gcds);

 private:
  // the set of keys we're analyzing
  vector<mpz_class>my_numbers;
  vector<string>my_strings;
  vector<pkey> propre;
};
