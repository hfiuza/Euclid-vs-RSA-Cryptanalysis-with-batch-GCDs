#include <cstddef>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>

//using namespace std;

int main(){
	mpz_t grande;
	mpz_init(grande);
	mpz_set_str(grande,"174811295983785623754328728733873272783873287337082048697",10);
	char str[50];
	mpz_get_str (str, 10, grande);
	std::cout<<str<<std::endl;

}
