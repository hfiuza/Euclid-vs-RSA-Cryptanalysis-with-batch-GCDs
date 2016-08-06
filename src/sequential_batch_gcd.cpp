#include <iostream>
#include <vector>
#include <stdio.h>
#include "Batch_GCD.hpp"

int main(){

   	const char* txtname="gen512-3233331101-0.dat";  
	Batch_GCD my_batch(txtname);
    	vector<mpz_class>gcds = my_batch.compute_gcds();
    	vector<mpz_class>factors = my_batch.gcds_treatment(gcds);
    	my_batch.saveResults(factors);
    	my_batch.printResults();
 }

