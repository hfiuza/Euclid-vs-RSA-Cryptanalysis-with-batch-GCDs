#include <iostream>
#include "Batch_GCD.hpp"

int main(){
   	const char* txtname="gen512-3233331101-0.dat";  
    Batch_GCD my_batch(txtname);
    vector<mpz_class>gcds = my_batch.compute_gcds();
    vector<mpz_class>factors = my_batch.gcds_treatment(gcds);
    vector<mpz_class>numbers = my_batch.get_my_numbers();
    for(int i=0; i<factors.size(); i++){
        if(factors[i]!=0){
            cout << i << " " << factors[i] << " || "<<numbers[i]/factors[i] <<endl;
        }
    }
    cout << endl; 
}

