#include <iostream>
#include <stdio.h>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include<vector>
#include<fstream>
#include <stdlib.h>
#include <time.h>
using namespace std;
 
vector<mpz_class> read_file(const char * txtname){  
    ifstream infile(txtname);
    string line;
    vector<string> vs;
    while(infile>>line){
        infile>>line;
        vs.push_back(line);
    }
    infile.close(); 
     
    vector <mpz_class> vec;
    mpz_class aux;
    for(int i=0; i<vs.size();i++){
        aux.set_str(vs[i],10);
        vec.push_back(aux);
    }
    return vec;
}
void print_keys(vector<mpz_class>&vec_keys, const char * txtname){
    ofstream outfile(txtname);
    int k = vec_keys.size();
    for(int i=0; i<k; i++)
        outfile << i << " " << vec_keys[i] << endl;
    outfile.close();
}
 
 
int main(){
    const char * txtname = "small_primes.txt";
    vector<mpz_class>vec_primes = read_file(txtname);
    vector<mpz_class>vec_keys;
    int k = vec_primes.size();
    int index1, index2, n_keys = 1000;
    int temp;
    for(int i=0; i<n_keys; i++){
	temp = rand();
        index1 = ((int)temp) % k;
	temp = rand();
        index2 = ((int)temp) % k;
	if(index1 < index2)
	        vec_keys.push_back(vec_primes[index1]*vec_primes[index2]);
    }
    print_keys(vec_keys, "keys.dat");    
     
}
