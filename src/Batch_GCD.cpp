#include <iostream>
#include <stdio.h>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include<vector>
#include<fstream>
#include<algorithm>
#include "Batch_GCD.hpp"
using namespace std;
typedef struct pkey pkey;

Batch_GCD::Batch_GCD(vector<mpz_class>set_of_numbers){
    my_numbers = set_of_numbers;
}
Batch_GCD::Batch_GCD(vector<string>set_of_strings){
    vector<string> vs2;
    vs2=nettoyer(set_of_strings);
    my_strings = vs2;
    my_numbers = commencer(vs2);
}
Batch_GCD::Batch_GCD(const char* txtname){
	ifstream infile(txtname);
	string line;
	vector<string> vs;
	vector<string> vs2;
	while(infile>>line){
        infile>>line;
		vs.push_back(line);
	}
	infile.close();
	vs2=nettoyer(vs);
	my_strings = vs2;
   	my_numbers =  commencer(vs2);
}
Batch_GCD::Batch_GCD(){

}
void Batch_GCD::choose_set_of_numbers(vector<string>set_of_strings){
	vector<string>vs2;
	vs2 = nettoyer(set_of_strings);
	my_strings = vs2;
	my_numbers = commencer(vs2);
}

//comparison functions for the struct pkey
bool compareString(const pkey &a, const pkey &b){//compares strings
		if(a.n.compare(b.n)<0)
			return true;
		else return false;
}
bool comparePos(const pkey &a, const pkey &b){//compares positions
		return a.pos<b.pos;
	}
vector<string> Batch_GCD::nettoyer (vector<string>& vs){
	pkey aux;
	vector<string> ordered;
	for(int i=0 ;i<vs.size(); i++){	
		aux.n=vs[i];
		aux.pos=i;
		aux.rep=false;
		propre.push_back((aux));
	}
	sort(propre.begin(), propre.end(),compareString);
	ordered.push_back(propre[0].n);
	for(int i=0; i<propre.size()-1; i++){
		if(propre[i].n.compare(propre[i+1].n)==0)
			propre[i+1].rep=true;
		else ordered.push_back(propre[i+1].n);
	}
	return ordered;	
}
vector<mpz_class> Batch_GCD::commencer (vector<string>& vs2){
	vector <mpz_class> vec;
	mpz_class aux;
	for(int i=0; i<vs2.size();i++){
		aux.set_str(vs2[i],10);
		vec.push_back(aux);
	}
	return vec;
}
void Batch_GCD::saveResults (vector<mpz_class>& results){
	int j=0;
	for(int i=0; i<results.size();i++){
		propre[j].result=results[i];
		j++;
		while(propre[j].rep){
			propre[j].result=results[i];
			j++;
		}
	}
}
void Batch_GCD::printResults(){
	 sort(propre.begin(),propre.end(),comparePos);
	 ofstream myfile;
	 myfile.open ("Results.txt");
	 mpz_class temp;
	 for (int i=0; i<propre.size(); i++){
		if(propre[i].result > 1){
			temp.set_str(propre[i].n, 10);
			myfile<<propre[i].pos;
			myfile<<" : ";
			myfile<<propre[i].n;
			myfile<<" = ";
			myfile<<propre[i].result;
			myfile<<" * ";
			myfile<<(temp/propre[i].result);
			myfile<<"\n";
		}

	}
	 myfile.close();
}
vector<mpz_class>Batch_GCD::get_my_numbers(){
    return my_numbers;
}
vector<string>Batch_GCD::get_my_strings(){
    return my_strings;
}
mpz_class Batch_GCD::product(){
    int k = my_numbers.size();
    if(k==0)
        return 1;
    else if(k==1)
        return my_numbers[0];
    else{
        vector<mpz_class>copy_of_my_numbers = my_numbers;
        while(k>1){
            for(int i=0; i<k; i+=2){
                if(i!=(k-1))
                    copy_of_my_numbers[i/2] = copy_of_my_numbers[i]*copy_of_my_numbers[i+1];
                else
                    copy_of_my_numbers[i/2] = copy_of_my_numbers[i];                
            }
            k = (k+1)/2;
            copy_of_my_numbers.resize(k);
        }
        return copy_of_my_numbers[0];
    }
}
vector<vector<mpz_class> > Batch_GCD::product_tree(){
    vector<vector<mpz_class> >result;
    vector<mpz_class>copy_of_my_numbers = my_numbers;
    result.push_back(copy_of_my_numbers);
    int k = copy_of_my_numbers.size();
    while(k>1){
         for(int i=0; i<k; i+=2){
            if(i!=(k-1))
                copy_of_my_numbers[i/2] = copy_of_my_numbers[i]*copy_of_my_numbers[i+1]; 
            else
                copy_of_my_numbers[i/2] = copy_of_my_numbers[i];
        }
        k = (k+1)/2;
        copy_of_my_numbers.resize(k);    
        result.push_back(copy_of_my_numbers); 
    }
    return result;
}
void Batch_GCD::print_product_tree(vector<vector<mpz_class> >p_T){
    int k  =p_T.size();
    for(int i=0; i<k; i++){
        cout << "[ ";
        for(int j=0; j<p_T[i].size(); j++){
            cout << p_T[i][j] << " ";
        }
        cout << "]" << endl;
    }
    cout << endl;
}
// builds the product tree of my_numbers and then compute the remainders of M when divided by its squares
vector<mpz_class> Batch_GCD::remainders_square(mpz_class M){
    vector<vector<mpz_class> >p_T = product_tree(); // 
    vector<mpz_class> result(1, M);
    int k = p_T.size();
    int len = 1; // represents the current size of result = size of p_T[i]
    for(int i=k-1; i>=0; i--){
        len =  p_T[i].size();
        result.resize(len);
        for(int j=len-1; j>=0; j--){
            result[j] = result[j/2] % (p_T[i][j]*p_T[i][j]);   
        }
    }
    return result;
}
void Batch_GCD::print_remainders(vector<mpz_class> rM){
    int k = rM.size();
    cout << "[ ";
    for(int i=0; i<k; i++){
        cout << rM[i] << " ";
    }
    cout << "]" << endl;
}
mpz_class Batch_GCD::gcd2 (mpz_class c1, mpz_class c2){
	mpz_t e1;
	mpz_init(e1);
	mpz_t e2;
	mpz_init(e2);
	mpz_t r;
	mpz_init(r);
	mpz_set(e1,c1.get_mpz_t());
	mpz_set(e2,c2.get_mpz_t());
	mpz_gcd(r,e1,e2);
	mpz_class l(r);
	return l;
}
vector<mpz_class> Batch_GCD::compute_gcds(){
    vector<mpz_class>rems = remainders_square(product());
    vector<mpz_class>gcds;
    int k = my_numbers.size();
    for(int i=0; i<k; i++){
        gcds.push_back(gcd2(my_numbers[i], rems[i]/my_numbers[i]));
    }   
    return gcds;
}
vector<mpz_class> Batch_GCD::compute_gcds_from_given_M(mpz_class M){
    vector<mpz_class>rems = remainders_square(M);
    vector<mpz_class>gcds;
    int k = my_numbers.size();
    for(int i=0; i<k; i++){
        gcds.push_back(gcd2(my_numbers[i], rems[i]/my_numbers[i]));
    }   
    return gcds;
}
vector<mpz_class> Batch_GCD::gcds_treatment(vector<mpz_class>&gcds){
    vector<mpz_class>factorizations;
    int k = my_numbers.size();
    for(int i=0; i<k; i++){
        if(gcds[i]==1) // we can't break the key
            factorizations.push_back(0);
        else if(gcds[i] < my_numbers[i]){
            if (gcds[i]*gcds[i] <= my_numbers[i])
                factorizations.push_back(gcds[i]);
            else
                factorizations.push_back(my_numbers[i]/gcds[i]);
        }
        else{
            bool is_repeated=false;
            for(int j=0; j<k; j++){
                if(my_numbers[j]==my_numbers[i] && i!=j){
                    is_repeated = true;
                    break;
                }
            }
            if(is_repeated)
                factorizations.push_back(0);
            else{
                bool pushed = false;
                for(int j=0; j<k; j++){
                    mpz_class temp = gcd2(my_numbers[i], my_numbers[j]);
                    if(j!=i && temp>1){
                        if(temp*temp <= my_numbers[i])
                            factorizations.push_back(temp);
                        else
                            factorizations.push_back(my_numbers[i]/temp);
                        pushed = true;
                        break;
                    }
                }
                if(!pushed){
                    cout << "ERROR ERROR ERROR"<< endl;
                }
            }
            
        }
    }
    return factorizations;
}


