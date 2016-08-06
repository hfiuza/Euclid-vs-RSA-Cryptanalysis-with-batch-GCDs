#include <iostream>
#include <vector>
//#include <string>
#include <stdio.h>
//#include <string.h>
#include "Batch_GCD.hpp"
#include "mpi.h"

string vecToString(vector<string> vs){
	string converted="";
	for(int i=0;i<vs.size();i++){
		converted=converted+vs[i];
		converted=converted+"a";
	}
	return converted;
}

vector<string> stringToVec(string str){
	vector<string> unconverted;
	int begin=-1;
	string aux;
	for(int i=0;i<str.size();i++){
		if(str.at(i)=='a'){
			aux=str.substr(begin+1,i-begin-1);
			begin=i;
			unconverted.push_back(aux);
		}
	}
	return unconverted;
}

char * stringToChar(string str){
	int len = str.length();
	char * s = new char[len+1];
	for(int i=0; i<len; i++)
		s[i] = str[i];
	s[len] = '\0';
	return  s;
}
string charToString(char * s){
	int len = strlen(s);
	string str = "";
	for(int i=0; i<len; i++)
		str += s[i];
	return str;
}


int main(int argc, char*argv[]){
	MPI_Init(&argc, &argv);

	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Request request[2];
	MPI_Status stats;
	vector<string> all_strings;
	string all;
	char * all2;
	int siz;
	Batch_GCD * batch_of_all;
	if(rank==0){
		// we read all numbers from a text file and we store the in a single string all (equivalent to the char * all2)
		const char * in_file_name = "gen512-3233331101-0.dat";
		ifstream in_file(in_file_name);
		string line;
		while(in_file>>line){
        		in_file>>line;
			all_strings.push_back(line);
		}
		in_file.close();
		batch_of_all = new Batch_GCD(all_strings);
		all_strings = batch_of_all->get_my_strings();
		// now all_strings is alphabetically sorted and doesn't contain repeated elements
		// we stress that the Batch_of_GCD object knows whether an element is repeated or not 
		all = vecToString(all_strings);
		all2 = stringToChar(all);
		siz = strlen(all2);
	}

	// we aim to broadcast the char * all2 that contains encrypts all numbers N_0, N_1, ..., N_{k-1}
	// so first we broadcast its size
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&siz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank!=0)
		all2 = new char[siz]; // here we reserve the necessary amount of memory at each processor
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(all2, siz, MPI_CHAR, 0, MPI_COMM_WORLD);
	all = charToString(all2);
	// this was the last call to all2. So we delete it to free memory
	delete[] all2;
	all_strings = stringToVec(all);
	// k is the total number of N_i's we're going to study
	int k = all_strings.size();	
	int q = k/nprocs;
	int r = k%nprocs;
	// the processors will contain respectively q+1, q+1, q+1, ... ,q+1, q, q, ..., q elements
	// so we build a vector of strings local_strings at each processor that contains the strings analyzed by this processor
	vector<string>::const_iterator first;
	vector<string>::const_iterator last;
	if(rank < r){
		first = all_strings.begin() + (q+1)*rank;
		last = all_strings.begin() + (q+1)*(rank+1);
	}
	else{
		first = all_strings.begin() + (q+1)*r + q*(rank-r);
		last = all_strings.begin()+(q+1)*r+q*(rank-r+1);
	}
	vector<string> local_strings(first, last);

	// we then create a Batch_GCD object containing the fraction of numbers studied by the processor

	Batch_GCD my_batch(local_strings);
	MPI_Barrier(MPI_COMM_WORLD);

	// our first operation is to calculate the local product and convert it to a string
	mpz_class local_product = my_batch.product();
	string local_product_string = local_product.get_str();
	char * local_product_char = stringToChar(local_product_string);
	int local_product_size = strlen(local_product_char);
	int * all_local_products_sizes;
	int * sum_of_all_local_products_sizes;
	if(rank==0)
		all_local_products_sizes = new int[nprocs];
	// We call MPI_Gather to gather the sizes of the local products strings in the root (0)
	MPI_Barrier(MPI_COMM_WORLD);


	MPI_Gather(&local_product_size, 1, MPI_INT, all_local_products_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	char * all_local_products_chars;
	string all_local_products_string;
	if(rank==0){
		// we then reserve sufficient memory at the root to receive the local products strings
		// this is some preprocessing required to use the function MPI_Gatherv
		sum_of_all_local_products_sizes = new int[nprocs];
		sum_of_all_local_products_sizes[0] = 0;
		for(int i=1; i<nprocs; i++)
			sum_of_all_local_products_sizes[i] = sum_of_all_local_products_sizes[i-1] + all_local_products_sizes[i-1];
		all_local_products_chars = new char [sum_of_all_local_products_sizes[nprocs-1]+all_local_products_sizes[nprocs-1] ];
	}
	// then we call MPI_Gatherv that allow us to gather arrays of different sizes	
	MPI_Barrier(MPI_COMM_WORLD);



 	MPI_Gatherv(local_product_char, local_product_size, MPI_CHAR, all_local_products_chars, all_local_products_sizes, sum_of_all_local_products_sizes, MPI_CHAR, 0, MPI_COMM_WORLD);	
	// M will contain the product of all numbers
	mpz_class M;
	int lenM;
	int all_local_products_chars_length;
	char * charM;
	string strM;
	if(rank==0){
		all_local_products_chars_length = strlen(all_local_products_chars);
		all_local_products_string = charToString(all_local_products_chars);
		vector<string>local_products_strings(nprocs);
		for(int i=0; i<nprocs; i++)
			local_products_strings[i] = all_local_products_string.substr(sum_of_all_local_products_sizes[i], all_local_products_sizes[i]);
		// we create a Batch_GCD object to calculate the product of all subproducts
		Batch_GCD temp_batch(local_products_strings);
		vector<mpz_class>numbers = temp_batch.get_my_numbers();
		M = temp_batch.product();
		//cout << " we got M = " << M << endl;
		strM = M.get_str();
		charM = stringToChar(strM);
		lenM = strlen(charM);	
	}
	// then we broadcast the string that represents M calculated in the processor 0 to all processors
	// first we broadcast its length
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&lenM, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	if(rank!=0) // then  we reserve enough memory space
		charM = new char[lenM];
	// then we call MPI_Bcast
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(charM, lenM, MPI_CHAR, 0, MPI_COMM_WORLD);
	strM = charToString(charM);
	M.set_str(strM, 10);
	// we call the special method from Batch_GCD class called compute_gcds_from_given_M	
	vector<mpz_class> gcds = my_batch.compute_gcds_from_given_M(M);
	
	
	// then we print for each processor the numbers whose corresponding gcd isn't neither 1 nor N_i
	vector<mpz_class> numbers = my_batch.get_my_numbers();
	
	// we will print all the keys we shall break
	vector<mpz_class> factors(numbers.size()); // will contain the smallest prime factor of numbers[i]
	for(int i=0; i<numbers.size(); i++){
		if(gcds[i]==1)
			factors[i]=1;
		else if(gcds[i]>1 && gcds[i] < numbers[i]){
			// in this case, gcds[i] is a prime factor of numbers[i]
			// store at factors the minimum of gcds[i]  and numbers[i]/gcds[i]
	        	//cout <<"rank: " << rank << ", " <<  i << " " << numbers[i] << " || "<<gcds[i] <<endl;	
			factors[i] = (gcds[i]*gcds[i] <= numbers[i]) ? gcds[i] : numbers[i]/gcds[i];
		}
        }
	//for(int i=0; i<numbers.size(); i++){
	//	cout << "at " << rank <<": " <<numbers[i]<<" "<<gcds[i] << endl; 
	//}

	// we treat the the cases when 1 < gcds[i] < numbers[i] first in all processors
	MPI_Barrier(MPI_COMM_WORLD);
	// the cases gcds[i] = numbers[i] can be much slower, so we only treat them after we finish the other cases

	bool entered_if_before = false;
	vector<string>local_products_strings;
	MPI_Bcast(&all_local_products_chars_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank!=0)
		all_local_products_chars = new char[all_local_products_chars_length];
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(all_local_products_chars, all_local_products_chars_length, MPI_CHAR, 0, MPI_COMM_WORLD);
	if(rank!=0){
		all_local_products_sizes = new int[nprocs];
		sum_of_all_local_products_sizes = new int[nprocs];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(all_local_products_sizes, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(sum_of_all_local_products_sizes, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank!=0){
		all_local_products_string = charToString(all_local_products_chars);
		local_products_strings.resize(nprocs);
		for(int k=0; k<nprocs; k++)
			local_products_strings[k] = all_local_products_string.substr(sum_of_all_local_products_sizes[k], all_local_products_sizes[k]);
	}
	//cout << "4 - hello from "<<rank << endl;

	for(int i=0; i<numbers.size(); i++){
		if(gcds[i] == numbers[i]){

		//cout << "here0.5 for numbers[i]="<<numbers[i] << endl;

			// this is the most complicated case
			// we need to compute the gcd of numbers[i] with every other N_i assigned to other processors
			// at first we check if there is a number from the current processor that has a common factor with numbers[i]
			mpz_class tempR = local_product % (numbers[i]*numbers[i]);
			mpz_class tempF = my_batch.gcd2(tempR/numbers[i], numbers[i]);
			mpz_class temp_temp;
			if(tempF==1);
			else if(tempF<numbers[i]){
		        	//cout <<"rank: " << rank << ", " <<  i << " " << numbers[i] << " || "<<tempF <<endl;				
				factors[i] = (tempF*tempF <= numbers[i]) ? tempF : numbers[i]/tempF;
				continue;
			}
			else{
				// in this case we need to compute gcd of numbers[i] with all numbers from this processor
				factors[i] = 1;
				for(int j=0 ; j<numbers.size(); j++){	
					if(j!=rank){
						temp_temp = my_batch.gcd2(numbers[i], numbers[j]);
						if(temp_temp>1){
							factors[i] = (temp_temp*temp_temp <= numbers[i]) ? temp_temp : numbers[i]/temp_temp;
							break;
						}
					}
				}
				if(factors[i]==1)
					cout << "ERROR ERROR ERROR ERROR 2" <<endl;
				continue;
			}

		//cout << "here1 for numbers[i]="<<numbers[i] << endl;


			// if no number in this processor has a common factor with numbers[i], then we need to compute gcd of numbers[i] with numbers from other processors
			// we first determine one processor that contains such a numbers[j]
			// we do so by calculating gcd2(local_product, numbers[i]) at each processor
			// if we get it to be greater than 1, then we know were to search numbers[j] with the desired property

			// now we are sure that local_products_strings[] has been initialized


			Batch_GCD batch_temp(local_products_strings);
			int j;
			mpz_class tempP;
			mpz_class tempGCD=1;
			for(j=0; j<nprocs; j++){
				if(j==i) continue;
				tempP.set_str(local_products_strings[j], 10);
				tempR = tempP%numbers[i];
				tempGCD = my_batch.gcd2(tempR, numbers[i]);
				if(tempGCD > 1)
					break; 
			}
			if(tempGCD==1){
				cout << "ERROR ERROR ERROR ERROR ERROR ERROR ERROR " << endl;
				break;
			}
			else if(tempGCD < numbers[i]){
				factors[i] = (tempGCD*tempGCD <= numbers[i]) ? tempGCD : numbers[i]/tempGCD;
				continue;
			}
			else{
				// in this case, tempGCD is equal to numbers[i]
				// so we know that the processor j contains an element that has a common factor with numbers[i]
				// these strings are stored at all_strings in a certain range
				// we do as we did in the beginning of the code
				if(j < r){
					first = all_strings.begin() + (q+1)*j;
					last = all_strings.begin() + (q+1)*(j+1);
				}
				else{
					first = all_strings.begin() + (q+1)*r + q*(j-r);
					last = all_strings.begin()+(q+1)*r+q*(j-r+1);
				}
				vector<string>searched_strings(first, last);		
				// in this case we need to compute gcd of numbers[i] with all numbers from the vector searched_strings
				factors[i] = 1;
				mpz_class temp_temp;
				mpz_class tempGCD=1;
				for(int k=0; k<searched_strings.size(); k++){
					temp_temp.set_str(searched_strings[k], 10);		
					tempGCD = my_batch.gcd2(numbers[i], temp_temp);
					if(tempGCD > 1){
						factors[i] = (tempGCD*tempGCD <= numbers[i]) ? tempGCD : numbers[i]/tempGCD;
						break;
					}
				}
				if(tempGCD==1)
					cout << "ERROR ERROR ERROR ERROR 2" <<endl;
			}
			//cout << "here3 for numbers[i]="<<numbers[i] << endl;


		}
	}
	
	for(int i=0; i<numbers.size(); i++){
		if(factors[i] != 1)
			cout <<"rank: " << rank << ", " <<  i << " " << numbers[i] << " = " <<factors[i] <<" * " << numbers[i]/factors[i] << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// now that we have computed factors, we gather these factors in the processor 0

	// the rest of the code consists of gathering factors[] values in the processor 0 and then printing numbers and factors in the original order
	// it might be slower than all the code above so we may consider to remove it	
	vector<string>factors_strings;
	vector<string>numbers_strings = my_batch.get_my_strings();
	for(int i=0; i<factors.size(); i++){
		string temp = factors[i].get_str();
		factors_strings.push_back(temp);
	}
	string local_factors_string = vecToString(factors_strings);
	char * local_factors_chars = stringToChar(local_factors_string);
	int local_factor_size = strlen(local_factors_chars);
	int * all_local_factors_strings_sizes;
	int * sum_of_all_local_factors_strings_sizes;
	//cout << "0 - hello from processor " << rank << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	/*for(int i=0; i<numbers.size(); i++){
		if(gcds[i]>1){
			if(numbers[i]%factors[i]!=0)
				cout << "ERROR: not a factor " << endl;
			else if(factors[i] == 1)
				cout << "ERROR: factors is equal to 1 while gcd is " << gcds[i] << "        and numbers is "<< numbers[i] << endl;
			else if(factors[i] >= numbers[i])
				cout << "ERROR: factors is greater or equal to gcd"<<endl;
		}

	}*/

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0)
		all_local_factors_strings_sizes = new int[nprocs];
	// We call MPI_Gather to gather the sizes of the local factors strings in the root (0)
	MPI_Gather(&local_factor_size, 1, MPI_INT, all_local_factors_strings_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	char * all_local_factors_chars;
	string all_local_factors_string;
	if(rank==0){
		// again, we then reserve sufficient memory at the root to receive the local factors strings
		// this is some preprocessing required to use the function MPI_Gatherv
		sum_of_all_local_factors_strings_sizes = new int[nprocs];
		sum_of_all_local_factors_strings_sizes[0] = 0;
		for(int i=1; i<nprocs; i++)
			sum_of_all_local_factors_strings_sizes[i] = sum_of_all_local_factors_strings_sizes[i-1] + all_local_factors_strings_sizes[i-1];
		all_local_factors_chars = new char[sum_of_all_local_factors_strings_sizes[nprocs-1] + all_local_products_sizes[nprocs-1] ];		
	}
	MPI_Gatherv(local_factors_chars, local_factor_size, MPI_CHAR, all_local_factors_chars, all_local_factors_strings_sizes, sum_of_all_local_factors_strings_sizes, MPI_CHAR, 0, MPI_COMM_WORLD);
	

	MPI_Barrier(MPI_COMM_WORLD);
	
	int all_local_factors_chars_length;
	if(rank==0){
		all_local_factors_chars_length = strlen(all_local_factors_chars);
		all_local_factors_string = charToString(all_local_factors_chars);
		vector<string>local_factors_strings(nprocs);

		for(int i=0; i<nprocs; i++)
			local_factors_strings[i] = all_local_factors_string.substr(sum_of_all_local_factors_strings_sizes[i], all_local_factors_strings_sizes[i]);
		vector<vector<string> >factors_strings_by_processor(nprocs, vector<string>());
		for(int i=0; i<nprocs; i++)
			factors_strings_by_processor[i] = stringToVec(local_factors_strings[i]);
		vector<mpz_class>all_factors;
		mpz_class tempFac;
		for(int i=0; i<nprocs; i++){
			for(int j=0; j<factors_strings_by_processor[i].size(); j++){
				tempFac.set_str(factors_strings_by_processor[i][j], 10);
				all_factors.push_back(tempFac);
			}
		}
		batch_of_all->saveResults(all_factors);
		batch_of_all->printResults();
	}

	MPI_Finalize();
	return 0;	


}



