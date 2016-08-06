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

/*
  These functions make the conversion between a vector of strings representing integers into a single string as well as the reverse operation
  A vector [1,2,43] is associated to the string 1a2a43
*/


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

int main(){
	vector<string> vec;
	vec.push_back("10");
	vec.push_back("100");
	vec.push_back("1000");
	vec.push_back("10000");
	string aux;
	aux=vecToString(vec);
	cout<<aux<<endl;
	vector<string>vec2;
	vec2=stringToVec(aux);
	for(int i=0;i<vec2.size();i++){
		cout<<vec2[i]<<endl;
	}
}
