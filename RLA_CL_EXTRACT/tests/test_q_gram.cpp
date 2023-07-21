

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <ctype.h>
#include <math.h>
#include <utility>
#include <iterator>
#include "test_q_gram.h"

set<string> generateKmers(string& str,int k){

	set<string> k_mer_list;

	for(int i = 0; i < (str.length()- (k-1));i++){
		string k_mer = str.substr(i,k);
		//cout << k_mer << endl;
		k_mer_list.insert(k_mer);
	}

	return k_mer_list;		

}

int power_mod(int base,int pow,int p){
	int res = 0;
	for (int i=2; i <= pow;i++){
	
		res += ((base * base) % p);	
	}
	return res % p;
}

std::vector<int> generateIntKmers(string& str,int k){

	std::vector<int> k_mer_list;
	int res = 0;
	int p = 97;
	/*
		The Integer value for first K- Mer
	*/
	
	for (int i= 0;i < k;i++){
		res += (int)str[i] * power_mod(26,(k-i+1),p);	
	}
	
	//cout << res << endl;
	k_mer_list.push_back(res);
	//cout << "Here:" << endl;
	
	int tmp_pw = pow(26,k-1);
	
	
	//cout << str << "String:" << endl;
	if (str.length() >= k+1) {
		for (int i=1; i <= (str.length()- k);i++){
			res = ((res - ((int)str[i-1] * tmp_pw) * 26 % p) + (int)str[(i-1) + k] );
			//cout << "Inside -i" << res<< str[i] <<endl;
			k_mer_list.push_back(res);
		}
	}
		
	return k_mer_list;
}




float calculateBasicQgram(string& str1, string& str2,int threshRem ){

	int len_str1,len_str_2;
	int k = 3;

	set<string> s1 = generateKmers(str1,k);
	set<string> s2 = generateKmers(str2,k);

	set<string> intersect;
	set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                inserter(intersect, intersect.begin()));
	
	set<string>::iterator itr;
    
  	// Displaying set elements
  	cout << endl;
  	cout << "Intersection:" << endl;
  	for (itr = intersect.begin();itr != intersect.end(); itr++) 
  	{
    		cout << *itr << " ";
  	} 
  	cout << endl;
	
	set<string> union_set;
	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
                 inserter(union_set, union_set.begin()));
	
	cout << "Union:" << endl;
  	for (itr = union_set.begin();itr != union_set.end(); itr++) 
  	{
    		cout << *itr << " ";
  	} 
	cout << endl;
	
	float q_gram_distance = (float)intersect.size() / (float)union_set.size();

	return q_gram_distance;

}

int calculateHammingDistance(string& str1,string& str2){
	int hammingDistance = 0;
	for (int i = 0; i < str1.length();i++){
		if (i > str2.length()){
		    break;
		}
		else{
		    if (str1[i] != str2[i]){
		    	hammingDistance += 1;	
		    }
		}
	}
	return hammingDistance;
}



float calculateBasicQgramInt(string& str1, string& str2,int threshRem ){

	int len_str1,len_str_2;
	int k = 3;

	set<int> s1 = generateIntKmers(str1,k);
	set<int> s2 = generateIntKmers(str2,k);

	set<int> intersect;
	set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                inserter(intersect, intersect.begin()));
	
	set<int>::iterator itr;
    
  	// Displaying set elements
  	cout << endl;
  	cout << "Intersection:" << endl;
  	for (itr = intersect.begin();itr != intersect.end(); itr++) 
  	{
    		cout << *itr << " ";
  	} 
  	cout << endl;
	
	set<int> union_set;
	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
                 inserter(union_set, union_set.begin()));
	
	cout << "Union:" << endl;
  	for (itr = union_set.begin();itr != union_set.end(); itr++) 
  	{
    		cout << *itr << " ";
  	} 
	cout << endl;
	
	float q_gram_distance = (float)intersect.size() / (float)union_set.size();

	return q_gram_distance;

}

double calculateBasicQgramInt(string& str1, string& str2,int threshRem ){

	int len_str1,len_str_2;
	int k ;
	int min_length_str = min(str1.length(),str2.length());
	//k = (min_length_str) * 0.35;	
	
	
	if (min_length_str < 4) 
		k = (min_length_str) * 0.35;
	else
		k = (min_length_str) * 0.75;

	std::vector<int> s1 = generateIntKmers(str1,k);
	std::vector<int> s2 = generateIntKmers(str2,k);
	
	std::sort(s1.begin(), s1.end());
    	std::sort(s2.begin(), s2.end());
	
	//s1.erase( unique( s1.begin(), s1.end() ), s1.end() );
	//s2.erase( unique( s2.begin(), s2.end() ), s2.end() );
	
	/*
	std::vector<char> v1(s1.size());
    	std::copy(s1.begin(), s1.end(), v1.begin());
	
	std::vector<char> s2(s2.size());
    	std::copy(s2.begin(), s2.end(), s2.begin());
	
	std::sort(v1.begin(), v1.end());
    	std::sort(v2.begin(), v2.end());
	*/
	
	vector<int> intersect;
	set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                inserter(intersect, intersect.begin()));
	

	double q_gram_distance = (double)intersect.size() / (double)min_length_str;
	//cout << "Q-gram Distance" << q_gram_distance << endl;
	return q_gram_distance;

}





double calculateBasicHausdorffDistance(string& str1,string& str2,int threshRem){
	
	int k = threshRem;
	std::set<std::string> s1 = generateKmers(str1,k);
	std::set<std::string> s2 = generateKmers(str2,k);
	
	set<std::string>::iterator itr,itr2;
	vector< int > min_distances;
	
	for (itr = s1.begin(); itr != s1.end(); itr++) 
	{
		vector<int> distance;
		 
		for (itr2 = s2.begin(); itr2 != s2.end(); itr2++){
			
			std::string tmp1 = *itr;
			std::string tmp2 = *itr2;
			cout << tmp1 << endl;
			cout << tmp2 << endl;
			cout << "Hamming distance" << calculateHammingDistance(tmp1,tmp2);
			
			distance.push_back(calculateHammingDistance(tmp1,tmp2));		
		}
		cout << "\nMin Element = "<< *min_element(distance.begin(), distance.end());
		min_distances.push_back(*min_element(distance.begin(), distance.end()));
	}
	int size_difference = (s1.size() - s2.size()) ; 
	
	if (size_difference < 0){
		size_difference = - size_difference;
	}
	
	
	return (double) (*max_element(min_distances.begin(), min_distances.end()) + size_difference);

}



int main() {
	
	// Test Case 1
	/*
	std::string s1 = "John";
	std::string s2 = "Johnson";
	
	set<string> output = generateKmers(s2,3);	
	set<string>::iterator itr;
    
  	cout << "Displaying K-mer:" << endl;
  	for (itr = output.begin();itr != output.end(); itr++) 
  	{
    		cout << *itr << " ";
  	} 
	float q_gram = calculateBasicQgram(s1,s2,2);
	std::cout <<"Q_Gram Distance"<< q_gram;
	double haus_dorff = calculateBasicHammingDistance(s1,s2,3);
	std::cout <<"Hausdorff Distance"<< haus_dorff;
	
	// Test Case 2 
	std::string s3 = "Nachiket";
	std::string s4 = "Ahmed";
	
	q_gram = calculateBasicQgram(s3,s4,2);
	std::cout <<"Q_Gram Distance"<< q_gram;
	haus_dorff = calculateBasicHammingDistance(s3,s4,3);
	std::cout <<"Hausdorff Distance"<< haus_dorff;
	
	
	
	// Test Case 3 
	std::string s5 = "Ahmedullah";
	std::string s6 = "Ahmed";
	
	set<string> output_3 = generateKmers(s5,3);	
    
  	cout << "Displaying K-mer:" << endl;
  	for (itr = output_3.begin();itr != output_3.end(); itr++) 
  	{
    		cout << *itr << " ";
  	} 
	q_gram = calculateBasicQgram(s5,s6,2);
	std::cout <<"Q_Gram Distance"<< q_gram;
	haus_dorff = calculateBasicHammingDistance(s5,s6,3);
	std::cout <<"Hausdorff Distance"<< haus_dorff;	
	*/
	
	set<int>::iterator itr;
	set<int>::iterator itr_2;
	
	std::string s1 = "John";
	std::string s2 = "Johnson";	
	std::set<int> k_mer_list = generateIntKmers(s1,3);
	std::set<int> k_mer_list_j = generateIntKmers(s2,3);	
	
	cout << "Displaying K-mer - 1:" << endl;
  	for (itr = k_mer_list.begin();itr != k_mer_list.end(); itr++) 
  	{
    		cout << *itr << " ";
  	}
  	
  	cout << "Displaying K-mer - 2:" << endl;
  	for (itr_2 = k_mer_list_j.begin();itr_2 != k_mer_list_j.end(); itr_2++) 
  	{
    		cout << *itr_2 << " ";
  	} 
	
	float data = calculateBasicQgramInt(s1,s2,3);	
	cout << "Q-Gram" << data;		
	return 1; 
} 
