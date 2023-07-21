
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
#include <vector>
#include <string>
#include <set>

using namespace std;


set<string> generateKmers(string& str,int k);
float calculateBasicQgram(string& str1, string& str2,int threshRem);
double calculateBasicHausdorffDistance(string& str1,string& str2,int threshRem);
int calculateHammingDistance(string& str1,string& str2);
 
