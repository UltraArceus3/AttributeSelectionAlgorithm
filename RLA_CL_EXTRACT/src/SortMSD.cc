/*
 * SortMSD.cpp
 *
 *  Created on: May 20, 2014
 *      Author: abdullah-al-mamun
 */

#include <iostream>
#include <vector>
#include <string>

#include "SortInsertion.h"
#include "SortMSD.h"

using namespace std;


const int R		= 40;
const int M 	= 5; // number of strings applicable for insertion sort
vector<StrPacket> aux;
char ALPHABET[256];

void initAlphabet();
int charAt(StrPacket s, int d);
void sortMSD2(vector<StrPacket> &a, int lo, int hi, int d);


 void initAlphabet()
{
	ALPHABET[97] = 0; // a
	ALPHABET[98] = 1;
	ALPHABET[99] = 2;
	ALPHABET[100] = 3;
	ALPHABET[101] = 4;
	ALPHABET[102] = 5;
	ALPHABET[103] = 6;
	ALPHABET[104] = 7;
	ALPHABET[105] = 8;
	ALPHABET[106] = 9;
	ALPHABET[107] = 10;
	ALPHABET[108] = 11;
	ALPHABET[109] = 12;
	ALPHABET[110] = 13;
	ALPHABET[111] = 14;
	ALPHABET[112] = 15;
	ALPHABET[113] = 16;
	ALPHABET[114] = 17;
	ALPHABET[115] = 18;
	ALPHABET[116] = 19;
	ALPHABET[117] = 20;
	ALPHABET[118] = 21;
	ALPHABET[119] = 22;
	ALPHABET[120] = 23;
	ALPHABET[121] = 24;
	ALPHABET[122] = 25; // z

	ALPHABET[48] = 26; // 0
	ALPHABET[49] = 27;
	ALPHABET[50] = 28;
	ALPHABET[51] = 29;
	ALPHABET[52] = 30;
	ALPHABET[53] = 31;
	ALPHABET[54] = 32;
	ALPHABET[55] = 33;
	ALPHABET[56] = 34;
	ALPHABET[57] = 35; // 9

	ALPHABET[35] = 36;
	ALPHABET[45] = 37;
	ALPHABET[46] = 38;
	ALPHABET[95] = 39;

}

 int charAt(StrPacket s, int d)
{
	if(d < s.str.length())
		return ALPHABET[(int) s.str.at(d)];
	else
		return -1;
}

vector<StrPacket> sortMSD(vector<StrPacket> a)
{
	//cout << "sortMSD()\n";
	initAlphabet();
	int N	= a.size();
	aux.resize(N);
	sortMSD2(a, 0, N - 1, 0);

	//for(int i = 0; i < a.size(); ++i)
		//cout << a[i].str << "\t";

	return a;
}

 void sortMSD2(vector<StrPacket> &a, int lo, int hi, int d)
{
	// cout << lo << ":" << hi << " d: " << d<< endl;
	if(hi <= lo + M)
	{
		sortInsertion(a, lo, hi, d);
		return;
	}

	int count[R + 2] = {0};
	for(int i = lo; i <= hi; ++i)
		count[charAt(a[i], d) + 2]++;

	for(int r = 0; r < R + 1; ++r)
		count[r + 1]	+= count[r];

	for(int i = lo; i <= hi; ++i)
		aux[count[charAt(a[i], d) + 1]++]	= a[i];

	for(int i = lo; i <= hi; ++i)
		a[i]	= aux[i - lo];

	for(int r = 0; r < R; ++r)
		sortMSD2(a, lo + count[r], lo + count[r + 1] - 1, d + 1);
}
