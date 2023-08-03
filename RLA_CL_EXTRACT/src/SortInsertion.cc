/*
 * SortInsertion.cpp
 *
 *  Created on: May 20, 2014
 *      Author: abdullah-al-mamun
 */

#include <vector>
#include <string>
#include "SortInsertion.h"
using namespace std;


int compareStr(string v, string w, int d)
{
	return v.substr(d).compare(w.substr(d));
}

void exch(vector<StrPacket> &a, int i, int j)
{
	StrPacket temp	= a[i];
	a[i]			= a[j];
	a[j]			= temp;
}

void sortInsertion(vector<StrPacket> &a, int lo, int hi, int d)
{
	//cout << lo << ":" << hi << ":" << d);
	for (int i = lo; i <= hi; ++i)
	{
		for (int j = i; j > lo && compareStr(a[j].str, a[j-1].str, d) < 0; --j)
			exch(a, j, j - 1);
	}
}


