/*
 * Cluster.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: abdullah-al-mamun
 */


// library
#include <iostream>
#include <vector>

#include "Cluster.h"


// namespace
using namespace std;

void Cluster::initCluster(float height, vector<string> item)
{
	this->height = height;
	itemArr.push_back(item);
}

void Cluster::initCluster(float height, Cluster c1, Cluster c2)
{
	this->height = height;
	itemArr = c1.itemArr;
	itemArr.insert(itemArr.end(), c2.itemArr.begin(), c2.itemArr.end());
}

vector<string> Cluster::getRepresentative()
{
	return itemArr.at(0);
}

void Cluster::printCluster()
{

}
