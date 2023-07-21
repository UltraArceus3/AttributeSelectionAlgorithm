/*
 * Cluster.h
 *
 *  Created on: Apr 4, 2014
 *      Author: abdullah-al-mamun
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <vector>
using namespace std;

class Cluster
{
	public:
		float height;
		vector<vector<string> > itemArr;

		void initCluster(float height, vector<string> item);
		void initCluster(float height, Cluster c1, Cluster c2);
		vector<string> getRepresentative();
		void printCluster();
};

#endif



