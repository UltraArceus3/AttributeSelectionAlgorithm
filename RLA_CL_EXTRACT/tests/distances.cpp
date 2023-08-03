//
// Created by Ahmed Soliman on 2/27/22.
//

#include "distances.h"
#include <cmath>
#include <cassert>
#include <string>
#include <vector>

#include <algorithm>
//Returns the Euclidean distance between two points
Distance distances::euclid_dist(Point a, Point b) {
    Distance d = 0;
    assert(a.size()==b.size());
    for(auto it_a=a.begin(), it_b=b.begin();it_a!=a.end();it_a++,it_b++) {
        d += pow((*it_a - *it_b),2);
    }
    return(sqrt(d));
}

Distance distances::q_gram_distance(std::string &str1, std::string &str2){

    int k = 3;

	std::set<std::string> s1 = distances::generateKmers(str1,k);
	std::set<std::string> s2 = distances::generateKmers(str2,k);

	std::set<std::string> intersect;
	std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                inserter(intersect, intersect.begin()));

	std::set<std::string> union_set;
	std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
                 inserter(union_set, union_set.begin()));

	double q_gram_distance = (double)intersect.size() / (double)union_set.size();

	return q_gram_distance;

}

std::set<std::string> generateKmers(std::string& str,int k){

	std::set<std::string> k_mer_list;

	for(int i = 0; i < (str.length() - k);i++){
		std::string k_mer = str.substr(i,k);
		k_mer_list.insert(k_mer);
	}
	return k_mer_list;

}


//Returns the Squared Euclidean distance between two points
Distance distances::squared_euclid_dist(Point a, Point b) {
    Distance d = 0;
    assert(a.size()==b.size());
    for(auto it_a=a.begin(), it_b=b.begin();it_a!=a.end();it_a++,it_b++) {
        d += pow((*it_a - *it_b),2);
    }
    return d;
}

//Returns the basic edit distance between two strings
int original_edit_dist(std::string &str1, std::string &str2)
{
    std::vector<std::vector<int>> matArr;
    int row		= {uint16_t(str1.length() + 1)};
    int col 	= {uint16_t(str2.length() + 1)};
    matArr.resize(row);
    for(int i = 0; i < row; ++i)
        matArr.at(i).resize(col, 0);
    for(int i = 0; i < row; i++)
    {
        for(int j = 0; j < col; j++)
        {
            if (i == 0)
                matArr[i][j] = j;
            else if (j == 0)
                matArr[i][j] = i;
            else
            {
                if((int)str1.at(i-1) == (int)str2.at(j-1))
                    matArr[i][j]	= std::min(std::min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
                else
                    matArr[i][j] 	= std::min(std::min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
            }
        }
    }

    return (matArr[row-1][col-1]);
}

//Returns the basic edit distance between two strings
inline int basic_edit_dist(std::string &str1, std::string &str2, int threshRem)
{
    int row;
    int col;
    int i;
    std::vector<std::vector<int>> matArr;

    row		= {uint16_t(str1.length() + 1)};
    col 	= {uint16_t(str2.length() + 1)};

    matArr.resize(row);
    for(i = 0; i < row; ++i)
        matArr.at(i).resize(col, 0);

    for(i = 0; i < row; i++)
    {
        for(int j = 0; j < col; j++)
        {
            if (i == 0)
                matArr[i][j] = j;
            else if (j == 0)
                matArr[i][j] = i;
            else
            {
                if((int)str1.at(i-1) == (int)str2.at(j-1))
                    matArr[i][j]	= std::min(std::min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
                else
                    matArr[i][j] 	= std::min(std::min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
            }

            if((row - col) == (i - j) && (matArr[i][j] > threshRem))
            {
                return distances::THRESHOLD_MAX;
            }
        }
    }

    return (matArr[row-1][col-1]);
}

//Returns the edit distance between two strings
int distances::edit_dist(std::string &str1, std::string &str2, int thresholdRem)
{
    if(thresholdRem == 0)
        return distances::THRESHOLD_MAX;
    else if((2 * thresholdRem + 1) >= (int) std::max(str1.length(), str2.length()))
        return basic_edit_dist(str1, str2, thresholdRem);
    else
    {
        std::string s1;
        std::string s2;
        int row;
        int col;
        int diagonal;
        int i;
        int j;
        std::vector<std::vector<int>> matArr;

        if (str1.length() > str2.length())
        {
            s1 = str2;
            s2 = str1;
        }
        else
        {
            s1 = str1;
            s2 = str2;
        }

        row	 		= {uint16_t(s1.length()+1)};
        col 		= 2 * thresholdRem + 1;
        diagonal 	= {uint16_t(thresholdRem + s2.length() - s1.length())};

        matArr.resize(row);
        for(i = 0; i < row; ++i)
            matArr.at(i).resize(col, 0);

        for(i = 0; i < thresholdRem + 1; i++)
        {
            for(j = thresholdRem - i; j < col; j++)
            {
                if (i == 0)
                    matArr[i][j]	= j - thresholdRem;
                else if(j == (thresholdRem - i))
                    matArr[i][j] 	= matArr[i - 1][j + 1] + 1;
                else if(j != (col - 1))
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j - (thresholdRem - i) - 1))
                        matArr[i][j]	= std::min(std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                    else
                        matArr[i][j] 	= std::min(std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                }
                else
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j - (thresholdRem - i) - 1))
                        matArr[i][j]	= std::min(matArr[i - 1][j], matArr[i][j - 1] + 1);
                    else
                        matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
                }

                if((j == diagonal) && matArr[i][j] > thresholdRem)
                    return distances::THRESHOLD_MAX;
            }
        }

        for(i = thresholdRem + 1; i < (int) s2.length() - thresholdRem + 1; i++)
        {
            for(j = 0; j < col; j++)
            {
                if(j == 0)
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j + (i - thresholdRem) - 1))
                        matArr[i][j]	= std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
                    else
                        matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
                }
                else if(j != (col - 1))
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j + (i - thresholdRem) - 1))
                        matArr[i][j] 	= std::min(std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                    else
                        matArr[i][j] 	= std::min(std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                }
                else
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j + (i - thresholdRem) - 1))
                        matArr[i][j] 	= std::min(matArr[i - 1][j], matArr[i][j - 1] + 1);
                    else
                        matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
                }
                if((j == diagonal) && (matArr[i][j] > thresholdRem))
                    return distances::THRESHOLD_MAX;
            }
        }

        for(i = {uint16_t(s2.length() - thresholdRem + 1)}; i < row; i++)
        {
            for(j = 0; j < col - i + (int) s2.length() - thresholdRem; j++)
            {
                if(j == 0)
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j + (i - thresholdRem) - 1))
                        matArr[i][j]	= std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
                    else
                        matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
                }
                else
                {
                    if((int)s1.at(i - 1) == (int)s2.at(j + (i - thresholdRem) - 1))
                        matArr[i][j] 	= std::min(std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                    else
                        matArr[i][j] 	= std::min(std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                }
                if((j == diagonal) && (matArr[i][j] > thresholdRem))
                    return distances::THRESHOLD_MAX;
            }
        }

        return matArr[row - 1][diagonal];

    }
}

