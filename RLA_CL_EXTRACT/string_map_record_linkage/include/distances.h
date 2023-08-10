//
// Created by ahs16102 on 2/27/22.
//

#ifndef SL_DISTANCES_H
#define SL_DISTANCES_H

#include <vector>
#include <iostream>
#include <set>
#include <string>


//Distance type
using Distance = double;

//Point type
using Point = std::vector<double>;

int original_edit_dist(std::string &str1, std::string &str2);


class distances {

public:
    static const int8_t THRESHOLD_MAX = INT8_MAX;
    static Distance euclid_dist(Point a, Point b);
    static Distance squared_euclid_dist(Point a, Point b);
    static int edit_dist(std::string &str1, std::string &str2, int thresholdRem) ;
    static Distance q_gram_distance(std::string &str1, std::string &str2);    
    //static std::set<std::string> generateKmers(std::string& str,int k);
};


#endif //SL_DISTANCES_H
