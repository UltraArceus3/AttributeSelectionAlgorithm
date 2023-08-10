//
// Created by Ahmed Soliman on 2/27/22.
//

#ifndef SL_STRING_MAP_H
#define SL_STRING_MAP_H

#include <vector>
#include <string>
#include "distances.h"
#include <random>
#include <iostream>
#include <chrono>
#include <ratio>
#include <functional>
#include <unordered_map>


using COUNT = unsigned long int;
using INDEX = unsigned long int;
using DIM = uint8_t;
using DistFuncCallback = std::function<int(std::string&, std::string&)>;
using milliseconds = std::chrono::duration<double, std::milli>;

/*The string_map class is used to map original objects (strings) into images in a d-dimensional Euclidean target space.
 * The dimension of the target space is user defined. */
class string_map {

private:
    std::vector<std::string> strings; // list of input strings to be mapped
    COUNT n_strings_input = 0;  // number of input strings
    COUNT n_strings = 0;  // number of strings after deduplication
    DIM n_dimensions = 0;  // dimensionality of the target space
    //pointer to a metric function that gets
    //original metric distance between two strings.
    DistFuncCallback distance_function;
    std::vector<std::vector<std::string>> pivots;  //pivot strings
    std::vector<Point> coord;  //image coordinates
    bool mapped = false;  //mapping status
    std::mt19937 generator;
    std::uniform_int_distribution<COUNT> distribution;
    milliseconds mapping_time = milliseconds::zero();
    std::unordered_map<std::string, Point> str2point;

    // Populate str2point map
    void populate_str2point();

public:
    //Constructor
    string_map(std::vector<std::string> const& strings,
               DistFuncCallback dist_function,
               DIM n_dimensions);

    //Load from file constructor
    explicit string_map(const std::string& file_name);

    int original_metric(INDEX ind1, INDEX ind2);

    //Get the dimensionality of the target space
    int get_n_dimensions() const;

    //Get the number of strings
    COUNT get_n_strings() const;

    //choose two distant strings as pivots.
    std::vector<COUNT> choose_pivots(int h);

    /* Get distance of two strings (indexed by a and b)
     * after projection on the first h-1 axes */
    Distance get_distance(COUNT a, COUNT b, int h);

    /* Get Euclidean distance of two strings (indexed by a and b)
     * based on the mapped coordinates */
    Distance get_euclidean_distance(COUNT a, COUNT b) const;

    /* Get Squared Euclidean distance of two strings (indexed by a and b)
     * based on the mapped coordinates */
    Distance get_squared_euclidean_distance(COUNT a, COUNT b) const;

    Distance get_q_gram_distance(std::string&  a, std::string&  b) const;
    

    // get the status of this mapper
    bool is_mapped() const;

    //Map the input strings (objects) into images in the target Euclidean space.
    void map();

    // get mapping time in milliseconds
    double get_mapping_time() const;

    // save mapping to file
    void save_mappings(const std::string& file_name);

    // Evaluate mappings
    void save_summary(const std::string &file_name, size_t n_record_pairs, bool squared_dist= true);

    // Empty constructor
    string_map();

    // Get coordinates for a string
    Point get_coordinates(const std::string& str);

};


#endif //SL_STRING_MAP_H
