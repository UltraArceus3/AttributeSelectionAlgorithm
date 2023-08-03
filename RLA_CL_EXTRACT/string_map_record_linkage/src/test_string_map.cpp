//
// Created by Ahmed Soliman on 2/27/22.
//
#define BOOST_LOG_DYN_LINK 1

#include "string_map.h"
#include "utils.h"
#include <vector>
#include <string>
#include <boost/log/trivial.hpp>


/* read strings from a user input file and map it using user
 * specified dimensionality */
string_map read_strings() {
    std::string filename(256, '\0');
    filename = "/home/nachiket/rla_cl_exact/data/ds11.1.fn";
    std::cout << "Please enter the path for the dataset containing strings to be mapped: " <<
    " (e.g. " << filename << ") : ";
    //std::cin >> filename;
    std::vector<std::string> strings = read_vector_string_from_file(filename);
    COUNT n_dim = 0;
    std::cout << "Please enter the number of dimensions: (e.g. 30) ";
    std::cin >> n_dim;
    auto sm = string_map(strings, original_edit_dist, n_dim);
    sm.map();
    return sm;
}

/* read strings from a user input file and map it using user
 * specified dimensionality */
string_map read_strings_unmapped() {
    std::string filename(256, '\0');
    filename = "/home/nachiket/rla_cl_exact/data/ds11.1.fn";
    std::cout << "Please enter the path for the dataset containing strings to be mapped: " <<
    " (e.g. " << filename << ") : ";
    //std::cin >> filename;
    std::vector<std::string> strings = read_vector_string_from_file(filename);
    COUNT n_dim = 0;
    std::cout << "Please enter the number of dimensions: (e.g. 30) ";
    std::cin >> n_dim;
    auto sm = string_map(strings, original_edit_dist, n_dim);
    //sm.map();
    return sm;
}


/* Run interactive summary statistics session on the input mappings and save the results
 * into a user specified file */
void summarize(string_map sm) {
    std::string filename(256, '\0');
    std::cout << "Please enter the prefix for the filename to be used for saving summary statistics: ";
    std::cin >> filename;
    filename.append("_").append(time_stamp()).append(".summary");
    size_t n_rp=1000000;
    std::cout << "Number of random record pairs (e.g. 60 000 000) : ";
    std::cin >> n_rp;
    sm.save_summary(filename, n_rp, true);
}


/* save the input mapping into a user specified file */
void save_mappings(string_map sm) {
    std::string filename(256, '\0');
    std::cout << "Please enter the prefix for the filename to be used for saving this mapping: ";
    std::cin >> filename;
    filename.append("_").append(time_stamp()).append(".strmap");
    sm.save_mappings(filename);
}

/* load pre-computed mappings from a user input file */
string_map load_mappings(){
    std::string filename(256, '\0');
    std::cout << "Please enter the path for a string mapping file (*.strmap): ";
    std::cin >> filename;
    auto sm = string_map(filename);
    return sm;
}

void display_point(Point p) {
    for(auto &it : p){
        std::cout << it << " ";
    }
    std::cout << std::endl;
}


int main() {
    BOOST_LOG_TRIVIAL(info) <<"Testing the string_map algorithm";
    std::string str(256, '\0');
    Point p;
    char choice='x';
    string_map sm;
    bool exit = false;
    do {
        std::cout << "(r)ead (m) read mapping (w)rite (s)ummarize (l)oad (g)et coordinates e(x)it : ";
        std::cin >> choice;
        switch (choice) {
            case 'r':
            case 'R':
                sm = read_strings();
                break;
            case 'm':
                sm = read_strings_unmapped();
                break;
            case 'w':
            case 'W':
                assert(sm.is_mapped());
                save_mappings(sm);
                break;
            case 's':
            case 'S':
                //assert(sm.is_mapped());
                summarize(sm);
                break;
            case 'l':
            case 'L':
                sm = load_mappings();
                break;
            case 'g':
            case 'G':
                std::cout << "Enter a string: ";
                std::cin >> str;
                p = sm.get_coordinates(str);
                display_point(p);
                break;
            case 'x':
            case 'X':
                exit = true;

        }
    } while(!exit);
    return 0;
}
