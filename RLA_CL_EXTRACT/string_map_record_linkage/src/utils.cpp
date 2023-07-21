//
// Created by Ahmed Soliman on 2/28/22.
//

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <ctime>
#include <iomanip>
#include <boost/log/trivial.hpp>


/* Read column of strings into a vector string up to <max_lines> */
std::vector<std::string> read_vector_string_from_file(const std::string& file_name, size_t max_lines = SIZE_MAX) {
    std::vector<std::string> strings;
    std::fstream fin;
    fin.open(file_name, std::ios::in);
    size_t n_lines = 0;
    while(fin.good()) {
        std::string str;
        fin >> str;
        if (str.size()>0) {
            std::string lower_str = boost::algorithm::to_lower_copy(str);
            strings.push_back(lower_str);
            n_lines++;
            if(n_lines>=max_lines) break;
        }
    }
    fin.close();
    BOOST_LOG_TRIVIAL(info) << "Read " << n_lines << " lines from file: " << file_name;
    return strings;
}

/* Return a timestamp string with the following format
 * YYYYMMDDHHmmss. */
std::string time_stamp(){
    time_t now = time(nullptr);
    tm ltm{};
    localtime_r(&now, &ltm);
    std::stringstream str;
    str << std::setfill('0')
                  << std::setw(4) << ltm.tm_year + 1900 << "_"
                  << std::setw(2) << ltm.tm_mon + 1 << "_"
                  << std::setw(2) << ltm.tm_mday << "_"
                  << std::setw(2) << ltm.tm_hour << "_"
                  << std::setw(2) << ltm.tm_min << "_"
                  << std::setw(2) << ltm.tm_sec;
    return str.str();
}

void parse_arguments(int argc, char **argv) {
    BOOST_LOG_TRIVIAL(debug) <<"argc = "<<argc;
    for(int i=0;i<argc;++i) {
        BOOST_LOG_TRIVIAL(debug) <<"argv["<<i<<"]= "<<argv[i];
    }
}
