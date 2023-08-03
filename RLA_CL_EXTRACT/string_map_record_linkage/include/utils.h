//
// Created by Ahmed Soliman on 2/28/22.
//

#ifndef SL_UTILS_H
#define SL_UTILS_H

#include <vector>
#include <string>

std::vector<std::string> read_vector_string_from_file(const std::string& file_name, size_t max_lines = SIZE_MAX);
std::string time_stamp();
void parse_arguments(int argc, char **argv);

#endif //SL_UTILS_H
