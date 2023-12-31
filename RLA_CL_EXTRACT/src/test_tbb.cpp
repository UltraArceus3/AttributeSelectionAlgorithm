#include <iostream>
#include <string>
#include <vector>
#include <algorithm> // std::max

#include "tbb/parallel_for.h"
#include "oneapi/tbb/blocked_range.h"


bool silent = false;

static const std::size_t N = 23;

class SubStringFinder {
    const std::string &str;
    std::vector<std::size_t> &max_array;
    std::vector<std::size_t> &pos_array;

public:
    void operator()(const oneapi::tbb::blocked_range<std::size_t> &r) const {
        for (std::size_t i = r.begin(); i != r.end(); ++i) {
            std::size_t max_size = 0, max_pos = 0;
            for (std::size_t j = 0; j < str.size(); ++j) {
                if (j != i) {
                    std::size_t limit = str.size() - (std::max)(i, j);
                    for (std::size_t k = 0; k < limit; ++k) {
                        if (str[i + k] != str[j + k])
                            break;
                        if (k > max_size) {
                            max_size = k;
                            max_pos = j;
                        }
                    }
                }
            }
            max_array[i] = max_size;
            pos_array[i] = max_pos;
        }
    }

    SubStringFinder(const std::string &s, std::vector<std::size_t> &m, std::vector<std::size_t> &p)
            : str(s),
              max_array(m),
              pos_array(p) {}
};

int main(int argc, char *argv[]) {
    // command line parsing
    utility::parse_cli_arguments(argc,
                                 argv,
                                 utility::cli_argument_pack()
                                     //"-h" option for displaying help is present implicitly
                                     .arg(silent, "silent", "no output"));

    std::string str[N] = { std::string("a"), std::string("b") };
    for (std::size_t i = 2; i < N; ++i)
        str[i] = str[i - 1] + str[i - 2];
    std::string &to_scan = str[N - 1];
    const std::size_t num_elem = to_scan.size();

    std::vector<std::size_t> max(num_elem);
    std::vector<std::size_t> pos(num_elem);

    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, num_elem),
                              SubStringFinder(to_scan, max, pos));

    for (std::size_t i = 0; i < num_elem; ++i)
        if (!silent)
            std::cout << " " << max[i] << "(" << pos[i] << ")"
                      << "\n";

    return 0;
}