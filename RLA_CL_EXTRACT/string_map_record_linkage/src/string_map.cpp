//
// Created by Ahmed Soliman on 2/27/22.
//

#include "string_map.h"
#include <cassert>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <boost/log/trivial.hpp>
#include <utility>
#include <unordered_map>
#include <map>
#include <unordered_set>


/* constructor */
string_map::string_map(std::vector<std::string> const& strings,
        DistFuncCallback dist_func,
        DIM n_dimensions) : strings(strings), n_strings_input(strings.size()),
                            n_dimensions(n_dimensions), distance_function(std::move(dist_func))
                            {
    /* Deduplicate the strings */
    std::unordered_set<std::string> s(this->strings.begin(), this->strings.end());
    this->strings.assign(s.begin(), s.end());
    this->n_strings = this->strings.size();
    /* Sort the strings */
    std::sort(this->strings.begin(), this->strings.end());
    this->pivots.reserve(n_dimensions);
    Point p;
    p.resize(this->n_dimensions);
    std::fill(p.begin(), p.end(), 0);
    this->coord.resize(this->n_strings);
    std::fill(coord.begin(), coord.end(), p);
    std::random_device rd;
    this->generator = std::mt19937(rd());
    this->distribution = std::uniform_int_distribution<COUNT>(0, this->n_strings - 1);
    BOOST_LOG_TRIVIAL(info) << "Number of input strings = " << this->n_strings_input;
    BOOST_LOG_TRIVIAL(info) << "Number of unique strings = " << this->n_strings;
    BOOST_LOG_TRIVIAL(info) << "Number of dimensions = " << this->get_n_dimensions();

}

/* load from file constructor */
/* load strings and corresponding coordinates from a file
 * also load all member variables so that the string mapping could be reused
 * incrementally and efficiently after offline computation.*/
string_map::string_map(const std::string& file_name) {
    std::fstream fs;
    std::string buffer(256,'\0');
    fs.open(file_name.c_str(), std::ios::in);
    if ( (fs.rdstate() & std::ifstream::failbit ) != 0 ) {
        BOOST_LOG_TRIVIAL(error) << "Can not open file '" << file_name.c_str() << "'" << std::endl;
        exit(1);
    }
    fs >> this->n_strings_input;
    fs >> this->n_strings;
    int x;
    fs >> x;
    this->n_dimensions = (DIM)x;
    double y;
    fs >> y;
    this->mapping_time = (std::chrono::duration<double, std::milli>)y;
    this->pivots.reserve(this->n_dimensions);
    std::string p1;
    std::string p2;
    COUNT i;
    for(i=0; i<this->n_dimensions; ++i) {
        fs >> p1;
        fs.get();
        fs >> p2;
        this->pivots.push_back(std::vector<std::string>({p1,p2}));
    }
    this->strings.reserve(this->n_strings);
    Point p;
    p.resize(this->n_dimensions);
    std::fill(p.begin(), p.end(), 0);
    this->coord.resize(this->n_strings);
    std::fill(coord.begin(), coord.end(), p);
    i=0;
    for(; i<this->n_strings; ++i) {
        std::string str;
        fs >> str;
        fs.get();
        this->strings.push_back(str);
        int j=0;
        double d=0.0;
        for(; j < this->n_dimensions-1; ++j) {
            fs >> std::skipws >> d;
            fs.get();
            this->coord[i][j] = d;
        }
        fs >> d;
        this->coord[i][j] = d;
    }
    fs.close();
    this->populate_str2point();
    std::random_device rd;
    this->generator = std::mt19937(rd());
    this->distribution = std::uniform_int_distribution<COUNT>(0, this->n_strings - 1);
    BOOST_LOG_TRIVIAL(info) << "Loading mappings from file " << file_name;
    BOOST_LOG_TRIVIAL(info) << "Number of input strings = " << this->n_strings_input;
    BOOST_LOG_TRIVIAL(info) << "Number of unique strings = " << this->n_strings;
    BOOST_LOG_TRIVIAL(info) << "Number of dimensions = " << this->get_n_dimensions();
    BOOST_LOG_TRIVIAL(info) << "Mapping took " << this->get_mapping_time() << " ms";
    BOOST_LOG_TRIVIAL(info) << "Loading done";
}


int string_map::original_metric(INDEX ind1, INDEX ind2) {
    int distance;
    distance = this->distance_function(this->strings.at(ind1), this->strings.at(ind2));
    return distance;
}

int string_map::get_n_dimensions() const{
    return this->n_dimensions;
}

COUNT string_map::get_n_strings() const{
    return this->n_strings;
}

/* choose two distant strings as pivots.
   these pivots define the h-th dimension in the target space.
   <h> is the current dimension.
   returns <a,b> the indices of the two pivot strings */
std::vector<COUNT> string_map::choose_pivots(int h) {
    assert( h < this->n_dimensions );
    /* get a random string index from the input strings */
    COUNT a = this->distribution(this->generator);
    COUNT b;
    for(int i=0;i<5;++i) {
        /* get <b>, the index of the farthest point from sa
         * where sa is the string indexed by <a> */
        Distance max = 0; // maximum distance encountered so far
        COUNT max_b = 0;  // index of the maximum distance
        for(b=0; b < this->n_strings; ++b) {
            Distance dist = this->get_distance(a, b, h);
            if(dist > max) {
                max = dist;
                max_b = b;
            }
        }
        b = max_b;
        /* get <a>, the index of the farthest point from sa
         * where sa is the string indexed by a */
        max = 0; // maximum distance encountered so far
        COUNT max_a = 0;  // index of the maximum distance
        for(a=0; a < this->n_strings; ++a) {
            Distance dist = this->get_distance(a, b, h);
            if(dist > max) {
                max = dist;
                max_a = a;
            }
        }
        a = max_a;
    }
    return std::vector<COUNT>({a,b});
}


Distance string_map::get_distance(COUNT a, COUNT b, int h) {
    assert( h < this->n_dimensions);
    /* get the two strings */
    std::string sa = this->strings[a];
    std::string sb = this->strings[b];
    /* get original metric distance */
    Distance distance = this->original_metric(a,b);
    for(int i=0;i<(h-1);++i) {
        /* get the difference in dimension i */
        Distance difference = this->coord[a][i] - this->coord[b][i];
        distance = sqrt(fabs(distance*distance - difference*difference));
    }
    return distance;
}


Distance string_map::get_euclidean_distance(COUNT a, COUNT b) const {
    assert(a < this->coord.size());
    assert(b < this->coord.size());
    Distance distance = distances::euclid_dist(coord[a], coord[b]);
    return distance;
}

Distance string_map::get_q_gram_distance(std::string&  a, std::string&  b) const {

    Distance distance = distances::q_gram_distance(a,b);
    return distance;
}


Distance string_map::get_squared_euclidean_distance(COUNT a, COUNT b) const {
    assert(a < this->coord.size());
    assert(b < this->coord.size());
    Distance distance = distances::squared_euclid_dist(coord[a], coord[b]);
    return distance;
}


bool string_map::is_mapped() const{
    return this->mapped;
}


/* Map the input strings (objects) into images in the target Euclidean space. For each axis, choose
   a pair of distant pivot strings, compute the distance between the pivot strings, if the distance is zero,
   set all coordinates in this axis to zeros, if not zero, compute coordinates of strings on current axis
   based on the distance between strings and each of the two pivots as well as the distance between the pivots.*/
void string_map::map() {
    BOOST_LOG_TRIVIAL(info) << "Mapping strings started";
    auto t1 = std::chrono::high_resolution_clock::now();
    for(int h=0; h < this->n_dimensions; ++h) {
        std::vector<COUNT> P = this->choose_pivots(h);
        /* store pivot strings */
        this->pivots.push_back({this->strings[P[0]], this->strings[P[1]]});
        Distance distance = this->get_distance(P[0], P[1], h);
        if(distance==0) {
            for(COUNT i=0;i < this->n_strings; ++i) {
                this->coord[i][h] = 0;
            }
            break;
        }
        /* compute coordinates of strings on current axis */
        for(COUNT i=0; i < this->n_strings; ++i) {
            Distance x = this->get_distance(i, P[0], h);
            Distance y = this->get_distance(i, P[1], h);
            this->coord[i][h] = (x*x + distance*distance - y*y) / (2*distance);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    this->mapping_time = t2 - t1;
    this->populate_str2point();
    this->mapped = true;
    BOOST_LOG_TRIVIAL(info) << "Mapping strings done in " << this->get_mapping_time() << " ms";
}


double string_map::get_mapping_time() const {
    return this->mapping_time.count();
}


/* save strings and corresponding coordinates into a file
 * also save all member variables so that the string mapping could be done
 * incrementally and offline.*/
void string_map::save_mappings(const std::string& file_name) {
    BOOST_LOG_TRIVIAL(info) << "Saving mappings to file " << file_name;
    std::fstream fs;
    fs.open(file_name.c_str(), std::ios::out);
    fs << this->n_strings_input << std::endl;
    fs << this->n_strings << std::endl;
    fs << (int)this->n_dimensions << std::endl;
    fs << this->mapping_time.count() << std::endl;
    COUNT i;
    for(i=0; i<this->n_dimensions; ++i){
        fs << this->pivots[i][0] << " "
           << this->pivots[i][1] << std::endl;
    }
    i=0;
    for(auto it=this->strings.begin(); it!=strings.end();++it) {
        // output the string
        fs << std::setw(15) << std::left << *it << " ";
        // output the coordinates
        int j=0;
        for(; j < this->n_dimensions-1; ++j) {
            fs << std::setw(6) << std::right << std::setprecision(5) << coord[i][j] << " ";
        }
        fs << std::setw(6) << std::right << std::setprecision(5) << coord[i][j] << std::endl;
        ++i;
    }
    fs.close();
}

/* for all pairs of strings (or up to at least a given number of record pairs)
 *   compute the Euclidean distance from the coordinates
 *   compute the Original edit distance
 *   insert the Euclidean distance into a map of Edit distance unique values
 *   save summary statistics to file */
void string_map::save_summary(const std::string &file_name, size_t n_record_pairs, bool squared_dist) {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::fstream fs;
    fs.open(file_name.c_str(), std::ios::out);
    std::unordered_map<int,std::vector<Distance>> M;  //Edit distance to Euclidean distances map
    COUNT rp = 0;  // record pair counter
    while(rp < n_record_pairs) {
        INDEX i = this->distribution(this->generator);
        INDEX j = this->distribution(this->generator);
        if(i==j) continue;
        std::string str1 = this->strings[i];
        std::string str2 = this->strings[j];
        if (str1.size() > 1 && str2.size() > 1){
            Distance euclidean_distance = squared_dist ?
            this->get_squared_euclidean_distance(i,j) :
            this->get_euclidean_distance(i, j);
            int edit_distance = original_edit_dist(str1, str2);
            Distance q_gram_distance = this-> get_q_gram_distance(str1,str2);
        
            auto it = M.find(edit_distance);
            if(it==M.end()) {
                M.insert({edit_distance, {q_gram_distance}});
            } else {
                it->second.push_back(q_gram_distance);
        }
        ++rp;
        }
       
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> t_summary = t2-t1;
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::endl;
    ss << "Number of input strings             = " << this->n_strings_input << std::endl;
    ss << "Number of unique strings            = " << this->get_n_strings() << std::endl;
    ss << "Number of dimensions                = " << this->get_n_dimensions() << std::endl;
    ss << "Mapping time                        = " << this->mapping_time.count() << " ms" << std::endl;
    ss << "Summary time                        = " << t_summary.count() << " ms" << std::endl;
    ss << "Summary time per 1M record pairs    = " << t_summary.count()*1e6/(double)rp << " ms" << std::endl;
    ss << "Number of record pairs investigated = " << n_record_pairs << std::endl << std::endl;
    ss << "================================================" << std::endl;
    ss << "Edit  Euclidean Euclidean Euclidean Number of   " << std::endl;
    ss << "Dist  Min       Avg       Max       Record paris" << std::endl;
    ss << "================================================" << std::endl;
    std::map<int,std::vector<Distance>> M2(M.begin(), M.end());
    for( const auto& it : M2) {
        ss << std::right << std::setw(4) << it.first << "  ";
        Distance min = *std::min_element(it.second.begin(), it.second.end());
        Distance max = *std::max_element(it.second.begin(), it.second.end());
        Distance sum = 0;
        for (auto& it2 : it.second) {
            sum += it2;
        }
        COUNT size = it.second.size();
        Distance avg = sum / (double) size;
        ss << std::setprecision(4) << std::setw(10) << std:: left << min;
        ss << std::setprecision(4) << std::setw(10) << std:: left << avg;
        ss << std::setprecision(4) << std::setw(10) << std:: left << max;
        ss << std::setprecision(4) << std::setw(10) << std:: left << size;
        ss << std::endl;
    }
    ss << "================================================" << std::endl;
    fs << ss.str();
    BOOST_LOG_TRIVIAL(info) << ss.str();
    Distance tau = 6.0;  // Euclidean distance threshold
    Distance theta = 1; // Edit distance threshold
    while(true) {
        BOOST_LOG_TRIVIAL(info) << "Enter the value of theta: ";
        std::cin >> theta;
        BOOST_LOG_TRIVIAL(info) << "Enter the value of tau (negative number to exit): ";
        std::cin >> tau;
        if (tau<0) break;
        BOOST_LOG_TRIVIAL(info) << "tau = " << tau;
        /* compute the pruning ratio */
        COUNT n_skips = 0; // How many times Euclidean distance > tau
        COUNT tp = 0; // How many times (Euclidean distance > tau) and (Edit distance > theta) (success)
        COUNT fp = 0; // How many times (Euclidean distance > tau) and (Edit distance <= theta) (failure)
        COUNT total = 0; // Total number of record pairs
        COUNT lttheta = 0; // Total number of record pairs where Edit distance <= theta
        for (const auto &it: M2) {
            auto edit_dist = (Distance) it.first;
            for (auto &it2: it.second) {
                auto euclid_dist = it2;
                total++;
                if (edit_dist <= theta) lttheta++;
                if (euclid_dist > tau) {
                    n_skips++;  // here we assume edit distance > theta
                    (edit_dist > theta) ? tp++ /* hit */: fp++ /*miss*/;
                }
            }
        }
        ss.str("");
        ss << std::endl;
        ss << "================================================" << std::endl;
        ss << "theta   = " << theta << std::endl;
        ss << "tau     = " << tau << std::endl;
        ss << "n_skips = " << n_skips << " / " << total << " ( " << (double) n_skips * 100.0 / (double) total << " % ) "
           << std::endl;
        ss << "tp      = " << tp << " / " << n_skips << " ( " << (double) tp * 100.0 / (double) n_skips << " % ) "
           << std::endl;
        ss << "fp      = " << fp << " / " << lttheta << " ( " << (double) fp * 100.0 / (double) lttheta << " % ) "
           << std::endl;
        ss << "================================================" << std::endl;
        fs << ss.str();
        BOOST_LOG_TRIVIAL(info) << ss.str();
    }
    fs.close();
}

/* empty constructor */
string_map::string_map() = default;

/* get coordinates for input string */
Point string_map::get_coordinates(const std::string& str) {
    auto it = this->str2point.find(str);
    Point p;
    if(it==str2point.end()) {
        p.resize(this->n_dimensions);
        std::fill(p.begin(), p.end(), 0);
    } else {
        p = it->second;
    }
    return p;
}

void string_map::populate_str2point() {
    COUNT i=0;
    for(auto it=this->strings.begin(); it!=strings.end();++it) {
        this->str2point.insert({*it,this->coord[i++]});
    }
}
