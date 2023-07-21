#ifndef PERFORMANCE_H_INCLUDED
#define PERFORMANCE_H_INCLUDED

#include "cluster.h"
#include <map>
#include <set>
#include <string>
#include <locale>

typedef unsigned long long int number_t;
typedef unsigned long long int count_t;
typedef long long int clusterID_t;

typedef std::string unique_id_t; //type for unique identifier used in performance evaluation

typedef std::map<unique_id_t,count_t>::value_type mapUIDtoCount_entry;
typedef std::map<unique_id_t,count_t> mapUIDtoCount_t;

typedef std::map<unique_id_t,std::set<std::pair<clusterID_t,count_t>>> mapUIDtoCIDs_t;
typedef mapUIDtoCIDs_t::value_type mapUIDtoCIDs_entry;

// structure that holds all statistics and performance metric values
typedef struct performance_struct{
    
    // Statistics about records in the dataset
    number_t nr=0;  // number of records in the dataset
    number_t nc=0;  // number of clusters in the dataset
    number_t ncr=0; // number of correct records in the dataset

    // Confusion matrix elements
    number_t a=0; //total number of true non-matches in the dataset
    number_t b=0; //total number of false matches in the dataset
    number_t c=0; //number of false non-matches
    number_t d=0; //total number of true matches in the dataset

    // Computed performance evaluation metrics
    double accuracy;
    double accuracy2;  //based on cluster ownership
    double precision;
    double recall;
    double f1_score;

} performance;

// initialize logging
void performance_logging_init(std::string log_file_name);

// evaluate the linkage performance for the dataset given
// ufi is the unique field identifier. It indicates the index
// of the field that is to be used as a gold standard
performance evaluateLinkagePerformance(clusters cs, uint8_t ufi, std::string log_file_name);

// print statistics and performance metrics
void printPerformance(performance p,std::ostream &os);

// load linkage output clustering data from file 
// for performance evaluations
clusters loadLinkageOutput(std::string fileLinkageOutput);

#endif
