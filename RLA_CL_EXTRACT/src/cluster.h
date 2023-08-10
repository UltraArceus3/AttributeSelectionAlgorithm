#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

#include <vector>
#include <string>

typedef unsigned int cluster_id;

typedef std::vector<std::string> record_t;
typedef std::vector<record_t> records_t;
typedef unsigned long long int count_t;

class cluster{

    public:
        cluster_id id; //cluster identifier
        records_t records; //records in this cluster
        cluster();
        cluster(cluster_id Id, records_t Records);
        count_t get_number_of_records();
        void clear();
};

typedef std::vector<cluster> clusters;

#endif