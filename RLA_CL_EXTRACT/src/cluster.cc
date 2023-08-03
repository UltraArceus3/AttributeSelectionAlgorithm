#include "cluster.h"

cluster::cluster(){
    this->id=0;
    this->records=records_t();
}

cluster::cluster(cluster_id Id, records_t Records){
    this->id=Id;
    this->records=records_t(Records);
}

count_t cluster::get_number_of_records(){
    return(this->records.size());
}

void cluster::clear(){
    this->id=0;
    this->records.clear();
}