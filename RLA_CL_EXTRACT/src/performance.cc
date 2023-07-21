//============================================================================
// Name        : performance.cpp
// Author      : Ahmed Soliman
// Version     : 1.0
// Copyright   :
// Description : functions for computing linkage performance metrics
//============================================================================

#include "performance.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/range/adaptor/reversed.hpp>

using namespace std;
using namespace boost;
namespace logging = boost::log;
namespace keywords = boost::log::keywords;
namespace src = boost::log::sources;

//initialize logging
void performance_logging_init(string log_file_name="detailed_performance.log")
{
    logging::add_file_log
    (
        keywords::file_name = log_file_name,
        keywords::format = "[%TimeStamp%]: %Message%"
    );
    logging::core::get()->set_filter
        (
            logging::trivial::severity >= logging::trivial::info
        );
}

//ufi is the unique field index
//ufi is used to identify truly matching records, hence computing 
//performance metrics such as precision and sensitivity/recall
performance evaluateLinkagePerformance(clusters cs, uint8_t ufi, string log_file_name="detailed_performance.log") {

    performance_logging_init(log_file_name);
    logging::add_common_attributes();

    using namespace logging::trivial;
    src::severity_logger< severity_level > lg;

    BOOST_LOG_SEV(lg, info) << "Starting linkage performance evaluation";

    performance p;

    p.accuracy=0.0;
    p.precision=0.0;
    p.recall=0.0;
    p.f1_score=0.0;

    number_t nc=0; // number of clusters in the dataset
    number_t nr=0; // number of records in the dataset

    number_t i,j;
    mapUIDtoCount_t muidtc; //map uid to count, where count is the number of records having uid in the current cluster.
    mapUIDtoCIDs_t muidtcids; //map uid to cluster ids

    unique_id_t uidcr; //unique ID for current record

    number_t nrcc=0; //number of records in current cluster
    number_t ntmcc=0; //number of total matches in current cluster
                      //(this reflects the number of links whether true or false)

    //intermediate statistics to evaluate the linkage performance
    number_t dcc=0; //number of true matches in current cluster
    number_t d=0; //total number of true matches in the dataset
    number_t bcc=0; //number of false matches in current cluster
    number_t b=0; //total number of false matches in the dataset
    number_t c=0; //number of false non-matches
    number_t ntl=0; //number of total links
    number_t ndl=0; //number of detected links
    number_t ntf=0; //number of times found

    mapUIDtoCount_t::iterator it1;
    mapUIDtoCIDs_t::iterator it2;

    number_t ncrcc=0; //number of correct records in current cluster
    number_t ncr=0; //total number of correct records in the dataset


    set<unique_id_t> suidowners; //set of unique IDs assigned as cluster owners
    set<pair<clusterID_t,count_t>> sPairCIDCnt;
    // map cluster ID to uid of the owner of the cluster
    // and count of records belonging to the owner
    map<clusterID_t, pair<unique_id_t, count_t>> mcto;
    map<unique_id_t, count_t> nroc;  // number of records in the owned cluster
    map<unique_id_t, clusterID_t> uidtcid;  // map unique id to its owned cluster ID
    cluster cc;

    nc=cs.size(); //get number of clusters
    
    //get total number of records
    nr=0;
    for(i=0;i<nc;++i){
        nr+=cs.at(i).get_number_of_records();            
    }
    BOOST_LOG_SEV(lg, info) << "Total number of records = " << nr;

    //calculate b,d hence precision
    //for each cluster in the cluster set
    for(i=0;i<nc;i++) {

        cc=cs.at(i);
        nrcc=cc.get_number_of_records();
        muidtc.clear();

        //for each record in the current cluster
        for(j=0;j<nrcc;j++) {
            uidcr=cc.records[j][ufi]; //get unique ID value
            muidtc[uidcr]++; //update count
        }
        //sort uids by frequency in non-increasing order
        multimap<count_t, unique_id_t> cnttuid;
        for(auto& it : muidtc) {
            cnttuid.insert({it.second, it.first});
        }
        muidtc.clear();
        for(auto& it : boost::adaptors::reverse(cnttuid)) {
            count_t count = it.first;
            unique_id_t uid = it.second;
            muidtc[uid]=count;
        }
        //decide the owner of the current cluster
        for(auto& it : muidtc) {
            unique_id_t uid = it.first;
            count_t cnt = it.second;
            if (suidowners.find(uid)!=suidowners.end()) {
                // current uid has already been labeled as an owner for another cluster C.
                if (cnt <= nroc[uid]) {
                    continue; // skip if cluster C has more than or equal uid records than in current cluster.
                }
                else { // disassociate uid from the cluster C and associate it with current cluster.
                    BOOST_LOG_SEV(lg, info) << "Change in cluster ownership: "
                                << "Owner ID: " << uid
                                << "  Old cluster: " << uidtcid[uid]
                                << "  New cluster: " << cc.id;
                    mcto.erase(uidtcid[uid]);
                    mcto[cc.id] = pair<unique_id_t, count_t>(uid, cnt);
                    uidtcid[uid] = cc.id;
                    nroc[uid] = cnt;
                }
            }
            // Assign cluster owner
            suidowners.insert(uid);   //add unique id for current cluster to the set of uids assigned as cluster owners
            uidtcid[uid] = cc.id;
            mcto[cc.id] = pair<unique_id_t, count_t>(uid, cnt);
            nroc[uid] = cnt;
            // log cluster owner
            BOOST_LOG_SEV(lg, info) << "Cluster " << cc.id << ": "
                                    << "Owner ID: " << uid;
        }

        dcc=0;
        for(it1=muidtc.begin();it1!=muidtc.end();++it1) {
            count_t cnt = it1->second;
            dcc+=0.5*(cnt*(cnt-1));
            muidtcids[it1->first].insert(pair<clusterID_t,count_t>(cc.id,cnt)); //i is the cluster ID or index
        }
        ntmcc=0.5*(nrcc*(nrcc-1));
        bcc=ntmcc-dcc;

        // log erroneous clusters
        if(bcc>0) {
            BOOST_LOG_SEV(lg, info) << "Error at cluster " << cc.id << ": ";
            for(it1=muidtc.begin();it1!=muidtc.end();++it1) {
                BOOST_LOG_SEV(lg, info) << it1->first << " => " << it1->second;
            }
        }

        //update running totals
        d+=dcc;
        b+=bcc;
    }

    // sum number of correct records based on cluster ownership
    for(auto& it: nroc){
        ncr += it.second;
    }

    //now scan the muidtcids map to compute c (number of false non-matches)
    //this statistic is required for computing sensitivity (recall)

    for(it2=muidtcids.begin();it2!=muidtcids.end();++it2) {
        sPairCIDCnt=it2->second;
        size_t s=sPairCIDCnt.size();
        if(s>1)
        {
            // log false non-matches count
            BOOST_LOG_SEV(lg, info) << "False non-matches detected for records with unique ID " << it2->first << ": ";
            BOOST_LOG_SEV(lg, info) << "  above records are scattered across " << s << " clusters as follows:";
            ntl=0;
            ndl=0;
            ntf=0;
            for(set<pair<clusterID_t,count_t>>::iterator it3=sPairCIDCnt.begin();it3!=sPairCIDCnt.end();++it3){
                // log false non-matches details
                BOOST_LOG_SEV(lg, info) << "   " << it3->second << " records in cluster " << it3->first;
                count_t cnt=it3->second;
                ndl+=0.5*(cnt*(cnt-1));
                ntf+=cnt;
            }
            ntl=0.5*(ntf*(ntf-1));
            // log false non-matches link counts
            BOOST_LOG_SEV(lg, info) << "  Total number of records = " << ntf;
            BOOST_LOG_SEV(lg, info) << "  Total number of links = " << ntl;
            BOOST_LOG_SEV(lg, info) << "  Total number of detected links = " << ndl;
            BOOST_LOG_SEV(lg, info) << "  Total number of missed links = " << ntl-ndl;
            c+=(ntl-ndl);
        }
    }

    p.nr=nr;
    p.nc=nc;
    p.ncr=ncr;
    p.b=b;
    p.c=c;
    p.d=d;
    count_t links_total=0.5*(nr*(nr-1));
    p.a=links_total-b-c-d;
    p.accuracy=(p.a+p.d)/((double)(p.a+p.d+p.b+p.c));

    p.precision=d/((double)(b+d));
    p.recall=d/((double)(c+d));
    p.accuracy2=ncr/(double)nr;
    p.f1_score=2*p.precision*p.recall/(p.precision+p.recall);
    
    return p;
}

void printPerformance(performance p,ostream& os){
        
    cout.imbue(std::locale("")); //for thousands comma separators

    string emphs="",emphe="";
    
    if(&os == &cout){
        emphs="\033[1m";
        emphe="\033[0m";
    }

    os << endl
       << emphs
       << "Statistics" << emphe << endl;
    
    os << "\t" << "nc  = " << p.nc << endl;
    os << "\t" << "nr  = " << p.nr << endl;
    os << "\t" << "ncr = " << p.ncr << endl;

    os << "\t" << "b   = " << p.b << endl;
    os << "\t" << "c   = " << p.c << endl;
    os << "\t" << "d   = " << p.d << endl;

    os << endl
       << emphs
       << "Linkage performance metrics" << emphe << endl;
    
    os << setprecision(4);
    os << "\t" << "Accuracy  = " << p.accuracy   << " = " << p.accuracy*100.0 << "%" << endl;
    os << "\t" << "Precision = " << p.precision  << " = " << p.precision*100.0 << "%" << endl;
    os << "\t" << "Recall    = " << p.recall     << " = " << p.recall*100.0 << "%" << endl;
    os << "\t" << "F1 Score  = " << p.f1_score   << " = " << p.f1_score*100.0 << "%" << endl;

    os << endl;

}

// load linkage output clustering data from file 
// for performance evaluations
clusters loadLinkageOutput(string fileLinkageOutput) {

    record_t r;  //record
    cluster c;   //cluster
    clusters cs; //clusters

    ifstream ifs;
    string lineStr;

    regex re("\\s+"); // \s matches a space, a tab, a line break, or a form feed.
    regex re_cluster("cluster\\s+(\\d+):\\s+");
    smatch what;

    ifs.open(fileLinkageOutput, ifstream::in);
    
    r.clear(); c.clear(); cs.clear();

    if(ifs.fail())
        cerr<<"Error: can not open the file '"<<fileLinkageOutput<<"'"<<endl;
    
    if(ifs.good())
        getline(ifs, lineStr);

    while(ifs.good())
    {
        if(regex_match(lineStr,re_cluster)){ 

            //read cluster id
            c.clear();
            if(regex_search(lineStr,what,re_cluster)){
                c.id=stol(what[1].str());
            }
            
            //fetch cluster records
            while(ifs.good())
            {
                getline(ifs, lineStr);
                if(lineStr=="" || regex_match(lineStr,re_cluster)){
                    cs.push_back(c); //end of cluster detected thus push the temp cluster into the set of clusters
                    break;
                }
                r.clear(); //clear the temporary record
                //populate the temporary record with current data from input file
                for(sregex_token_iterator it = sregex_token_iterator(lineStr.begin(), lineStr.end(), re, -1); it != sregex_token_iterator(); ++it){    
                    string temp=(string)*it;
                    boost::algorithm::to_lower(temp);
                    r.push_back(temp);
                }
                //push the new temp record into the temp cluster
                c.records.push_back(r);
            } 
        }
        else
        {
            if(ifs.good())
                getline(ifs, lineStr);
        }
        
    }
    ifs.close();
    return cs; //return the clusters loaded from the input file
}