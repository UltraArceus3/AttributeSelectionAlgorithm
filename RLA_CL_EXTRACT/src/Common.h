/*
 * Common.h
 *
 *  Created on: Jul 30, 2014
 *      Author: abdullah-al-mamun
 */

#ifndef COMMON_H_
#define COMMON_H_


typedef struct
{
	std::string str;
	int ind;
}StrPacket;

typedef struct
{
    int n_datasets = 0;
    std::string datasets_file_names = std::string(256, '\0');
    long n_records_total_read = 0;
    int l_mer = 0;
    int threshold = 0;
    double linkage_time_in_seconds = -1;
    double accuracy_type[4] = {0.0};
	long n_records_type[4] = {0};
    double accuracy = -0.1;
    double precision = -0.1;
    double recall = -0.1;
}LinkParam;

#endif /* COMMON_H_ */
