#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <map>
#include <omp.h>


int numThreads = 6;

std::vector<std::string> generate_kmer(const std::string& str_d,int k){
    int length = str_d.size();

    std::vector<std::string> k_mer_list;

    if (k >= length){
        k_mer_list.push_back(str_d);
    }
    else {
        for(int i = 0; i < (length - (k-1));i++){
		    k_mer_list.push_back(str_d.substr(i,k));
	    }
    }
    return k_mer_list;
}

int calculateBasicED2(std::string& str1, std::string& str2, int threshRem,int threshold)
{
	int row, col, i, j;
	std::vector < std::vector < int > > matArr;

	row		= str1.length() + 1;
	col 	= str2.length() + 1;

	matArr.resize(row);
	for(i = 0; i < row; ++i)
		matArr[i].resize(col, 0);

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if((int)str1[i-1] == (int)str2[j-1])
					matArr[i][j]	= std::min(std::min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] 	= std::min(std::min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
			}

			if((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return threshold + 1;
			}
		}
	}

	return (matArr[row-1][col-1]);
}

// calculates edit distance between two string
// takes two strings and a threshold value as input
// returns global variable threshold + 1 if distance exceeds theshold param
// returns edit distance
// core mechanism is a DP algo 

int calculateBasicED(std::string& str1, std::string& str2, int threshRem,int threshold)
{
	int dist	= threshRem;

	if(abs((int)(str1.length() - str2.length())) > dist)
		return threshold + 1;
	else if(str1.compare(str2) == 0)
		return 0;
	else if(dist == 0)
		return threshold + 1;
	else if((2 * dist + 1) >= std::max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist,threshold);
	else
	{
		std::string s1, s2;
		int row, col, diagonal;
		int i, j;
		std::vector<std::vector<int> > matArr;

		if (str1.length() > str2.length())
		{
			s1 = str2;
			s2 = str1;
		}
		else
		{
			s1 = str1;
			s2 = str2;
		}

		row	 		= s1.length() + 1;
		col 		= 2 * dist + 1;
		diagonal 	= dist + s2.length() - s1.length();

		matArr.resize(row);
		for(i = 0; i < row; ++i)
			matArr[i].resize(col, 0);


		for(i = 0; i < dist + 1; i++)
		{
			for(j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j]	= j - dist;
				else if(j == (dist - i))
					matArr[i][j] 	= matArr[i - 1][j + 1] + 1;
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= std::min(std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= std::min(std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= std::min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if((j == diagonal) && matArr[i][j] > dist)
					return threshold + 1;
			}
		}

		for(i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for(j = 0; j < col; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= std::min(std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= std::min(std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= std::min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		for(i = s2.length() - dist + 1; i < row; i++)
		{
			for(j = 0; j < col - i + s2.length() - dist; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= std::min(std::min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= std::min(std::min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		return matArr[row - 1][diagonal];

	}
}





// Function to calculate Levenshtein distance
int LevenshteinDistance(const std::string& a, const std::string& b) {
    int m = a.size();
    int n = b.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            dp[i][j] = std::min({dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost});
        }
    }

    return dp[m][n];
}



std::vector<int> EdgeGeneration(std::vector<std::vector<std::string>>& deDuplicate,
                                const std::vector<int>& samp,
                                int blockingAttribute,
                                std::map<std::string,std::vector<int>> dictBlocks,
                                int threshold) {
    std::vector<int> recordPreOut;

    // std::cout << "We are here!" << std::endl;
    // for (int i=0; i < dictBlocks_key.size();i++) {
    //     //std::cout << "What's happening here! " << i << std::endl;
    //     for (int j = 0; j < dictBlocks_val[i].size();j++) {
            
    //         if (dictBlocks_val[i][j] > deDuplicate.size()) {
    //             std::cout << dictBlocks_val[i][j] << std::endl;
    //         }
    //     }

    //     if (dictBlocks_val[i].size() > 0) {
    //         dictBlocks[dictBlocks_key[i]] = dictBlocks_val[i];
    //     }
    //     else {
    //         dictBlocks[dictBlocks_key[i]] = std::vector<int>();
    //     }
        
    // }

    //std::cout << samp.size() <<" C++ Samp Size" << std::endl;

    //std::cout << "We are till here! Thanks god" << std::endl;
    #pragma omp parallel num_threads(numThreads) 
	{
        std::vector<int> tmp_out;
        #pragma omp for schedule(dynamic,1) nowait
        for (int i = 0; i < samp.size(); i++) {

            std::vector<std::string>& record = deDuplicate[samp[i]];
            std::vector<std::string> kmerList = generate_kmer(record[blockingAttribute], 3);
            std::set<int> kmerIndex;

            //std::cout << "Works till C++ call!" << std::endl;

            for (std::string& kmer : kmerList) {
                transform(kmer.begin(), kmer.end(), kmer.begin(), ::tolower);
                if (dictBlocks.find(kmer) != dictBlocks.end()) {
                    
                    auto it = dictBlocks.find(kmer);
                    //std::cout << it->first;
                    std::vector<int> tmp = it->second;
                    // if (kmer == "wli") {
                    //     std::cout << tmp[i] << std::endl;
                    // }
                    for (int j = 0; j < tmp.size();j++){
                        kmerIndex.insert(tmp[j]);
                    }
                }
            }
            //std::cout << "Works till First for loop" << std::endl;
            
            for (int index : kmerIndex) {
                

                std::vector<std::string> tmpRecord = deDuplicate[index];
                // std::cout << "ran dedup" << std::endl;
                // std::cout << index << std::endl;
                // std::cout << deDuplicate.size() << std::endl;
                
                int lDist = 0;

                for (int j = 1; j < tmpRecord.size(); j++) {
                    std::string& a = tmpRecord[j];
                    std::string& b = record[j];
					

                    //lDist += LevenshteinDistance(a, b);
                    lDist += calculateBasicED(a,b,threshold,threshold);
                    if (lDist > threshold) {
                        lDist = threshold + 1;
                        break;
                    }
                }

                if (lDist <= threshold) {
                    tmp_out.push_back(index);
                   // recordPreOut.push_back(index);
                }
            }
        }

        #pragma omp critical
		{
				//set<std::pair<int,int>>::iterator itr;
				//edgeList.insert(edgeList.end(),set_of_edges.begin(), set_of_edges.end());
                recordPreOut.insert(recordPreOut.end(),tmp_out.begin(),tmp_out.end());
                /*
					DEBUGGING PURPOSES
					cout << set_of_edges.size() << endl;
				*/
				

		}

    
    }

    return recordPreOut;
}

