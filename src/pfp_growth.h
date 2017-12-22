#ifndef PFP_GROWTH_H
#define PFP_GROWTH_H

#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <map>
#include <vector>
#include <set>
#include <climits>      // For INT_MAX
#include <omp.h>
#include <chrono>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <thread>

#include "fp_tree.h"

#define N_THREADS       2
#define OMP_FOR_SCHEDULE static
#define BUFFER_SIZE     8192
#define MASTER          0
#define NO_TAG          0
#define END             "***END***"
#define ESCAPE_CHAR     string("#")

#define ERR(s){cerr << "ERROR " << s << endl; exit(EXIT_FAILURE);}
#define ECHO(a){cout << a << endl;}

using namespace std;
using namespace std::chrono;

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b);
bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b);
void share_dataset(char *dataset_path);
void split_dataset(char *dataset_path);
void read_dataset(char *filename, vector<vector<string>> &transactions, map<string, int> &support_count);
void read_transactions(vector<vector<string>> &transactions, map<string, int> &support_count);
void read_chunk(vector<vector<string>> &transactions, map<string, int> &support_count);
string parse_result(const vector<string> &itemset, const int &support);
void send_result(vector<string> itemset, const int &support, const int &threshold);
void send_end();
void get_results(const int &mpi_num_processes, map<string, int> &results);
void pfp_growth(const vector<vector<string>> &transactions, const map<string, int> &support_count, int threshold);
void pfp_growthR(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold);
void pfp_growth_local(const vector<vector<string>> &transactions, const map<string, int> &support_count, int threshold, map<string, int> &results);
void pfp_growthR_local(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold, map<string, int> &results);



#endif // PFP_GROWTH_H
