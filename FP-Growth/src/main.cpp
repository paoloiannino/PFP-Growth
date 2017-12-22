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
#include <thread>

#include "fp_tree.h"

//#define PRINT_RESULT
#define ERR(s){cerr << "ERROR " << s << endl; exit(EXIT_FAILURE);}
#define ECHO(a){cout << a << endl;}

using namespace std;
using namespace std::chrono;

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b);
void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold);

vector<pair<vector<string>, int>> results;

int main(int argc, char **argv)
{
    int threshold (0);
    string line;
    string item;
    ifstream dataset;
    istringstream iss;
    vector<string> tmp;
    vector<string> pattern;
    vector<vector<string>> transactions;
    steady_clock::time_point begin_all;
    steady_clock::time_point begin_fp_only;
    steady_clock::time_point end;
    duration<double> fp_duration;

    if(argc != 3 || (threshold = stoi(argv[2])) < 0)
        ERR("Usage: ./fp-growth file threshold");

    begin_all = steady_clock::now();
    dataset.open(argv[1]);
    if(!dataset.is_open())
        ERR("Opening file@main");

    while(getline(dataset, line)){
        iss.str(line);
        tmp.clear();
        while(iss >> item)
            tmp.push_back(item);

        transactions.push_back(tmp);
        iss.clear();
    }

    pattern.push_back(string(""));
    begin_fp_only = steady_clock::now();
    fp_growth(transactions, pattern, threshold);
    end = steady_clock::now();

#ifdef PRINT_RESULT
    for(pair<vector<string>, int> result : results){
        for(string item : result.first)
            cout << item << " ";

        cout << "#" << result.second << endl;
    }
#endif

    cout.unsetf (ios::floatfield);
    cout.precision(5);
    fp_duration = duration_cast<duration<double>>(end - begin_all);
    cout << "TOTAL TIME: " << fp_duration.count() << endl;
    fp_duration = duration_cast<duration<double>>(end - begin_fp_only);
    cout << "FP TIME: " << fp_duration.count() << endl;

    return 0;
}

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b) {
    return a.second < b.second;
}

struct support_compare{
    support_compare(const map<string, int> &support_reference)
	: support_reference(support_reference) {}

    bool operator() (const string &a, const string &b){
	return support_reference.at(a) > support_reference.at(b);
    }

    map<string, int> support_reference;
};

void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold){
    map<string, int> support_count;
    vector<pair<string, int>> header_table;
    shared_ptr<vector<string>> tmp_vector;
    shared_ptr<fp_tree> ft = make_shared<fp_tree>();
    shared_ptr<vector<vector<string>>> tmp_transactions;

    for(vector<string> transaction : transactions)
        for(string item : transaction)
            if(support_count.find(item) == support_count.end())
            support_count[item] = 1;
            else support_count[item]++;

    for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++){
        header_table.push_back({it->first, it->second});
        tmp_vector = make_shared<vector<string>>(pattern);
        tmp_vector->push_back(it->first);

       if(it->second > threshold)
           results.push_back({*tmp_vector, it->second});
    }

    sort(header_table.begin(), header_table.end(), pair_compare);

    for(vector<string> transaction : transactions){
        sort(transaction.begin(), transaction.end(), support_compare(support_count));
        ft->insert_transaction(transaction.begin(), transaction.end());
    }

    for(pair<string, int> header_row : header_table){
        if(header_row.second > threshold){
            tmp_vector = make_shared<vector<string>>(pattern);
            tmp_vector->push_back(header_row.first);
            tmp_transactions = ft->get_transaction(header_row.first);
            fp_growth(*tmp_transactions, *tmp_vector, threshold);
        }
    }
}
