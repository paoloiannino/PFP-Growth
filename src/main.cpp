#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <map>
#include <vector>
#include <set>
#include <climits>      // For INT_MAX

#include "fp_tree.h"

using namespace std;

vector<pair<vector<string>, int>> result;

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b) {
    return a.second < b.second;
}

bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b) {
    return a.second > b.second;
}

struct support_compare{
    support_compare(const map<string, int> &support_reference)
        : support_reference(support_reference) {};

    bool operator() (const string &a, const string &b){
        return support_reference.at(a) > support_reference.at(b);
    }

    map<string, int> support_reference;
};


void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold){
    map<string, int> support_count;
    vector<pair<string, int>> header_table;
    vector<string> *tmp_vector;
    fp_tree *ft = new fp_tree();
    shared_ptr<vector<vector<string>>> tmp_transactions;
    vector<pair<vector<string>, int>> *tmp_result = &result;


    for(vector<string> transaction : transactions)
        for(string item : transaction)
            if(support_count.find(item) == support_count.end())
                support_count[item] = 1;
            else support_count[item]++;

    for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++){
        header_table.push_back({it->first, it->second});
        tmp_vector = new vector<string>(pattern);
        tmp_vector->push_back(it->first);

        if(it->second > threshold)
            result.push_back({*tmp_vector, it->second});
    }

    sort(header_table.begin(), header_table.end(), pair_compare);

    for(vector<string> transaction : transactions){
        sort(transaction.begin(), transaction.end(), support_compare(support_count));
        ft->insert_transaction(transaction.begin(), transaction.end());
    }

    support_count.clear();
    for(pair<string, int> header_row : header_table){
        tmp_vector = new vector<string>(pattern);
        tmp_vector->push_back(header_row.first);
        tmp_transactions = ft->get_transaction(header_row.first);
        fp_growth(*tmp_transactions, *tmp_vector, threshold);
    }

    delete ft;

    return;
}

int main(int argc, char **argv)
{
    int threshold (0);
    int word_per_line (INT_MAX);
    int shortest_itemset (INT_MAX);
    int i (0);
    ifstream dataset;
    istringstream iss;
    string item;
    string line;
    vector<vector<string>> transactions;
    vector<string> *tmp;

    if(argc == 5
            && (threshold = stoi(argv[2])) >=0
            && (word_per_line = stoi(argv[3])) > 0
            && (shortest_itemset = stoi(argv[4])) > 0){
        dataset.open(argv[1]);

        if(!dataset.is_open())
            cerr << "ERROR: opening file" << endl;

    } else {
        cout << "Usage: ./fpgrowth file threshold wordPerLine shortestItemset" << endl;
    }

    while(getline(dataset, line)){
        iss.str(line);
        tmp = new vector<string>;
        i = 0;
        while(iss >> item){
            tmp->push_back(item);
            i++;

            if(i == word_per_line) break;
        }

        transactions.push_back(*tmp);
        iss.clear();
    }

    vector<string> pattern;
    pattern.push_back(string());

    fp_growth(transactions, pattern, threshold);
    sort(result.begin(), result.end(), pair_compare_result);

    for(pair<vector<string>, int> pattern : result){
        if(pattern.second < threshold
                || pattern.first.size() < (shortest_itemset + 1)) continue;

        sort(pattern.first.begin(), pattern.first.end());
        for(string item : pattern.first)
            cout << item << " ";

        cout << "SUPPORT: " << pattern.second << endl;
    }

    return 0;
}
