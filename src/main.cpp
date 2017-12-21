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
#define BUFFER_SIZE     8192
#define MASTER          0
#define NO_TAG          0
#define END             "***END***"
#define ESCAPE_CHAR     string("#")
//#define DEBUG

#define ERR(s){cerr << "ERROR " << s << endl; exit(EXIT_FAILURE);}
#define ECHO(a){cout << a << endl;}

using namespace std;
using namespace std::chrono;

vector<pair<vector<string>, int>> fp_result;

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b);
bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b);
void share_dataset(char *dataset_path);
void split_dataset(char *dataset_path);
void read_transactions(vector<vector<string>> &transactions, map<string, int> &support_count);
void read_chunk(vector<vector<string>> &transactions, map<string, int> &support_count);
void send_result(vector<string> itemset, const int &support);
void send_end();
void get_results(const int &mpi_num_processes, map<string, int> &results);
void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold);
void pfp_growth(const vector<vector<string>> &transactions, const map<string, int> &support_count, int threshold);
void pfp_growthR(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold);

template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flip_pair<A,B>);
    return dst;
}


int main(int argc, char **argv)
{
    int mpi_rank;
    int mpi_num_processes;
    int threshold (0);
    int word_per_line (INT_MAX);
    unsigned int shortest_itemset (INT_MAX);
    vector<vector<string>> transactions;
    map<string, int> support_count;
    map<string, int> results;
    multimap<int, string> top100;

#ifdef DEBUG
    this_thread::sleep_for(milliseconds(5000));
#endif

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if(argc != 5
            || (threshold = stoi(argv[2])) < 0
            || (word_per_line = stoi(argv[3])) <= 0
            || (shortest_itemset = stoi(argv[4])) <= 0)
            ERR("Usage: ./fpgrowth file threshold wordPerLine shortestItemset");

    if(mpi_rank == MASTER){
        if(threshold == 0) split_dataset(argv[1]);
        else share_dataset(argv[1]);

        get_results(mpi_num_processes, results);
    } else {
        if(threshold == 0) read_chunk(transactions, support_count);
        else read_transactions(transactions, support_count);
        pfp_growth(transactions, support_count, threshold);
    }

    if(MPI_Finalize() != MPI_SUCCESS)
        ERR("MPI_Finalize@main");

    for(pair<string, int> result : results)
        if(result.second > 1)
            cout << result.first << " #" << result.second << endl;

    return 0;
}

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b) {
    return a.second < b.second;
}

bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b) {
    return a.second > b.second;
}

struct CharStream : std::streambuf
{
    CharStream(char* begin, char* end) {
        this->setg(begin, begin, end);
    }
};

struct support_compare{
    support_compare(const map<string, int> &support_reference)
        : support_reference(support_reference) {}

    bool operator() (const string &a, const string &b){
        return support_reference.at(a) > support_reference.at(b);
    }

    map<string, int> support_reference;
};

void share_dataset(char *dataset_path){
    int nbytes;
    char buffer[BUFFER_SIZE];
    ifstream dataset;

    dataset.open(dataset_path);
    if(!dataset.is_open())
        ERR("Opening file@share_dataset");

    while(!dataset.eof()){
        if(dataset.fail())
            ERR("Exceeded buffer size@share_dataset");

        memset(buffer, 0, BUFFER_SIZE);
        dataset.getline(buffer, BUFFER_SIZE - 1);
        nbytes = strlen(buffer);

        if(MPI_Bcast(&nbytes, 1, MPI_INT, MASTER, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Bcast@share_dataset");

        if(MPI_Bcast(buffer, nbytes, MPI_CHAR, MASTER, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Bcast@share_dataset");
    }

    if(dataset.bad())
        ERR("Reading file@share_dataset");

    dataset.close();
}

void split_dataset(char *dataset_path){
    int process_id = MASTER + 1;
    int mpi_size;
    int nbytes;
    char buffer[BUFFER_SIZE];
    ifstream dataset;

    dataset.open(dataset_path);
    if(!dataset.is_open())
        ERR("Opening file@split_dataset");

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    while(!dataset.eof()){
        if(dataset.fail())
            ERR("Exceeded buffer size@split_dataset");

        memset(buffer, 0, BUFFER_SIZE);
        dataset.getline(buffer, BUFFER_SIZE - 1);
        nbytes = strlen(buffer);
        if(nbytes == 0){
            for(int i = MASTER + 1; i < mpi_size; i++)
                if(MPI_Send(&nbytes, 1, MPI_INT, i, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
                    ERR("MPI_Send@split_dataset");

            break;
        }

        if(MPI_Send(&nbytes, 1, MPI_INT, process_id, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Send@split_dataset");

        if(MPI_Send(buffer, nbytes, MPI_CHAR, process_id, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Send@split_dataset");

        process_id = (process_id + 1) % mpi_size;
        if(process_id == MASTER) process_id = MASTER + 1;
    }

    if(dataset.bad())
        ERR("Reading file@split_dataset");

    dataset.close();
}

void read_transactions(vector<vector<string>> &transactions, map<string, int> &support_count){
    int nbytes = INT_MAX;
    char buffer[BUFFER_SIZE];
    istringstream string_stream;
    string item;
    string line;
    vector<string> transaction;

    while(nbytes > 0){
        if(MPI_Bcast(&nbytes, 1, MPI_INT, MASTER, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Bcast@read_transactions");

        memset(buffer, 0, BUFFER_SIZE);
        if(MPI_Bcast(buffer, nbytes, MPI_CHAR, MASTER, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Bcast@read_transactions");

        CharStream char_stream(buffer, buffer + strlen(buffer));
        istream in_stream(&char_stream);
        while(getline(in_stream, line)){
            string_stream.str(line);
            transaction.clear();
            while(string_stream >> item){
                transaction.push_back(item);
                if(support_count.find(item) == support_count.end())
                    support_count[item] = 1;
                else support_count[item]++;
            }

            transactions.push_back(transaction);
            string_stream.clear();
        }
    }
}

void read_chunk(vector<vector<string>> &transactions, map<string, int> &support_count){
    int nbytes = INT_MAX;
    char buffer[BUFFER_SIZE];
    istringstream string_stream;
    string item;
    string line;
    vector<string> transaction;

    while(true){
        if(MPI_Recv(&nbytes, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            ERR("MPI_Recv@read_chunk");

        if(nbytes == 0) break;

        memset(buffer, 0, BUFFER_SIZE);
        if(MPI_Recv(buffer, nbytes, MPI_CHAR, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            ERR("MPI_Recv@read_chunk");

        CharStream char_stream(buffer, buffer + strlen(buffer));
        istream in_stream(&char_stream);
        while(getline(in_stream, line)){
            string_stream.str(line);
            transaction.clear();
            while(string_stream >> item){
                transaction.push_back(item);
                if(support_count.find(item) == support_count.end())
                    support_count[item] = 1;
                else support_count[item]++;
            }

            transactions.push_back(transaction);
            string_stream.clear();
        }
    }
}

void send_result(vector<string> itemset, const int &support){
    int string_length;
    string result;

    sort(itemset.begin(), itemset.end());
    for(string item : itemset)
        result += item + " ";
    result += " " + ESCAPE_CHAR + to_string(support);
    string_length = result.length() + 1;
    if(MPI_Send(&string_length, 1, MPI_INT, MASTER, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
                ERR("MPI_Send@send_result");

    if(MPI_Send(result.c_str(), string_length, MPI_CHAR, MASTER, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
                ERR("MPI_Send@send_result");
}

void send_end(){
    char end[] = END;
    int string_length = strlen(end) + 1;

    if(MPI_Send(&string_length, 1, MPI_INT, MASTER, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
                ERR("MPI_Send@send_result");

    if(MPI_Send(end, string_length, MPI_CHAR, MASTER, NO_TAG, MPI_COMM_WORLD) != MPI_SUCCESS)
                ERR("MPI_Send@send_result");
}

void get_results(const int &mpi_num_processes, map<string, int> &results){
    int num_processes_ended = 0;
    int result_length;
    char result[BUFFER_SIZE];
    string result_string;
    int ended[mpi_num_processes];
    size_t escape_pos;
    string itemset;
    int support;

    memset(ended, false, mpi_num_processes);
    while(num_processes_ended < mpi_num_processes - 1){
        for(int process_id = MASTER + 1; process_id < mpi_num_processes; process_id++){
            if(ended[process_id] == true) continue;

            if(MPI_Recv(&result_length, 1, MPI_INT, process_id, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
                ERR("MPI_Recv@get_results");

            if(MPI_Recv(result, result_length, MPI_CHAR, process_id, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
                ERR("MPI_Recv@get_results");

            if(strcmp(result, END) == 0){
                num_processes_ended++;
                ended[process_id] = true;
            } else {
                result_string = string(result);
                escape_pos = result_string.find_last_of(ESCAPE_CHAR);
                itemset = result_string.substr(0, escape_pos);
                support = stoi(result_string.substr(escape_pos + 1));
                if(results.find(itemset) == results.end()) results[itemset] = support;
                else results[itemset] += support;

                // cout << "Process " << process_id << ": " << result << endl;
            }
        }
    }
}

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
            fp_result.push_back({*tmp_vector, it->second});
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

void pfp_growth(const vector<vector<string>> &transactions, const map<string, int> &support_count, int threshold){
    vector<pair<string, int>> header_table;
    shared_ptr<vector<string>> tmp_vector;
    shared_ptr<fp_tree> ft = make_shared<fp_tree>();
    vector<string> pattern;
    vector<string> tmp_transaction;
    steady_clock::time_point master_begin;
    steady_clock::time_point master_end;
    duration<double> master_duration;

    pattern.push_back(string());

    master_begin = steady_clock::now();

    int item_id = 1;
    int my_rank;
    int mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    mpi_size = (threshold == 0)? 1 : mpi_size - 1;
    int my_item = (threshold == 0)? 1 : my_rank;

    for(map<string, int>::const_iterator it = support_count.begin(); it != support_count.end(); it++, item_id++){
        if(item_id == my_item){
            my_item += mpi_size;
            if(it->second > threshold){
//            tmp_vector = make_shared<vector<string>>();
//            tmp_vector->push_back("");
//            tmp_vector->push_back(it->first);
//            send_result(*tmp_vector, it->second);
            header_table.push_back({it->first, it->second});
            }
        }
    }

    if(header_table.size() == 0){
        send_end();
        return;
    }

#pragma omp parallel for num_threads(N_THREADS) private(tmp_transaction) shared(ft, transactions, support_count)
    for(vector<vector<string>>::const_iterator transaction = transactions.begin(); transaction < transactions.end(); transaction++){
        tmp_transaction = *transaction;
        sort(tmp_transaction.begin(), tmp_transaction.end(), support_compare(support_count));
#pragma omp critical
        ft->insert_transaction(tmp_transaction.begin(), tmp_transaction.end());
    }
    master_end = steady_clock::now();

    master_duration = duration_cast<duration<double>>(master_end - master_begin);

    cout.unsetf (ios::floatfield);
    cout.precision(2);

    pfp_growthR(ft, header_table, pattern, threshold);

    send_end();
}

void pfp_growthR(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold){
    map<string, int> support_count;
    vector<pair<string, int>> new_header_table;
    shared_ptr<vector<string>> new_pattern;
    shared_ptr<fp_tree> new_ft;
    shared_ptr<vector<string>> tmp_vector;
    shared_ptr<vector<vector<string>>> tmp_transactions;

#pragma omp parallel for num_threads(N_THREADS) private(support_count, new_header_table, new_pattern, new_ft, tmp_vector, tmp_transactions) shared(ft, header_table, pattern, threshold)
    for(vector<pair<string, int>>::const_iterator header_row = header_table.begin(); header_row < header_table.end(); header_row++){
            new_pattern = make_shared<vector<string>>(pattern);
            new_pattern->push_back(header_row->first);
            tmp_transactions = ft->get_transaction(header_row->first);

            for(vector<string> transaction : *tmp_transactions)
                for(string item : transaction)
                    if(support_count.find(item) == support_count.end())
                        support_count[item] = 1;
                    else support_count[item]++;

            for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++){
                if(it->second > threshold){
                    tmp_vector = make_shared<vector<string>>(*new_pattern);
                    tmp_vector->push_back(it->first);
#pragma omp critical
                    send_result(*tmp_vector, it->second);
                    new_header_table.push_back({it->first, it->second});
               }
            }

            sort(new_header_table.begin(), new_header_table.end(), pair_compare);

            new_ft = make_shared<fp_tree>();
            for(vector<string> transaction : *tmp_transactions){
                sort(transaction.begin(), transaction.end(), support_compare(support_count));
                new_ft->insert_transaction(transaction.begin(), transaction.end());
            }

            pfp_growthR(new_ft, new_header_table, *new_pattern, threshold);
            new_header_table.clear();
            support_count.clear();
    }
}


