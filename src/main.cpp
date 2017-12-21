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

#define N_THREADS       8
#define BUFFER_SIZE     8192
#define MASTER          0
#define NO_TAG          0
#define END             "***END***"
#define DEBUG

#define ERR(s){cerr << "ERROR " << s << endl; exit(EXIT_FAILURE);}
#define ECHO(a){cout << a << endl;}

using namespace std;
using namespace std::chrono;

vector<pair<vector<string>, int>> fp_result;
vector<pair<vector<string>, int>> pfp_result;

bool pair_compare (const pair<string, int> &a, const pair<string, int> &b);
bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b);
void share_dataset(char *dataset_path);
void read_transactions(vector<vector<string>> &transactions);
void send_result(const vector<string> &itemset, const int &support);
void send_end();
void get_results(const int &mpi_num_processes);
void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold);
void master(const vector<vector<string>> &transactions, int threshold);
void pfp_growth(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold);


int main(int argc, char **argv)
{
    int mpi_rank;
    int mpi_num_processes;
    int threshold (0);
    int word_per_line (INT_MAX);
    unsigned int shortest_itemset (INT_MAX);
//    int i (0);
    ifstream dataset;
    istringstream iss;
    string item;
    string line;
    vector<vector<string>> transactions;
    vector<string> tmp;
    steady_clock::time_point fp_begin;
    steady_clock::time_point fp_end;
    steady_clock::time_point pfp_begin;
    steady_clock::time_point pfp_end;
    duration<double> fp_duration;
    duration<double> pfp_duration;

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
            share_dataset(argv[1]);
            get_results(mpi_num_processes);
    } else {
        read_transactions(transactions);
        master(transactions, threshold);
    }

    if(MPI_Finalize() != MPI_SUCCESS)
        ERR("MPI_Finalize@main");

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
        : support_reference(support_reference) {};

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

void read_transactions(vector<vector<string>> &transactions){
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
            while(string_stream >> item)
                transaction.push_back(item);

            transactions.push_back(transaction);
            string_stream.clear();
        }
    }
}

void send_result(const vector<string> &itemset, const int &support){
    int string_length;
    string result;
    for(string item : itemset)
        result += " " + item;

    result += " [" + to_string(support) + "]";
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

void get_results(const int &mpi_num_processes){
    int num_processes_ended = 0;
    int result_length;
    char result[BUFFER_SIZE];
    int ended[mpi_num_processes];

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
            } else cout << "Process " << process_id << ": " << result << endl;
        }
    }
}

void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold){
    map<string, int> support_count;
    vector<pair<string, int>> header_table;
    shared_ptr<vector<string>> tmp_vector;
    shared_ptr<fp_tree> ft = make_shared<fp_tree>();
    shared_ptr<vector<vector<string>>> tmp_transactions;
    //vector<pair<vector<string>, int>> *tmp_result = &result;

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

void master(const vector<vector<string>> &transactions, int threshold){
    map<string, int> support_count;
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
    for(vector<string> transaction : transactions)
        for(string item : transaction)
            if(support_count.find(item) == support_count.end())
                support_count[item] = 1;
            else support_count[item]++;

    int item_id = 1;
    int my_rank;
    int mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    mpi_size -= 1;
    int my_item = my_rank;
    for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++, item_id++){
        if(item_id == my_item){
            my_item += mpi_size;
            if(it->second > threshold){
            tmp_vector = make_shared<vector<string>>();
            tmp_vector->push_back("");
            tmp_vector->push_back(it->first);
            //pfp_result.push_back({*tmp_vector, it->second});
            send_result(*tmp_vector, it->second);
            header_table.push_back({it->first, it->second});
            }
        }
    }

    if(header_table.size() == 0){
        send_end();
        return;
    }

    sort(header_table.begin(), header_table.end(), pair_compare);

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

    pfp_growth(ft, header_table, pattern, threshold);

    send_end();
}

void pfp_growth(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold){
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
                    //pfp_result.push_back({*tmp_vector, it->second});
                    new_header_table.push_back({it->first, it->second});
               }
            }

            sort(new_header_table.begin(), new_header_table.end(), pair_compare);

            new_ft = make_shared<fp_tree>();
            for(vector<string> transaction : *tmp_transactions){
                sort(transaction.begin(), transaction.end(), support_compare(support_count));
                new_ft->insert_transaction(transaction.begin(), transaction.end());
            }

            pfp_growth(new_ft, new_header_table, *new_pattern, threshold);
            new_header_table.clear();
            support_count.clear();
    }
}


