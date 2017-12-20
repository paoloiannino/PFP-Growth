#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <map>
#include <vector>
#include <set>
#include <climits>      // For INT_MAX
//#include <omp.h>
#include <chrono>
#include <mpi.h>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi.hpp>

//#include "fp_tree.h"

#define N_THREADS 8

using namespace std;
using namespace std::chrono;
using namespace boost::mpi;

// /* This structure contains the payload of the package broadcasted to all nodes with MPI */
// struct package {
//  fp_tree ft;
//  vector<pair<string, int>> header_table;
//  vector<string> pattern;
//  int threshold;
//  int num_machines;
// } ;

// vector<pair<vector<string>, int>> fp_result;
// vector<pair<vector<string>, int>> pfp_result;
//
// vector<vector<string>> read_file(char* path, int word_per_line);
// bool pair_compare (const pair<string, int> &a, const pair<string, int> &b);
// bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b);
// void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold);
// void master(const vector<vector<string>> &transactions, int threshold);
// void slave();
// void pfp_growth(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold);
// void mpi_pfp_growth ( shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table, const vector<string> &pattern, int threshold );

void share_dataset()

int main(int argc, char **argv)
{
    int threshold (0);
     int word_per_line (INT_MAX);
     unsigned int shortest_itemset (INT_MAX);
     int i (0);
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

   int num_processes, id;

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
        tmp.clear();
        i = 0;
        while(iss >> item){
            tmp.push_back(item);
            i++;

            if(i == word_per_line) break;
        }

        transactions.push_back(tmp);
        iss.clear();
    }

   MPI_Init(&argc, &argv); /* MPI initialization, all the processes copy of the current variables */
   MPI_Comm_size(MPI_COMM_WORLD, &num_processes); /* getting number of processes */
   MPI_Comm_rank(MPI_COMM_WORLD, &id); /* Get the rank (id) of the process */

  // if (world.rank() == 0) {

       //int size_payload;
       //MPI_Recv(&size_payload, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   //     vector<string> pattern;
   //     pattern.push_back(string());
   //
   //     fp_begin = steady_clock::now();
   //     fp_growth(transactions, pattern, threshold);
   //     fp_end = steady_clock::now();
   //     fp_duration = duration_cast<duration<double>>(fp_end - fp_begin);
   //
   //     pfp_begin = steady_clock::now();
   //     master(transactions, threshold); /* this method send the data */
   //
   //     for (int i=1; i <= num_processes; i++)
   //     {
   //         //RECEIVE Size
   //         //RECEIVE actual data
   //     }
   //
   //     pfp_end = steady_clock::now();
   //     pfp_duration = duration_cast<duration<double>>(pfp_end - pfp_begin);
   //
   //     sort(fp_result.begin(), fp_result.end(), pair_compare_result);
   //     sort(pfp_result.begin(), pfp_result.end(), pair_compare_result);
   //
   //     cout.unsetf (ios::floatfield);
   //     cout.precision(2);
   //     cout << "FP-Growth: " << fp_duration.count() << " seconds." << endl;
   //     for(pair<vector<string>, int> pattern : fp_result){
   //         if(pattern.second < threshold
   //                 || pattern.first.size() < (shortest_itemset + 1)) continue;
   //
   //         sort(pattern.first.begin(), pattern.first.end());
   //         for(string item : pattern.first)
   //             cout << item << " ";
   //
   //         cout << "SUPPORT: " << pattern.second << endl;
   //     }
   //
   //     cout << endl;
   //
   //     cout << "PFP-Growth [" << N_THREADS << " threads]: " << pfp_duration.count() << " seconds." << endl;
   //
   //     for(pair<vector<string>, int> pattern : pfp_result){
   //         if(pattern.second < threshold
   //                 || pattern.first.size() < (shortest_itemset + 1)) continue;
   //
   //         sort(pattern.first.begin(), pattern.first.end());
   //         for(string item : pattern.first)
   //             cout << item << " ";
   //
   //         cout << "SUPPORT: " << pattern.second << endl;
   //     }
   //
   //     cout << endl;
   //     cout << "PFP-Growth is " << ((fp_duration.count() - pfp_duration.count()) / fp_duration.count()) * 100 << "% faster" << endl;
   //     cout << endl;
   //
   //
   // } else {
   //     slave(); //this method will receive the data from the master node and do the computation
   //
   // }

   cout << "Process: " << id << endl;
   for(vector<string> transaction : transactions){
        for(string item : transaction)
            cout << item << " ";

        cout << endl << endl;
    }

   MPI_Finalize(); /* Finalize the MPI environment*/

   return 0;

}

// vector<vector<string>> read_file(char* path, int word_per_line){
//
//    istringstream iss;
//    ifstream dataset;
//    vector<vector<string>> transactions;
//    vector<string> tmp;
//    string item;
//    string line;
//    int i (0);
//
//    dataset.open(path);
//
//    if(!dataset.is_open()){
//        cerr << "ERROR: opening file" << endl;
//        return transactions;
//    }
//
//
//    while(getline(dataset, line)){
//        iss.str(line);
//        tmp.clear();
//        i = 0;
//        while(iss >> item){
//            tmp.push_back(item);
//            i++;
//
//            if(i == word_per_line) break;
//        }
//
//        transactions.push_back(tmp);
//        iss.clear();
//    }
//
//    dataset.close();
//
//    return transactions;
// }
//
// bool pair_compare (const pair<string, int> &a, const pair<string, int> &b) {
//    return a.second < b.second;
// }
//
// bool pair_compare_result (pair<vector<string>, int> a, pair<vector<string>, int> b) {
//    return a.second > b.second;
// }
//
// struct support_compare{
//    support_compare(const map<string, int> &support_reference)
//        : support_reference(support_reference) {};
//
//    bool operator() (const string &a, const string &b){
//        return support_reference.at(a) > support_reference.at(b);
//    }
//
//    map<string, int> support_reference;
// };
//
// void fp_growth(const vector<vector<string>> &transactions, const vector<string> &pattern, int threshold){
//    map<string, int> support_count;
//    vector<pair<string, int>> header_table;
//    shared_ptr<vector<string>> tmp_vector;
//    shared_ptr<fp_tree> ft = make_shared<fp_tree>();
//    shared_ptr<vector<vector<string>>> tmp_transactions;
//    //vector<pair<vector<string>, int>> *tmp_result = &result;
//
//    for(vector<string> transaction : transactions)
//        for(string item : transaction)
//            if(support_count.find(item) == support_count.end())
//                support_count[item] = 1;
//            else support_count[item]++;
//
//    for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++){
//        header_table.push_back({it->first, it->second});
//        tmp_vector = make_shared<vector<string>>(pattern);
//        tmp_vector->push_back(it->first);
//
//        if(it->second > threshold)
//            fp_result.push_back({*tmp_vector, it->second});
//    }
//
//    sort(header_table.begin(), header_table.end(), pair_compare);
//
//    for(vector<string> transaction : transactions){
//        sort(transaction.begin(), transaction.end(), support_compare(support_count));
//        ft->insert_transaction(transaction.begin(), transaction.end());
//    }
//
//    for(pair<string, int> header_row : header_table){
//        if(header_row.second > threshold){
//            tmp_vector = make_shared<vector<string>>(pattern);
//            tmp_vector->push_back(header_row.first);
//            tmp_transactions = ft->get_transaction(header_row.first);
//            fp_growth(*tmp_transactions, *tmp_vector, threshold);
//        }
//    }
//
//    return;
// }
//
// void master(const vector<vector<string>> &transactions, int threshold){
//    map<string, int> support_count;
//    vector<pair<string, int>> header_table;
//    shared_ptr<vector<string>> tmp_vector;
//    shared_ptr<fp_tree> ft = make_shared<fp_tree>();
//    vector<string> pattern;
//    pattern.push_back(string());
//
//    for(vector<string> transaction : transactions)
//        for(string item : transaction)
//            if(support_count.find(item) == support_count.end())
//                support_count[item] = 1;
//            else support_count[item]++;
//
//    for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++){
//        if(it->second > threshold){
//            tmp_vector = make_shared<vector<string>>();
//            tmp_vector->push_back("");
//            tmp_vector->push_back(it->first);
//            pfp_result.push_back({*tmp_vector, it->second});
//            header_table.push_back({it->first, it->second});
//        }
//    }
//
//    sort(header_table.begin(), header_table.end(), pair_compare);
//
//    for(vector<string> transaction : transactions){
//        sort(transaction.begin(), transaction.end(), support_compare(support_count));
//        ft->insert_transaction(transaction.begin(), transaction.end());
//    }
//    /* MPI parallelization on the wrapper function mpi_pfp_growth */
//    size_t size_ht = 0;
//    for (auto e : header_table){
//        size_ht += e.first.size();
//    }
//
//    size_ht += sizeof(vector<pair<string, int>>);
//    size_ht += sizeof(pair<string, int>) * header_table.size();
//    size_ht += sizeof(int) * header_table.size();
//
//    MPI_Send(&header_table, size_ht, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
//
//    cout << size_ht << endl;
//
//
//    //send_mpi(ft, header_table, pattern, threshold);
//
// }
//
// void slave(){
//    /* Split the header table in world_size parts and get the id^th one */
//    /* WARNING: in this case we should use the world_size-1 bacause world_size contains master process */
//    //slave() //this method will receive the data from the master node and do the computation
//    /*size_t const my_header_size = header_table.size() / (num_processes - 1);
//    vector<pair<string, int>> my_header_table ( header_table.begin() + (my_header_size * (id - 1)), header_table.begin() + (my_header_size * id) );
//    for (auto e : header_table){
//
//    }
//    for (const auto& i : my_header_table) {
//        std::cout << i.first << " - " << i.second << std::endl;
//    }
//    pfp_growth (ft, my_header_table, pattern, threshold);*/
//    /* Returning the result
//    MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);*/
//    cout << "Slave" << endl;
// }
//
// void pfp_growth(shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table,  const vector<string> &pattern, int threshold){
//    map<string, int> support_count;
//    vector<pair<string, int>> new_header_table;
//    shared_ptr<vector<string>> new_pattern;
//    shared_ptr<fp_tree> new_ft;
//    shared_ptr<vector<string>> tmp_vector;
//    shared_ptr<vector<vector<string>>> tmp_transactions;
//
//
// //    for(pair<string, int> header_row : header_table){
// //            new_pattern = make_shared<vector<string>>(pattern);
// //            new_pattern->push_back(header_row.first);
// //            tmp_transactions = ft->get_transaction(header_row.first);
//
//
// //#pragma omp parallel for num_threads(N_THREADS) private(support_count, new_header_table, new_pattern, new_ft, tmp_vector, tmp_transactions) shared(ft, header_table, pattern, threshold)
//    for(vector<pair<string, int>>::const_iterator header_row = header_table.begin(); header_row < header_table.end(); header_row++){
//            new_pattern = make_shared<vector<string>>(pattern);
//            new_pattern->push_back(header_row->first);
//            tmp_transactions = ft->get_transaction(header_row->first);
//
//            for(vector<string> transaction : *tmp_transactions)
//                for(string item : transaction)
//                    if(support_count.find(item) == support_count.end())
//                        support_count[item] = 1;
//                    else support_count[item]++;
//
//            for(map<string, int>::iterator it = support_count.begin(); it != support_count.end(); it++){
//                if(it->second > threshold){
//                    tmp_vector = make_shared<vector<string>>(*new_pattern);
//                    tmp_vector->push_back(it->first);
// //#pragma omp critical
//                    pfp_result.push_back({*tmp_vector, it->second}); // NB. send mpi NO write on report
//                    new_header_table.push_back({it->first, it->second});
//               }
//            }
//
//            sort(new_header_table.begin(), new_header_table.end(), pair_compare);
//
//            new_ft = make_shared<fp_tree>();
//            for(vector<string> transaction : *tmp_transactions){
//                sort(transaction.begin(), transaction.end(), support_compare(support_count));
//                new_ft->insert_transaction(transaction.begin(), transaction.end());
//            }
//
//            pfp_growth(new_ft, new_header_table, *new_pattern, threshold);
//            new_header_table.clear();
//            support_count.clear();
//    }
// }
//
//
// void mpi_pfp_growth ( shared_ptr<fp_tree> ft, const vector<pair<string, int>> &header_table, const vector<string> &pattern, int threshold ){
//    int num_processes, id;
//    MPI_Init(NULL, NULL); /* MPI initialization, all the processes copy of the current variables */
//    MPI_Comm_size(MPI_COMM_WORLD, &num_processes); /* getting number of processes */
//    MPI_Comm_rank(MPI_COMM_WORLD, &id); /* Get the rank (id) of the process */
//
//    if (id == 0) {
//        /*
//            We don't need to broadcast anything:
//            - All the processes will get the copy of the variables in the moment of MPI_Init
//            - Each process should execute the same function and return the result to the main
//            - The main machine will only be in charge of receive the partial result of the machines, merge and sort them
//        */
//
//        /* Need to define a well-known variable that contain the results */
//        int size_payload;
//        MPI_Recv(&size_payload, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//    } else {
//        /* Split the header table in world_size parts and get the id^th one */
//        /* WARNING: in this case we should use the world_size-1 bacause world_size contains master process */
//        size_t const my_header_size = header_table.size() / (num_processes - 1);
//        vector<pair<string, int>> my_header_table ( header_table.begin() + (my_header_size * (id - 1)), header_table.begin() + (my_header_size * id) );
//        for (const auto& i : my_header_table) {
//            std::cout << i.first << " - " << i.second << std::endl;
//        }
//        pfp_growth (ft, my_header_table, pattern, threshold);
//        /* Returning the result
//        MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);*/
//
//    }
//
//    MPI_Finalize(); /* Finalize the MPI environment*/
// }
