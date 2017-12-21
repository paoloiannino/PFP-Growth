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

#include "pfp_growth.h"

//#define PRINT_RESULT

using namespace std;
using namespace std::chrono;

int main(int argc, char **argv)
{
    int mpi_rank;
    int mpi_num_processes;
    int threshold (0);
    vector<vector<string>> transactions;
    map<string, int> support_count;
    map<string, int> results;
    steady_clock::time_point begin;
    steady_clock::time_point end;
    duration<double> pfp_duration;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if(argc != 3 || (threshold = stoi(argv[2])) < 0)
            ERR("Usage: ./pfp-growth file threshold");

    if(mpi_rank == MASTER){
        begin = steady_clock::now();
        if(threshold == 0) split_dataset(argv[1]);
        else share_dataset(argv[1]);

        get_results(mpi_num_processes, results);
        end = steady_clock::now();
    } else {
        if(threshold == 0) read_chunk(transactions, support_count);
        else read_transactions(transactions, support_count);
        pfp_growth(transactions, support_count, threshold);
    }

#ifdef PRINT_RESULT
    for(pair<string, int> result : results)
        if(result.second > 1)
            cout << result.first << " #" << result.second << endl;
#endif

    if(mpi_rank == MASTER){
        pfp_duration = duration_cast<duration<double>>(end - begin);
        cout.unsetf (ios::floatfield);
        cout.precision(5);
        cout << "TIME: " << pfp_duration.count() << endl;
    }

    if(MPI_Finalize() != MPI_SUCCESS)
        ERR("MPI_Finalize@main");

    return 0;
}




