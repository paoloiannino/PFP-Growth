#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <map>
#include <vector>
#include <set>
#include <string.h>
#include <stdio.h>
#include <climits>      // For INT_MAX
#include <omp.h>
#include <chrono>
#include <mpi.h>

#define N_THREADS       8
#define BUFFER_SIZE     8192
#define MASTER          0
#define NO_TAG          0

#define ERR(s){cerr << "ERROR " << s << endl; exit(EXIT_FAILURE);}
#define ECHO(a){cout << a << endl;}
#define OK cout << "OK" << endl;

using namespace std;
using namespace std::chrono;

struct CharStream : std::streambuf
{
    CharStream(char* begin, char* end) {
        this->setg(begin, begin, end);
    }
};

void share_dataset(char *dataset_path);

void read_transactions(vector<vector<string>> &transactions);

int main(int argc, char **argv)
{
    int mpi_rank, mpi_num_processes;
    int threshold (0);
    int word_per_line (INT_MAX);
    unsigned int shortest_itemset (INT_MAX);
    vector<vector<string>> transactions;
    MPI_Init(&argc, &argv); /* MPI initialization, all the processes copy of the current variables */
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_processes); /* getting number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); /* Get the rank (id) of the process */

    if(argc != 5
        || (threshold = stoi(argv[2])) < 0
        || (word_per_line = stoi(argv[3])) <= 0
        || (shortest_itemset = stoi(argv[4])) <= 0)
        ERR("Usage: ./fpgrowth file threshold wordPerLine shortestItemset");

    if(mpi_rank == MASTER)
        share_dataset(argv[1]);
    else read_transactions(transactions);

    //cout << "Process: " << mpi_rank << endl;

    MPI_Finalize(); /* Finalize the MPI environment*/

    return 0;
}

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

    return;
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
            ERR("MPI_Recv@share_dataset");

        memset(buffer, 0, BUFFER_SIZE);
        if(MPI_Bcast(buffer, nbytes, MPI_CHAR, MASTER, MPI_COMM_WORLD) != MPI_SUCCESS)
            ERR("MPI_Recv@share_dataset");

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

    return;
}
