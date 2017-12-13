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

#include "fp_tree.h"

#define N_THREADS 8

using namespace std;
using namespace std::chrono;


int main(int argc, char **argv)
{
    int num_processes, id;
    MPI_Init(&argc, &argv); /* MPI initialization, all the processes copy of the current variables */
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes); /* getting number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &id); /* Get the rank (id) of the process */

    if (id == 0) {
        cout << "master" << endl;
    } else {
        cout << "slave" << endl;
    }

    MPI_Finalize(); /* Finalize the MPI environment*/

    return 0;

}
