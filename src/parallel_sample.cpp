# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>
# include <vector>

using namespace std;

int main () {

    MPI_Init(NULL, NULL); /* MPI initialization */

    /* getting number of processes */
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* Get the rank (id) of the process */
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);


    /* Send / Receive pattern */
    /* Send / Receive pattern */
    /*
        This pattern in our case should be implemented using
        1. Fix data exchanges
        2. Same computation for each node
    */

    if (id == 0) {
        int v1[] = {1,2,3,4,5,6,7,8,9};
        /* Parameters: Data, Data Size, Data Type, Destination, Tag, Communicator */
        /* NB: Tag can contain additional information, in our case "Special Key" */
        MPI_Send(v1, sizeof(v1), MPI_INT, 1, 0, MPI_COMM_WORLD);
        cout << "Process " << id << " sent: " << endl;
        for (const auto& i : v1) {
            std::cout << i << std::endl;
        }
        int remote_sum;
        MPI_Recv(&remote_sum, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Process " << id << " REMOTE sum: " << remote_sum <<endl;
    } else if (id == 1) {
        int v1[9];
        /* Parameters: Data, Data Size, Data Type, Source, Tag, Communicator, Status */
        /* NB: hardcoded 9 as size, should be changed. */
        MPI_Recv(v1, sizeof(v1), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Process " << id << " received:" << endl;
        for (const auto& i : v1) {
            std::cout << i << std::endl;
        }
        /* Computation can be done for the data received */
        int sum = 0;
        for (const auto& i : v1) {
            sum += i;
        }

        cout << "Process " << id << " computed sum: " << sum <<endl;

        /* Returning the result */
        MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    }


    /* Finalize the MPI environment*/
    MPI_Finalize();
    return 0;
}
