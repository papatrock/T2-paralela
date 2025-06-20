#include <mpi.h>
#include <string.h>
#include <stdio.h>

#define STD_TAG 0
int main(int argc, char **argv) {
    int i, my_rank, n_procs; char msg[100]; MPI_Status status;
    printf("isso sera imprimido duas vezes\n"); // isso sera imprimido duas vezes

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
        printf("e isso só uma\n");
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    if (my_rank != 0) {
        sprintf(msg, "I’m alive!");
    MPI_Send(msg, strlen(msg) + 1, MPI_CHAR, 0, STD_TAG, MPI_COMM_WORLD);
    } else {
        for (i = 1; i < n_procs; i++) {
            MPI_Recv(msg, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
            printf("Proc %d: %s \n", status.MPI_SOURCE, msg);
        }
    }

    MPI_Finalize();

    return 0;
}
