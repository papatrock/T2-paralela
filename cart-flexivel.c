#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
    int rank, procs;
    MPI_Comm cart_comm;
    int reorder;
    int coord[2], id;
    int up, down, left, right;
    // ...
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    int ndims = 2;
    int dims[2] = {0, 0}; // Zeros significam "deixe o MPI decidir"
    int periodical[2] = {0, 0};
    reorder = 1;

    // Pede ao MPI para encontrar a melhor grade 2D para 'procs' processos
    MPI_Dims_create(procs, ndims, dims);

    // Se procs=12, dims será [4, 3] (ou [3,4], etc.)
    // Se procs=2, dims será [2, 1] (ou [1,2])
    if (rank == 0) {
        printf("Para %d processos, o MPI criou uma grade de %d x %d\n", procs, dims[0], dims[1]);
    }

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodical, reorder, &cart_comm);
    // ... resto do código ...
}