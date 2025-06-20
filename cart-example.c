#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
    int rank, procs;
    MPI_Comm cart_comm;
    int reorder;
    int coord[2], id;
    // As variaveis up, down, etc. não foram usadas, podem ser removidas por enquanto
    // int up, down, left, right;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    // Esta parte agora está correta para uma execução com 2 processos
    if(procs != 2) {
        if (rank == 0) { // Apenas o rank 0 imprime o erro
            printf("Erro: Este programa requer 2 processos. Finalizando.\n");
        }
        MPI_Finalize(); // Finaliza corretamente mesmo em caso de erro
        exit(1);
    }

    int dim[2];
    dim[0] = 1; //linhas
    dim[1] = 2; //colunas

    int periodical[2] = {1,1};
    reorder = 0;

    MPI_Cart_create(MPI_COMM_WORLD,2,dim,periodical, reorder, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coord);

    // Você poderia adicionar um printf aqui para ver as coordenadas de cada rank
    printf("Sou o rank %d e minhas coordenadas na grade são (%d, %d)\n", rank, coord[0], coord[1]);

    MPI_Comm_free(&cart_comm); // Boa prática liberar o comunicador


    MPI_Finalize();

    return 0;
}