#include <mpi.h>
#include <stdio.h>
#include <unistd.h> // <--- Adicione este cabeçalho para a função gethostname

int main(int argc, char **argv) {
    int my_rank, n_procs;
    char hostname[256]; // Buffer para guardar o nome da máquina

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    // Pega o nome da máquina atual
    gethostname(hostname, sizeof(hostname));

    // Imprime o rank E o nome da máquina
    printf("Sou o rank %d e estou rodando no host: %s\n", my_rank, hostname);

    // ... (resto da sua lógica MPI, se tiver) ...

    MPI_Finalize();
    return 0;
}