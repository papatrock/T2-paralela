#include <stdio.h>
#include <mpi.h>

#define STD_TAG 0

static long num_steps = 100000;
double step;

int main(int argc, char *argv[])
{

    int i;
    double x, pi, sum = 0.0;
    int my_rank, n_procs;
    // ideia, vetor de somas parciais de tamanho n_procs
    // cada um soma no seu indice
    MPI_Status status;

    MPI_Init(&argc, &argv);
    // devolve em (myrank) a posição do proc
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // devolve em (nprocs) o total de proces
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    double minha_soma = 0.0;

    step = 1.0 / (double)num_steps;

    for (i = my_rank; i < num_steps; i += n_procs)
    {
        x = (i + 0.5) * step;
        minha_soma += 4.0 / (1.0 + x * x);
    }
    // rank 0 espera o resultado das outras
    if (my_rank == 0)
    {
        double soma_total = minha_soma;
        double soma_recebida;

        for (i = 1; i < n_procs; ++i)
        {
            // MPI_Recv(buffer, contagem, tipo, origem, tag, comunicador, status)
            MPI_Recv(&soma_recebida, 1, MPI_DOUBLE, i, STD_TAG, MPI_COMM_WORLD, &status); // status?
            soma_total += soma_recebida;
            printf("Mestre recebeu %f do processo %d\n", soma_recebida, i);
        }
        pi = step * soma_total;
        printf("Pi : %lf\n", pi);
    }
    else
    {
        // MPI_Send(buffer, contagem, tipo, destino, tag, comunicador)
        MPI_Send(&minha_soma, 1, MPI_DOUBLE, 0, STD_TAG, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}