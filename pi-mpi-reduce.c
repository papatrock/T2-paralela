#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    long num_steps = 100000;
    double step;
    int i, my_rank, n_procs;
    double x, pi, minha_soma = 0.0, soma_total = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    // Cada processo calcula sua parte do trabalho (igual a antes)
    step = 1.0 / (double)num_steps;
    for (i = my_rank; i < num_steps; i += n_procs)
    {
        x = (i + 0.5) * step;
        minha_soma += 4.0 / (1.0 + x * x);
    }

    // Comunicação: usa MPI_Reduce para somar todas as 'minha_soma'
    // e guardar o resultado em 'soma_total' no processo 0.
    // MPI_Reduce(&valor_de_envio, &valor_de_recebimento, contagem, tipo, operação, raiz, comunicador)
    MPI_Reduce(&minha_soma, &soma_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Apenas o processo raiz (0) calcula e imprime o resultado final
    if (my_rank == 0)
    {
        pi = step * soma_total;
        printf("Pi é aproximadamente %.16f\n", pi);
    }

    MPI_Finalize();
    return 0;
}