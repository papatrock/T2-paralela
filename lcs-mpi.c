#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char *read_seq(char *fname)
{
    FILE *fseq = fopen(fname, "rt");
    if (fseq == NULL)
    {
        perror("Erro ao abrir o arquivo");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    fseek(fseq, 0L, SEEK_END);
    long size = ftell(fseq);
    rewind(fseq);

    char *seq = (char *)calloc(size + 1, sizeof(char));
    if (seq == NULL)
    {
        perror("Erro ao alocar memoria para sequencia");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int i = 0;
    int c;
    while ((c = fgetc(fseq)) != EOF)
    {
        if (c != '\n' && c != '\r')
        {
            seq[i++] = c;
        }
    }
    seq[i] = '\0';

    fclose(fseq);
    return seq;
}

// substiruido pelo alocar matriz local
mtype **allocateScoreMatrix(int sizeA, int sizeB)
{
    int i;
    // Allocate memory for LCS score matrix
    mtype **scoreMatrix = (mtype **)malloc((sizeB + 1) * sizeof(mtype *));
    for (i = 0; i < (sizeB + 1); i++)
        scoreMatrix[i] = (mtype *)malloc((sizeA + 1) * sizeof(mtype));
    return scoreMatrix;
}

// TODO modificar alocação da matriz para o metodo mazieiro
mtype **alocar_matriz_local(int n_rows, int n_cols)
{
    // n_rows e n_cols + 2 : coluna e linhas que vão receber (ou que são 0)
    mtype **matrix = (mtype **)malloc((n_rows + 2) * sizeof(mtype *));
    if (matrix == NULL)
    {
        perror("Falha ao  alocar matriz local");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Aloca a memória para cada linha da matriz
    for (int i = 0; i < n_rows + 2; i++)
    {

        matrix[i] = (mtype *)calloc((n_cols + 2), sizeof(mtype));
        if (matrix[i] == NULL)
        {
            perror("Falha ao  alocar matriz local");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
    return matrix;
}

void initScoreMatrix(mtype **scoreMatrix, int sizeA, int sizeB)
{
    int i, j;
    // Fill first line of LCS score matrix with zeroes
    for (j = 0; j < (sizeA + 1); j++)
        scoreMatrix[0][j] = 0;

    // Do the same for the first collumn
    for (i = 1; i < (sizeB + 1); i++)
        scoreMatrix[i][0] = 0;
}

mtype LCS(MPI_Comm comm_world, int my_rank, int n_procs, char *seqA, int sizeA, char *seqB, int sizeB)
{
    // setup chart
    MPI_Comm comm_cart;
    int ndims = 2, dims[2] = {2, 3}, periods[2] = {0, 0}, reorder = 1;
    MPI_Cart_create(comm_world, ndims, dims, periods, reorder, &comm_cart);
    if (comm_cart == MPI_COMM_NULL)
        return 0;

    int my_cart_rank, my_coords[2];
    MPI_Comm_rank(comm_cart, &my_cart_rank);
    MPI_Cart_coords(comm_cart, my_cart_rank, ndims, my_coords);

    int rank_up, rank_down, rank_left, rank_right;
    MPI_Cart_shift(comm_cart, 0, -1, &rank_down, &rank_up);
    MPI_Cart_shift(comm_cart, 1, -1, &rank_right, &rank_left);

    // alocacao local
    int size_A_borda = sizeA + 1;
    int size_B_borda = sizeB + 1;
    int start_col = (my_coords[1] * size_A_borda) / dims[1];
    int end_col = ((my_coords[1] + 1) * size_A_borda) / dims[1] - 1;
    int local_cols = end_col - start_col + 1;
    int start_row = (my_coords[0] * size_B_borda) / dims[0];
    int end_row = ((my_coords[0] + 1) * size_B_borda) / dims[0] - 1;
    int local_rows = end_row - start_row + 1;
    mtype **matriz_local = alocar_matriz_local(local_rows, local_cols);

    // int total_diagonals = dims[0] + dims[1] - 1;

    for (int bi = 0; bi < dims[0]; bi++)
    {
        for (int bj = 0; bj < dims[1]; bj++)
        {
            // Verifica se EU sou o processo responsável por este bloco
            if (my_coords[0] == bi && my_coords[1] == bj)
            {
                // 1. RECEBE dados se não estiver na borda global
                if (rank_up != MPI_PROC_NULL)
                {
                    MPI_Recv(&matriz_local[0][1], local_cols, MPI_UNSIGNED_SHORT, rank_up, 0, comm_cart, MPI_STATUS_IGNORE);
                }
                if (rank_left != MPI_PROC_NULL)
                {
                    mtype temp_col_buffer[local_rows];
                    MPI_Recv(temp_col_buffer, local_rows, MPI_UNSIGNED_SHORT, rank_left, 1, comm_cart, MPI_STATUS_IGNORE);
                    for (int i = 0; i < local_rows; ++i)
                        matriz_local[i + 1][0] = temp_col_buffer[i];
                }

                // 2. CALCULA o bloco local
                for (int i = 1; i <= local_rows; i++)
                {
                    for (int j = 1; j <= local_cols; j++)
                    {
                        int global_j = start_col + j - 1;
                        int global_i = start_row + i - 1;
                        if (global_j < sizeA && global_i < sizeB && seqA[global_j] == seqB[global_i])
                        {
                            matriz_local[i][j] = matriz_local[i - 1][j - 1] + 1;
                        }
                        else
                        {
                            matriz_local[i][j] = max(matriz_local[i - 1][j], matriz_local[i][j - 1]);
                        }
                    }
                }

                // 3. ENVIA dados se não estiver na borda global
                if (rank_down != MPI_PROC_NULL)
                {
                    MPI_Send(&matriz_local[local_rows][1], local_cols, MPI_UNSIGNED_SHORT, rank_down, 0, comm_cart);
                }
                if (rank_right != MPI_PROC_NULL)
                {
                    mtype temp_col_buffer[local_rows];
                    for (int i = 0; i < local_rows; ++i)
                        temp_col_buffer[i] = matriz_local[i + 1][local_cols];
                    MPI_Send(temp_col_buffer, local_rows, MPI_UNSIGNED_SHORT, rank_right, 1, comm_cart);
                }
            }
        }
    }
    // ... O resto da função (coleta de resultado e limpeza) continua igual ...
    mtype final_score = 0;
    if (my_rank == 0)
    {
        if (my_coords[0] == dims[0] - 1 && my_coords[1] == dims[1] - 1)
        {
            final_score = matriz_local[local_rows][local_cols];
        }
        else
        {
            int rank_do_ultimo_proc;
            int coords_do_ultimo[2] = {dims[0] - 1, dims[1] - 1};
            MPI_Cart_rank(comm_cart, coords_do_ultimo, &rank_do_ultimo_proc);
            MPI_Recv(&final_score, 1, MPI_UNSIGNED_SHORT, rank_do_ultimo_proc, 999, comm_world, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        if (my_coords[0] == dims[0] - 1 && my_coords[1] == dims[1] - 1)
        {
            final_score = matriz_local[local_rows][local_cols];
            MPI_Send(&final_score, 1, MPI_UNSIGNED_SHORT, 0, 999, comm_world);
        }
    }

    for (int i = 0; i < local_rows + 2; i++)
        free(matriz_local[i]);
    free(matriz_local);
    MPI_Comm_free(&comm_cart);
    return final_score;
}

void printMatrix(char *seqA, char *seqB, mtype **scoreMatrix, int sizeA,
                 int sizeB)
{
    int i, j;

    // print header
    printf("Score Matrix:\n");
    printf("========================================\n");

    // print LCS score matrix allong with sequences

    printf("    ");
    printf("%5c   ", ' ');

    for (j = 0; j < sizeA; j++)
        printf("%5c   ", seqA[j]);
    printf("\n");
    for (i = 0; i < sizeB + 1; i++)
    {
        if (i == 0)
            printf("    ");
        else
            printf("%c   ", seqB[i - 1]);
        for (j = 0; j < sizeA + 1; j++)
        {
            printf("%5d   ", scoreMatrix[i][j]);
        }
        printf("\n");
    }
    printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB)
{
    int i;
    for (i = 0; i < (sizeB + 1); i++)
        free(scoreMatrix[i]);
    free(scoreMatrix);
}

// mpicc -o lcs-mpi lcs-mpi.c
int main(int argc, char **argv)
{
    int my_rank, n_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    // sequence pointers for both sequences
    char *seqA, *seqB;

    // sizes of both sequences
    int global_sizes[2];

    // só o rank 0 le os dados e depois transmite
    if (my_rank == 0)
    {
        seqA = read_seq(argv[1]);
        seqB = read_seq(argv[2]);

        global_sizes[0] = strlen(seqA);
        global_sizes[1] = strlen(seqB);
    }

    MPI_Bcast(global_sizes, 2, MPI_INT, 0, MPI_COMM_WORLD);
    int sizeA = global_sizes[0];
    int sizeB = global_sizes[1];

    if (my_rank != 0)
    {
        seqA = (char *)malloc((sizeA + 1) * sizeof(char));
        seqB = (char *)malloc((sizeB + 1) * sizeof(char));
    }

    MPI_Bcast(seqA, sizeA + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(seqB, sizeB + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    // allocate LCS score matrix
    // agora é alocado por cada nó
    // mtype **scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

    // initialize LCS score matrix
    // initScoreMatrix(scoreMatrix, sizeA, sizeB);

    // fill up the rest of the matrix and return final score (element locate at the last line and collumn)
    // int LCS(MPI_Comm comm_world, int my_rank, int n_procs, char *seqA, int sizeA, char *seqB, int sizeB)
    double start_time = MPI_Wtime();
    mtype score = LCS(MPI_COMM_WORLD, my_rank, n_procs, seqA, sizeA, seqB, sizeB);
    double end_time = MPI_Wtime();

    // print score
    // double tempoFinal = fim - inicio;
    if (my_rank == 0)
    {
        // printf("\nScore: %d\n", score);
        printf("%0.8f\n", end_time - start_time);
    }
    free(seqA);
    free(seqB);
    MPI_Finalize();
    return EXIT_SUCCESS;
}