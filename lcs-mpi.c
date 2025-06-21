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
    // file pointer
    FILE *fseq = NULL;
    // sequence size
    long size = 0;
    // sequence pointer
    char *seq = NULL;
    // sequence index
    int i = 0;

    // open file
    fseq = fopen(fname, "rt");
    if (fseq == NULL)
    {
        printf("Error reading file %s\n", fname);
        exit(1);
    }

    // find out sequence size to allocate memory afterwards
    fseek(fseq, 0L, SEEK_END);
    size = ftell(fseq);
    rewind(fseq);

    // allocate memory (sequence)
    seq = (char *)calloc(size + 1, sizeof(char));
    if (seq == NULL)
    {
        printf("Erro allocating memory for sequence %s.\n", fname);
        exit(1);
    }

    // read sequence from file
    while (!feof(fseq))
    {
        seq[i] = fgetc(fseq);
        if ((seq[i] != '\n') && (seq[i] != EOF))
            i++;
    }
    // insert string terminator
    seq[i] = '\0';

    // close file
    fclose(fseq);

    // return sequence pointer
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

int LCS(MPI_Comm comm_world, int my_rank, int n_procs, char *seqA, int sizeA, char *seqB, int sizeB)
{

    // setup cart
    MPI_Comm comm_cart;
    int ndims = 2, dims[2] = {2, 3}, periods[2] = {0, 0}, reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);

    if (comm_cart == MPI_COMM_NULL)
        return 0; // TODO verificar esse retorno

    int my_cart_rank, my_coords[2];
    MPI_Comm_rank(comm_cart, &my_cart_rank);
    MPI_Cart_coords(comm_cart, my_cart_rank, ndims, my_coords);

    // Descobre vizinhos
    int rank_up, rank_down, rank_left, rank_right;
    MPI_Cart_shift(comm_cart, 0, 1, &rank_up, &rank_down);
    MPI_Cart_shift(comm_cart, 1, 1, &rank_left, &rank_right);

    // aloca dados locais
    // assim nao precisa alocar uma matriz tao grande, cada processo aloca um pedaço
    int start_col = (my_coords[1] * (sizeA + 1)) / dims[1];
    int end_col = ((my_coords[1] + 1) * (sizeA + 1)) / dims[1] - 1;
    int local_cols = end_col - start_col + 1;

    int start_row = (my_coords[0] * (sizeB + 1)) / dims[0];
    int end_row = ((my_coords[0] + 1) * (sizeB + 1)) / dims[0] - 1;
    int local_rows = end_row - start_row + 1;

    if (my_rank == 0)
    { // debug
        printf("Tamanho da matriz global: %d x %d\n", sizeB + 1, sizeA + 1);
    }
    printf("Tamanho do bloco local por processo: %d x %d\n", local_rows, local_cols);

    mtype **matriz_local = alocar_matriz_local(local_rows, local_cols);

    mtype top_buffer[local_cols];
    mtype left_buffer[local_rows];

    int total_diagonals = dims[0] + dims[1] - 1;

    // processos da borda nao esperam dados
    if (rank_up != MPI_PROC_NULL)
    {
        MPI_Recv(&matriz_local[0][1], local_cols, MPI_UNSIGNED_SHORT, rank_up, 0, comm_cart, MPI_STATUS_IGNORE);
    }

    if (rank_left != MPI_PROC_NULL)
    {
        // dataType novo pra receber a coluna
        MPI_Datatype col_type;
        MPI_Type_vector(local_rows, 1, local_cols + 2, MPI_UNSIGNED_SHORT, &col_type);
        MPI_Type_commit(&col_type);
        MPI_Recv(&matriz_local[1][0], 1, col_type, rank_left, 0, comm_cart, MPI_STATUS_IGNORE);
        MPI_Type_free(&col_type);
    }

    // calcula bloco local
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
    // envia dados, se nao for a borda
    if (rank_down != MPI_PROC_NULL)
    {
        MPI_Send(&matriz_local[local_rows][1], local_cols, MPI_UNSIGNED_SHORT, rank_down, 0, comm_cart);
    }

    if (rank_right != MPI_PROC_NULL)
    {
        // Prepara o buffer com a última coluna calculada
        MPI_Datatype col_type;
        MPI_Type_vector(local_rows, 1, local_cols + 2, MPI_UNSIGNED_SHORT, &col_type);
        MPI_Type_commit(&col_type);
        MPI_Send(&matriz_local[1][local_cols], 1, col_type, rank_right, 0, comm_cart);
        MPI_Type_free(&col_type);
    }

    mtype final_score = 0;

    // O processo no canto inferior direito tem o resultado final
    if (rank_down == MPI_PROC_NULL && rank_right == MPI_PROC_NULL)
    {
        final_score = matriz_local[local_rows][local_cols];
        // Se este processo não for o rank 0, envie o resultado para ele
        if (my_rank != 0)
        {
            MPI_Send(&final_score, 1, MPI_UNSIGNED_SHORT, 0, 0, comm_world);
        }
    }

    // O Rank 0 (global) recebe o resultado final se ele não for o último processo.
    if (my_rank == 0 && (rank_down != MPI_PROC_NULL || rank_right != MPI_PROC_NULL))
    {
        int rank_do_ultimo_proc;
        int coords_do_ultimo[2] = {dims[0] - 1, dims[1] - 1};
        MPI_Cart_rank(comm_cart, coords_do_ultimo, &rank_do_ultimo_proc);
        MPI_Recv(&final_score, 1, MPI_UNSIGNED_SHORT, rank_do_ultimo_proc, 0, comm_world, MPI_STATUS_IGNORE);
    }

    // limpeza
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
        printf("\nScore: %d\n", score);
        printf("Tempo de execução: %f segundos\n", end_time - start_time);
    }
    free(seqA);
    free(seqB);
    MPI_Finalize();
    return 0;
}