#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define GRID_SIZE 10 // Tamanho da grade global (sem contar as bordas fixas)
#define LOCAL_SIZE (GRID_SIZE / 2) // Tamanho da sub-grade de cada processo (10/2=5)
#define TIMESTEPS 100 // Quantos passos de tempo simular

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int world_rank, procs;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    // Este exemplo foi feito especificamente para 4 processos
    if (procs != 4) {
        if (world_rank == 0) {
            printf("Por favor, rode este exemplo com 4 processos (-np 4).\n");
        }
        MPI_Finalize();
        return 1;
    }

    // 1. CRIAÇÃO DA TOPOLOGIA CARTESIANA 2x2
    MPI_Comm comm_cart;
    int dims[2] = {2, 2};
    int periods[2] = {0, 0}; // 0 = não-periódico (bordas fixas)
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_cart);

    int my_rank, my_coords[2];
    MPI_Comm_rank(comm_cart, &my_rank);
    MPI_Cart_coords(comm_cart, my_rank, 2, my_coords);

    // 2. ENCONTRAR VIZINHOS
    int rank_up, rank_down, rank_left, rank_right;
    MPI_Cart_shift(comm_cart, 0, 1, &rank_up, &rank_down);
    MPI_Cart_shift(comm_cart, 1, 1, &rank_left, &rank_right);

    // 3. ALOCAÇÃO DAS GRADES (LOCAL + CÉLULAS FANTASMAS)
    // Alocamos (N+2) x (N+2) para ter a moldura de 1 célula
    double local_grid[LOCAL_SIZE + 2][LOCAL_SIZE + 2] = {0.0};
    double new_local_grid[LOCAL_SIZE + 2][LOCAL_SIZE + 2] = {0.0};

    // 4. INICIALIZAÇÃO (Condições de Contorno)
    // Processos na borda de cima da grade global (linha 0)
    if (my_coords[0] == 0) {
        for (int i = 0; i < LOCAL_SIZE + 2; i++) local_grid[0][i] = 100.0;
    }
    // Processos na borda de baixo da grade global (linha 1)
    if (my_coords[0] == 1) {
        for (int i = 0; i < LOCAL_SIZE + 2; i++) local_grid[LOCAL_SIZE + 1][i] = 100.0;
    }
    // (Bordas da esquerda e direita são 0.0, já inicializado)

    // Tipos de dados MPI para enviar/receber colunas (um vetor não-contíguo)
    MPI_Datatype col_type;
    MPI_Type_vector(LOCAL_SIZE, 1, LOCAL_SIZE + 2, MPI_DOUBLE, &col_type);
    MPI_Type_commit(&col_type);


    // 5. LAÇO PRINCIPAL DA SIMULAÇÃO
    for (int t = 0; t < TIMESTEPS; t++) {

        // --- TROCA DE BORDAS COM VIZINHOS (A DANÇA DA COMUNICAÇÃO) ---
        // MPI_Sendrecv é ótimo aqui pois envia e recebe ao mesmo tempo, evitando deadlocks

        // Envia para cima, recebe de baixo
        MPI_Sendrecv(&local_grid[1][1], LOCAL_SIZE, MPI_DOUBLE, rank_up, 0,
                     &local_grid[LOCAL_SIZE + 1][1], LOCAL_SIZE, MPI_DOUBLE, rank_down, 0,
                     comm_cart, MPI_STATUS_IGNORE);

        // Envia para baixo, recebe de cima
        MPI_Sendrecv(&local_grid[LOCAL_SIZE][1], LOCAL_SIZE, MPI_DOUBLE, rank_down, 0,
                     &local_grid[0][1], LOCAL_SIZE, MPI_DOUBLE, rank_up, 0,
                     comm_cart, MPI_STATUS_IGNORE);

        // Envia para esquerda, recebe da direita (usando o tipo de dado de coluna)
        MPI_Sendrecv(&local_grid[1][1], 1, col_type, rank_left, 0,
                     &local_grid[1][LOCAL_SIZE + 1], 1, col_type, rank_right, 0,
                     comm_cart, MPI_STATUS_IGNORE);

        // Envia para direita, recebe da esquerda
        MPI_Sendrecv(&local_grid[1][LOCAL_SIZE], 1, col_type, rank_right, 0,
                     &local_grid[1][0], 1, col_type, rank_left, 0,
                     comm_cart, MPI_STATUS_IGNORE);


        // --- COMPUTAÇÃO LOCAL ---
        // Agora que as células fantasmas estão preenchidas, calculamos o interior
        for (int i = 1; i <= LOCAL_SIZE; i++) {
            for (int j = 1; j <= LOCAL_SIZE; j++) {
                new_local_grid[i][j] = (local_grid[i - 1][j] + local_grid[i + 1][j] +
                                        local_grid[i][j - 1] + local_grid[i][j + 1]) / 4.0;
            }
        }

        // --- ATUALIZAÇÃO DA GRADE PARA A PRÓXIMA ITERAÇÃO ---
        for (int i = 1; i <= LOCAL_SIZE; i++) {
            for (int j = 1; j <= LOCAL_SIZE; j++) {
                local_grid[i][j] = new_local_grid[i][j];
            }
        }

        // Sincroniza todos os processos antes do próximo passo de tempo
        MPI_Barrier(comm_cart);
    }

    // Imprime um valor do centro da grade local para ver o resultado
    if (my_rank == 0) {
        printf("Simulação finalizada após %d passos.\n", TIMESTEPS);
        printf("Valor aproximado no centro da minha grade (rank 0): %f\n", local_grid[LOCAL_SIZE/2][LOCAL_SIZE/2]);
    }

    // 6. LIMPEZA
    MPI_Type_free(&col_type);
    MPI_Comm_free(&comm_cart);
    MPI_Finalize();

    return 0;
}