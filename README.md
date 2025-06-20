# Resumo das Principais Funções MPI em C (Open MPI)

Este é um guia de referência rápida para as funções mais comuns usadas na programação paralela com MPI em C.

## 1. Gerenciamento do Ambiente MPI

Funções para inicializar e finalizar o ambiente MPI. Devem ser as primeiras e últimas chamadas MPI em seu programa.

| Função | Descrição | Sintaxe em C |
| :--- | :--- | :--- |
| `MPI_Init` | Inicializa o ambiente de execução do MPI. | `int MPI_Init(int *argc, char ***argv);` |
| `MPI_Finalize` | Finaliza o ambiente MPI. Nenhuma outra função MPI pode ser chamada após esta. | `int MPI_Finalize(void);` |

## 2. Consulta de Informações do Comunicador

Funções para obter o tamanho do grupo de processos e o identificador (rank) de um processo específico.

| Função | Descrição | Sintaxe em C |
| :--- | :--- | :--- |
| `MPI_Comm_size` | Obtém o número total de processos em um comunicador. | `int MPI_Comm_size(MPI_Comm comm, int *size);` |
| `MPI_Comm_rank` | Obtém o rank (ID de 0 a `size-1`) do processo atual. | `int MPI_Comm_rank(MPI_Comm comm, int *rank);` |

## 3. Comunicação Ponto a Ponto

Usadas para enviar e receber mensagens entre dois processos específicos.

| Função | Descrição | Sintaxe em C |
| :--- | :--- | :--- |
| `MPI_Send` | **Envio bloqueante:** O processo remetente espera até que o destinatário inicie o recebimento. | `int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);` |
| `MPI_Recv` | **Recebimento bloqueante:** O processo fica bloqueado até que a mensagem seja completamente recebida. | `int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);` |
| `MPI_Isend` | **Envio não bloqueante:** Retorna imediatamente, permitindo que o processo continue executando. | `int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);` |
| `MPI_Irecv` | **Recebimento não bloqueante:** Inicia o recebimento e retorna imediatamente. | `int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);` |
| `MPI_Wait` | Espera pela conclusão de uma operação não bloqueante (`MPI_Isend`/`MPI_Irecv`). | `int MPI_Wait(MPI_Request *request, MPI_Status *status);` |

## 4. Comunicação Coletiva

Envolvem a comunicação entre **todos** os processos de um comunicador.

| Função | Descrição | Sintaxe em C |
| :--- | :--- | :--- |
| `MPI_Bcast` | **Broadcast:** Envia dados de um processo raiz para todos os outros. | `int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);` |
| `MPI_Scatter` | **Espalhar:** Distribui um array da raiz em partes iguais para todos os processos. | `int MPI_Scatter(const void *sendbuf, int sendcount, ..., void *recvbuf, int recvcount, ..., int root, MPI_Comm comm);` |
| `MPI_Gather` | **Juntar:** Coleta dados de todos os processos e os agrupa em um processo raiz. | `int MPI_Gather(const void *sendbuf, int sendcount, ..., void *recvbuf, int recvcount, ..., int root, MPI_Comm comm);` |
| `MPI_Allgather`| **Juntar para Todos:** Coleta dados de todos e distribui o resultado para todos. | `int MPI_Allgather(const void *sendbuf, int sendcount, ..., void *recvbuf, int recvcount, ..., MPI_Comm comm);` |
| `MPI_Reduce` | **Redução:** Realiza uma operação (soma, max, min) nos dados de todos e armazena o resultado na raiz. | `int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);` |
| `MPI_Allreduce`| **Redução para Todos:** Similar ao `MPI_Reduce`, mas o resultado é distribuído para todos. | `int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);` |
| `MPI_Barrier` | **Barreira:** Sincroniza todos os processos. Nenhum avança até que todos cheguem neste ponto. | `int MPI_Barrier(MPI_Comm comm);` |

---


| Função | Descrição | Sintaxe em C |
| :--- | :--- | :--- |
| **`MPI_Cart_create`** | A função mais importante. Cria um **novo comunicador** com uma topologia cartesiana (grade). Você define o número de dimensões e o tamanho de cada uma. | `int MPI_Cart_create(MPI_Comm old_comm, int ndims, const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart);` |
| **`MPI_Cart_coords`** | Traduz o **rank** de um processo no comunicador cartesiano para suas **coordenadas** lógicas na grade (ex: rank 5 -> coordenadas `[1, 2]`). | `int MPI_Cart_coords(MPI_Comm comm_cart, int rank, int maxdims, int coords[]);` |
| **`MPI_Cart_rank`** | O inverso de `MPI_Cart_coords`. Traduz as **coordenadas** da grade para o **rank** do processo correspondente (ex: coordenadas `[1, 2]` -> rank 5). | `int MPI_Cart_rank(MPI_Comm comm_cart, const int coords[], int *rank);` |
| **`MPI_Cart_shift`** | **Essencial para comunicação.** Encontra os ranks dos processos vizinhos ao longo de uma dimensão específica (ex: "quem é meu vizinho para cima e para baixo?"). | `int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest);` |
| **`MPI_Cart_sub`** | (Avançado) Particiona a grade para criar sub-comunicadores. Muito útil para fazer uma operação coletiva (como um `MPI_Bcast`) apenas em uma linha ou coluna da grade. | `int MPI_Cart_sub(MPI_Comm comm_cart, const int remain_dims[], MPI_Comm *new_comm);` |
| **`MPI_Comm_free`** | Libera um comunicador que não será mais usado (como o `comm_cart` criado por `MPI_Cart_create`). É uma boa prática de gerenciamento de memória. | `int MPI_Comm_free(MPI_Comm *comm);` |

---

```c
#include <mpi.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // 1. DEFINIR A GRADE (Ex: grade 2D para 4 processos, 2x2)
    int ndims = 2;
    int dims[2] = {2, 2}; // 2 processos na dimensão 0 (linhas), 2 na dimensão 1 (colunas)
    int periods[2] = {0, 0}; // Não-periódico (sem "wrap-around")
    int reorder = 1; // Permite ao MPI reordenar os ranks para otimização

    MPI_Comm comm_cart; // Novo comunicador cartesiano
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);


    // 2. OBTER INFORMAÇÕES SOBRE O PROCESSO ATUAL NA GRADE
    int my_cart_rank;
    int my_coords[2]; // Coordenadas [linha, coluna]
    MPI_Comm_rank(comm_cart, &my_cart_rank); // Pega o rank DENTRO do novo comunicador
    MPI_Cart_coords(comm_cart, my_cart_rank, ndims, my_coords);


    // 3. ENCONTRAR OS VIZINHOS
    int rank_up, rank_down, rank_left, rank_right;

    // Vizinhos na dimensão 0 (para cima e para baixo)
    MPI_Cart_shift(comm_cart, 0, 1, &rank_up, &rank_down);

    // Vizinhos na dimensão 1 (para esquerda e para direita)
    MPI_Cart_shift(comm_cart, 1, 1, &rank_left, &rank_right);
    // Nota: Se um vizinho não existe, a função retorna MPI_PROC_NULL.


    // 4. EXECUTAR A LÓGICA DE COMUNICAÇÃO
    // if (rank_up != MPI_PROC_NULL) MPI_Recv(..., rank_up, ...);
    // ... calcular o bloco local ...
    // if (rank_down != MPI_PROC_NULL) MPI_Send(..., rank_down, ...);


    // 5. LIMPAR
    // Libera o comunicador que criamos quando não precisarmos mais dele
    MPI_Comm_free(&comm_cart);

    MPI_Finalize();
    return 0;
}
```