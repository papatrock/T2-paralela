// --- Início do Programa MPI ---
MPI_Init(...);

// 1. Criar uma grade 2D de processos
MPI_Cart_create(...);
// Descobrir as coordenadas (meu_i, meu_j) do processo atual na grade
MPI_Cart_coords(...);
// Descobrir quem são os vizinhos de cima, baixo, esquerda, direita
MPI_Cart_shift(...);

// 2. Alocar memória APENAS para o meu bloco local
bloco_local = alocar_meu_bloco(...);
// Inicializar as bordas do meu bloco se eu estiver na borda da matriz global

// --- Laço Principal da Wavefront ---
// (A lógica exata do loop pode variar, mas a ideia é esta)

// 3. Receber dados dos vizinhos de cima e da esquerda (se eles existirem)
if (vizinho_de_cima != MPI_PROC_NULL) {
    MPI_Recv(borda_de_cima, ..., vizinho_de_cima, ...);
}
if (vizinho_da_esquerda != MPI_PROC_NULL) {
    MPI_Recv(borda_da_esquerda, ..., vizinho_da_esquerda, ...);
}

// 4. Calcular o meu bloco local usando os dados recebidos
calcular_bloco(bloco_local, borda_de_cima, borda_da_esquerda);

// 5. Enviar meus resultados para os vizinhos de baixo e da direita (se eles existirem)
if (vizinho_de_baixo != MPI_PROC_NULL) {
    MPI_Send(minha_borda_de_baixo, ..., vizinho_de_baixo, ...);
}
if (vizinho_da_direita != MPI_PROC_NULL) {
    MPI_Send(minha_borda_da_direita, ..., vizinho_da_direita, ...);
}

// 6. O resultado final estará no processo do último bloco (canto inferior direito)
// Ele pode enviar o resultado para o processo 0 no final.

MPI_Finalize();