#!/bin/sh

# Nome do job
#SBATCH --job-name=sequencial_pdol

# Numero total de processos MPI
#SBATCH --ntasks=1

# Nome do arquivo de saida. %j sera o ID do seu job.
#SBATCH --output=/dev/null

./lcs-seq entradas/40k_A.in entradas/40k_B.in >> saidas/tempo_sequencial_40k.out