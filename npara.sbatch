#!/bin/sh

# Nome do job
#SBATCH --job-name=mpi_dopol_test

# Numero de nós
#SBATCH --nodes=1

# Numero total de processos MPI
#SBATCH --ntasks=1

#SBATCH --ntasks-per-node=1

# Nome do arquivo de saida. %j sera o ID do seu job.
#SBATCH --output=resultado_%j.out

echo "Job iniciado em: $(date)"
echo "Nos alocados: $SLURM_JOB_NODELIST"
echo "Total de tasks: $SLURM_NTASKS"

# O mpirun usa a variavel $SLURM_NTASKS que o SLURM cria automaticamente.
# O --map-by garante a distribuicao de 1 processo por no.
./lcs entradas/40k_A.in entradas/40k_B.in

echo "Job finalizado em: $(date)"