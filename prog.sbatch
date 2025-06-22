#!/bin/sh

# Nome do job
#SBATCH --job-name=mpi_dopol_test

# Numero de nós
#SBATCH --nodes=2

# Numero total de processos MPI
#SBATCH --ntasks=6

#SBATCH --ntasks-per-node=3

# Nome do arquivo de saida. %j sera o ID do seu job.
#SBATCH --output=saida_%j.out


#echo "Job iniciado em: $(date)"
#echo "Nos alocados: $SLURM_JOB_NODELIST"
#echo "Total de tasks: $SLURM_NTASKS"

# O mpirun usa a variavel $SLURM_NTASKS que o SLURM cria automaticamente.
# O --map-by garante a distribuicao de 1 processo por no.
mpirun -np $SLURM_NTASKS ./lcs-mpi entradas/5k_A.in entradas/5k_B.in

#echo "Job finalizado em: $(date)"