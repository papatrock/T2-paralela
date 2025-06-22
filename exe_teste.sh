#!/bin/bash

NUM_EXECUCOES=20
#SBATCH=prog.sbatch
SBATCH=sequencial.sbatch

# rm -f $ARQUIVO_SAIDA_SEQ
# touch $ARQUIVO_SAIDA_SEQ

echo "Iniciando $NUM_EXECUCOES execuções"

for i in $(seq 1 $NUM_EXECUCOES)
do
    echo "--> Submetendo e AGUARDANDO a execução $i de $NUM_EXECUCOES..."

    sbatch --wait $SBATCH
done

echo "Todos os $NUM_EXECUCOES jobs foram executados em sequência."