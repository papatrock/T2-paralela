#!/bin/bash

NUM_EXECUCOES=20
ARQUIVO_A="entradas/40k_A.in"
ARQUIVO_B="entradas/40k_B.in"
ARQUIVO_SAIDA_SEQ="saidas/tempos_sequencial_40k.txt"

# rm -f $ARQUIVO_SAIDA_SEQ
# touch $ARQUIVO_SAIDA_SEQ

echo "Iniciando $NUM_EXECUCOES execuções para o código sequencial (UMA DE CADA VEZ)..."

for i in $(seq 1 $NUM_EXECUCOES)
do
    echo "--> Submetendo e AGUARDANDO a execução $i de $NUM_EXECUCOES..."

    sbatch --wait prog.sbatch # $ARQUIVO_A $ARQUIVO_B $ARQUIVO_SAIDA_SEQ
done

echo "Todos os $NUM_EXECUCOES jobs foram executados em sequência."
# echo "Resultados acumulados em '$ARQUIVO_SAIDA_SEQ'."