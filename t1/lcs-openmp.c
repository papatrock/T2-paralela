#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#define NUM_THREADS 16

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

mtype **allocateScoreMatrix(int sizeA, int sizeB)
{
  int i;
  // Allocate memory for LCS score matrix
  mtype **scoreMatrix = (mtype **)malloc((sizeB + 1) * sizeof(mtype *));
  for (i = 0; i < (sizeB + 1); i++)
    scoreMatrix[i] = (mtype *)malloc((sizeA + 1) * sizeof(mtype));
  return scoreMatrix;
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

int calculateBlockSize(int sizeA, int sizeB, int numThreads)
{
  int totalElements = sizeA * sizeB;
  int elementsPerThread = totalElements / numThreads;
  int blockSize = (int)sqrt(elementsPerThread);

  if (blockSize > sizeA || blockSize > sizeB)
  {
    blockSize = (sizeA < sizeB) ? sizeA : sizeB;
  }

  // if (blockSize < 16)
  //   blockSize = 16;

  return (blockSize > 0) ? blockSize : 1; // sem tamanho minimo

  // return blockSize;
}

int LCS(mtype **scoreMatrix, int sizeA, int sizeB, char *seqA, char *seqB)
{
  int bi, bj, i, j;

  int numThreads = omp_get_max_threads();
  int blockSize = calculateBlockSize(sizeA, sizeB, numThreads);

  int numBlockRows = (sizeB + blockSize - 1) / blockSize; // (7+4-1) / 4  = 2
  int numBlockCols = (sizeA + blockSize - 1) / blockSize; // (10 + 4 - 1) / 4 = 3
  int totalDiagonals = numBlockRows + numBlockCols - 1;   // 4 + 10 - 1

  for (int d = 0; d < totalDiagonals; ++d) // percorre as diagonais
  {
#pragma omp parallel for private(bi, bj, i, j) schedule(dynamic)
    for (bi = 0; bi <= d; ++bi) // uma thread pega o 0, outra a 1 .... até d
    {

      bj = d - bi;

      if (bi < numBlockRows && bj < numBlockCols)
      {
        int startRow = bi * blockSize + 1;
        int endRow = (bi + 1) * blockSize;
        int startCol = bj * blockSize + 1;
        int endCol = (bj + 1) * blockSize;

        if (endRow > sizeB)
          endRow = sizeB;
        if (endCol > sizeA)
          endCol = sizeA;

#ifdef DEBUGMATRIX
        printf("thread:%d processando o bloco: bi=%d bj=%d, linhas %d–%d, colunas %d–%d\n", omp_get_thread_num(), bi, bj, startRow, endRow, startCol, endCol);
#endif

        for (i = startRow; i <= endRow; ++i)
        {
          for (j = startCol; j <= endCol; ++j)
          {
#ifdef DEBUGMATRIX
            printf("%c == %c ?\n", seqA[j - 1], seqB[i - 1]);
#endif
            if (seqA[j - 1] == seqB[i - 1])
              scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
            else
              scoreMatrix[i][j] = (scoreMatrix[i - 1][j] > scoreMatrix[i][j - 1])
                                      ? scoreMatrix[i - 1][j]
                                      : scoreMatrix[i][j - 1];
          }
        }
      }
    }
  }

  return scoreMatrix[sizeB][sizeA];
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
#pragma omp parallel for
  for (i = 0; i < (sizeB + 1); i++)
    free(scoreMatrix[i]);
  free(scoreMatrix);
}

// gcc -O3 -DDEBUGMATRIX -fopenmp lcs-paralelizado.c -o lcs-para -lm
int main(int argc, char **argv)
{
  // sequence pointers for both sequences
  char *seqA, *seqB;

  // sizes of both sequences
  int sizeA, sizeB;

  // read both sequences
  int num_threads = NUM_THREADS;
  char *fileA = "fileA.in";
  char *fileB = "fileB.in";

  if (argc > 1)
    fileA = argv[1];
  if (argc > 2)
    fileB = argv[2];
  if (argc > 3)
    num_threads = atoi(argv[3]);

  omp_set_num_threads(num_threads);

  seqA = read_seq(fileA);
  seqB = read_seq(fileB);

  // find out sizes
  sizeA = strlen(seqA);
  sizeB = strlen(seqB);
  // allocate LCS score matrix
  mtype **scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

  // initialize LCS score matrix
  initScoreMatrix(scoreMatrix, sizeA, sizeB);

  // fill up the rest of the matrix and return final score (element locate at the last line and collumn)
  double start = omp_get_wtime();
  mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);
  double end = omp_get_wtime();

  /* if you wish to see the entire score matrix,
   for debug purposes, define DEBUGMATRIX. */
#ifdef DEBUGMATRIX
  printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
#endif

  // print score
  double time = end - start;
  printf("\nScore: %d tempo:%0.8f\n", score, time);
  printf("%0.8f", time);
  // free score matrix
  freeScoreMatrix(scoreMatrix, sizeB);

  return EXIT_SUCCESS;
}
