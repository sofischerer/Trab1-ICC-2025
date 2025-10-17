#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

void multiplicar_matrizes(real_t *A, real_t *B, real_t *C, int m, int n, int p) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            C[i * p + j] = 0.0f;
            for (int k = 0; k < n; ++k) {
                C[i * p + j] += A[i * n + k] * B[k * p + j];
            }
        }
    }
}

void print_matriz(real_t* A, int lin, int col){
    for(int i = 0; i < lin; i++){
        for(int j = 0; j < col; j++){
            printf("%13.8g ", A[(i*col)+j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}


/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t *A, real_t *B)
{
  int metade = (k/2);
  for (int i = 0; i<n; i++){
    for (int j=max(i-metade, 0); j<=min(n-1, i+metade); j++){
      A[(n*i)+j] = generateRandomA(i, j, k);
    }
  }
  for (int i = 0; i<n; i++){
      B[i] = generateRandomB(k);
  }
}

void calcular_AtA(real_t* A, real_t* AtA, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            real_t soma = 0.0;
            for (int k = 0; k < n; k++) {
                soma += A[i * n + k] * A[j * n + k];
            }
            AtA[i * n + j] = soma;
        }
    }
}

void calcular_Atb(real_t* A, real_t* b, real_t* Atb, int n) {
    for (int i = 0; i < n; i++) {
        real_t soma = 0.0;
        for (int k = 0; k < n; k++) {
            soma += A[n*k + i] * b[k];
        }
        Atb[i] = soma;
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, 
			  real_t *ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();

  calcular_AtA(A, ASP, n);
  calcular_Atb(A, b, bsp, n);
  
  *tempo = timestamp() - *tempo;
 
}


void geraDLU (real_t *A, int n, int k,
	      real_t *D, real_t *L, real_t *U, rtime_t *tempo)
{
  *tempo = timestamp();

  int metade = k/2;
  for(int i=0; i<n; i++){
    for(int j=max(0, i-metade); j<i; j++){
      L[(i*n)+j] = A[(i*n)+j];
    }

    D[(i*n)+i] = A[(i*n)+i];

    for(int j=i+1; j<min(i+metade+1, n); j++){
      U[(i*n)+j] = A[(i*n)+j];
    }    
  }

  *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k,
		 real_t *M, rtime_t *tempo)
{
  *tempo = timestamp();
  if(w==-1){
    for(int i=0; i<n;i++){
      M[(i*n)+i] = 1;
    }
  }
  else if(w==0.0){
    for(int i=0; i<n; i++){
      M[(i*n)+i] = 1.0/D[(i*n)+i];
    }
  }else if(w==1.0){
    //IMPLEMENTAR ERRO
  }else if(w>1.0 && w<=2.0){
    //IMPLEMENTAR ERRO
  }else{
    //IMPLEMENTAR ERRO
  }


  *tempo = timestamp() - *tempo;
}


void calcResiduoSL (real_t *A, real_t *b, real_t *X, real_t* r,
		      int n, int k, rtime_t *tempo)
{
  *tempo = timestamp();

  real_t tmp[n];
  multiplicar_matrizes(A, X, tmp, n, n, 1);
  print_matriz(tmp, 1, n);
  for (int i=0; i<n; i++){
    r[i] = b[i] - tmp[i];
  }

  *tempo = timestamp() - *tempo;
}


