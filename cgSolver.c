#include "sislin.h"
#include "utils.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>


void print_matriz(real_t** A, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%13.8g ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}

void print_vetor(real_t* B, int n){
    for (int i = 0; i<n; i++){
        printf("%13.8g ", B[i]);
    }
    printf("\n");
}

real_t* transfere_diagonais(real_t** matriz, int k, int n){
    int d = (k-1)/2;
    int total_elementos = (n*k)-d*(d+1);
    int metade = (k/2);
    int pos = 0;
    real_t* vetor = calloc(total_elementos, sizeof(real_t));
    for (int i=0; i<n; i++){
        for(int j=i-metade;j<i+metade+1; j++){
            if(j>=0 && j<n){
                vetor[pos] = matriz[i][j];
                pos++;
            }
        }
    }
    return vetor;
}

void alocar_elementos(){

}

int main(){
    srandom(20252);
    
    int k = 5;
    int n = 10;
    real_t** matriz = calloc(k*n, sizeof(real_t*));
    real_t* B = calloc(n, sizeof(real_t));
    real_t* X = calloc(n, sizeof(real_t));
    criaKDiagonal(n, k, matriz, B);
    real_t* A = transfere_diagonais(matriz, k, n);
    
    print_matriz(matriz, n);
    // print_vetor(A, (n*k)-((k-1)/2)*(((k-1)/2)+1));
    print_vetor(B, n);
}