#include "sislin.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

void alocar_elementos(real_t** A, real_t** B, real_t** X, real_t** ASP, real_t** bsp, int n){
    *A = calloc(n*n, sizeof(real_t));
    *B = calloc(n, sizeof(real_t));
    *X = calloc(n, sizeof(real_t));
    *ASP = calloc(n*n, sizeof(real_t));
    *bsp = calloc(n, sizeof(real_t));
    return;
}

void test_prints(real_t* A, real_t* b, real_t* X, real_t* ASP, real_t* bsp, int n){

    printf("\n");
    print_matriz(A, n, n);
    print_matriz(b, 1, n);
    print_matriz(X, 1, n);
    print_matriz(ASP, n, n);
    print_matriz(bsp, 1, n);
    return;
}

int main(){

    srandom(20252);
    
    real_t* A;
    real_t* b;
    real_t* X;
    real_t* ASP;
    real_t* bsp;
    real_t tempo = 0;
    int k, n, w, maxiter;
    real_t erro;

    scanf("%d\n%d\n%d\n%d\n%lf", &n, &k, &w, &maxiter, &erro);
    alocar_elementos(&A, &b, &X, &ASP, &bsp, n);
    criaKDiagonal(n, k, A, b);
    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tempo);
    test_prints(A, b, X, ASP, bsp, n);
}