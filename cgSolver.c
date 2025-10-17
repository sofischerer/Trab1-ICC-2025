#include "sislin.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

void test_prints(real_t* A, real_t* b, real_t* X, real_t* ASP, real_t* bsp, real_t* D, real_t* L, real_t* U, real_t* r, int n){

    printf("\n");
    printf("---------------------------------------------\nA:\n");
    print_matriz(A, n, n);
    printf("---------------------------------------------\nb:\n");
    print_matriz(b, 1, n);
    printf("---------------------------------------------\nX:\n");
    print_matriz(X, 1, n);
    printf("---------------------------------------------\nASP:\n");
    print_matriz(ASP, n, n);
    printf("---------------------------------------------\nbsp:\n");
    print_matriz(bsp, 1, n);
    printf("---------------------------------------------\nD:\n");
    print_matriz(D, n, n);
    printf("---------------------------------------------\nL:\n");
    print_matriz(L, n, n);
    printf("---------------------------------------------\nU:\n");
    print_matriz(U, n, n);
    printf("---------------------------------------------\nr:\n");
    print_matriz(r, 1, n);

    return;
}

real_t calc_err(real_t* x, real_t* x_old, int n){
    real_t max = fabs(x[0] - x_old[0]);
    real_t tmp;
    for(int i=1; i<n; i++){
        tmp = fabs(x[i] - x_old[i]);
        if(tmp > max) max = tmp;
    }
    return max;
}

real_t escalar(real_t* A, int n){
    real_t sum = 0.0;
    for(int i=0; i<n; i++){
        sum += A[i]*A[i];
    }
    return sum;
}

void copiar_vetor(real_t* origem, real_t* destino, int n){
    for(int i=0; i<n; i++){
        destino[i] = origem[i];
    }
}


int main(){

    srandom(20252);    
    real_t tempo = 0;
    int k, n, w, maxiter;
    real_t max_err;

    scanf("%d %d %d %d %lf", &n, &k, &w, &maxiter, &max_err);

    //Feio mas hey as long as it works
    real_t *A = calloc(n*n, sizeof(real_t));
    real_t *b = calloc(n, sizeof(real_t));
    real_t *X = calloc(n, sizeof(real_t));
    real_t *ASP = calloc(n*n, sizeof(real_t));
    real_t *bsp = calloc(n, sizeof(real_t));
    real_t *D = calloc(n*n, sizeof(real_t));
    real_t *L = calloc(n*n, sizeof(real_t));
    real_t *U = calloc(n*n, sizeof(real_t));
    real_t *r = calloc(n, sizeof(real_t));
    real_t *M = calloc(n*n, sizeof(real_t));
    real_t *x_old = calloc(n, sizeof(real_t));
    real_t *r_old = calloc(n, sizeof(real_t));
    real_t *p = calloc(n, sizeof(real_t));
    real_t *tmp = calloc(n, sizeof(real_t));    
    real_t alpha;
    real_t beta;
    double err = 0.0;
    int iter = 0;

    criaKDiagonal(n, k, A, b);
    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tempo);
    geraDLU(A, n, k, D, L, U, &tempo);
    geraPreCond(D, L, U, w, n, k, M, &tempo);    
    
    //APLICAR M EM A e b
    multiplicar_matrizes(M, ASP, A, n, n, n);
    multiplicar_matrizes(M, bsp, b, n, n, 1);

    calcResiduoSL(A, b, X, r, n, k, &tempo);
    copiar_vetor(r, p, n);
    do{
        copiar_vetor(X, x_old, n);
        copiar_vetor(r, r_old, n);
        alpha = escalar(r, n);
        multiplicar_matrizes(A, p, tmp, n, n, 1);
        alpha /= escalar(tmp, n);
        for(int i=0; i<n;i++){
            X[i] = X[i] + alpha*p[i];
            r[i] = r[i] - alpha*tmp[i];
        }
        err = calc_err(X, x_old, n);
        if(err < max_err){
            break;
        }
        beta = escalar(r, n)/escalar(r_old, n);
        for(int i=0; i<n; i++){
            p[i] = r[i] + beta*p[i];
        }
        iter++;
    }while(iter < maxiter);
    test_prints(A, b, X, ASP, bsp, D, L, U, r, n);
    printf("%d\n", iter);
}