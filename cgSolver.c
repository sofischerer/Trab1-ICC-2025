#include "pcgc.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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

real_t residuo_euc(real_t* A, real_t n, rtime_t* tempo){
    
    *tempo = timestamp();
    real_t sum = 0.0;
    for(int i=0;i<n;i++){
        sum += A[i]*A[i];
    }
    *tempo = timestamp() - *tempo;
    return sqrt(sum);
}


int main(){

    srandom(20252);    
    int k, n, maxiter;
    real_t w;
    real_t max_err;

    scanf("%d %d %lf %d %lf", &n, &k, &w, &maxiter, &max_err);

    if (n<10){
        fprintf(stderr, "Erro: N muito pequeno (menor que 10)");
        return -1;
    }
    if (k<=1){
        fprintf(stderr, "Erro: K muito pequeno (menor que 3)");
        return -1;
    }
    if (k%2 == 0){
        fprintf(stderr, "Erro: K par");
        return -1;
    }
    if(w!=-1 && w!=0){
        fprintf(stderr, "Erro: w inválido");
        return -1;
    }

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
    rtime_t tempo_pc, tempo_iter, tempo_residuo, temp;
    double err = 0.0;
    tempo_iter = 0.0;
    tempo_pc = 0.0;
    int iter = 0;
    //Podia limpar isso mas não estou com tempo para isso, para o T2 será limpo e organizado

    criaKDiagonal(n, k, A, b);
    genSimetricaPositiva(A, b, n, k, ASP, bsp, &temp);
    tempo_pc += temp;
    geraDLU(A, n, k, D, L, U, &temp);
    tempo_pc += temp;
    geraPreCond(D, L, U, w, n, k, M, &temp);    
    tempo_pc += temp;
    
    //APLICAR M EM A e b
    multiplicar_matrizes(M, ASP, A, n, n, n);
    multiplicar_matrizes(M, bsp, b, n, n, 1);

    calcResiduoSL(A, b, X, r, n, k);
    copiar_vetor(r, p, n);
    temp = 0.0;
    do{
        temp = timestamp();
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
        temp = timestamp() - temp;
        tempo_iter += temp;
    }while(iter < maxiter);
    if (iter >= maxiter){
        fprintf(stderr, "Método não convergiu no limite de iterações (%d)", iter);
        return -1;
    }
    tempo_iter /= iter;
    printf("%d\n", n);
    print_matriz( X, 1, n);
    printf("%.8g\n", err);
    printf("%.16g\n", residuo_euc(r, n, &tempo_residuo));
    printf("%.8g\n", tempo_pc);
    printf("%.8g\n", tempo_iter);
    printf("%.8g\n", tempo_residuo);
}