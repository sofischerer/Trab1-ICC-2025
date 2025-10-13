#include "sislin.h"
#include "utils.h"
#include <time.h>
#include <stdlib.h>

void print_matriz(real_t** A, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}

int main(){
    srandom(20252);
    
    int k = 5;
    int n = 10;
    real_t** A = calloc(n, sizeof(real_t*));
    real_t** B = NULL;
    criaKDiagonal(n, k, A, B);


}