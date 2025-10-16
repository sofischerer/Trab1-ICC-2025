#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"

void criaKDiagonal(int n, int k, double *A, double *B);
void print_matriz(double* A, int lin, int col);

void genSimetricaPositiva(double *A, double *b, int n, int k, double *ASP, double *bsp, double *tempo);
void geraDLU (double *A, int n, int k, double **D, double **L, double **U, double *tempo);
void geraPreCond(double *D, double *L, double *U, double w, int n, int k, double **M, double *tempo);
double calcResiduoSL (double *A, double *b, double *X, int n, int k, double *tempo);

#endif // __SISLIN_H__

