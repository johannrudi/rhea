#ifndef ADJOINT_MATH_H
#define ADJOINT_MATH_H

#include <stdlib.h>
#include <stdio.h>

int
int_pow (int i, int j);

void
print_square_matrix (double **A, int n);

double **
allocate_square_matrix (int n);

void
destroy_square_matrix (double **a, int n);

void
sub_cofactor_matrix (double **subcofA, double **A, int p, int q, int n);

double
determinant (double **A, int n);

void
cofactor_matrix (double **cofA, double **A, int n);

void
transpose_matrix (double **transA, double **A, int n);

void
adjoint_matrix (double **adjA, double **A, int n);

void
inverse_matrix (double **invA, double **A, int n);

void
matrix_multiplies_vector (double *b, double **A, double *x, int n);

void
solve_x_from_Ax_b (double *x, double **A, double *b, int n);

#endif
