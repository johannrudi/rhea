#include <adjoint_math.h>

int
int_pow(int base, int exp)
{
  int result = 1;
  while (exp)
  {
    if (exp & 1)
      result *= base;
      exp /= 2;
      base *= base;
  }
  return result;
}

void
print_square_matrix (double **A, int n)
{
  int i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf("%10.4f", A[i][j]);
    }
    printf("\n");
  }
}

double **
allocate_square_matrix (int n)
{
  double **A;

  A = (double **) malloc(n * sizeof(double *));
  for (int i = 0; i < n; i++)
    A[i] = (double *) malloc (n * sizeof(double));

  return (A);
}

void
destroy_square_matrix (double **A, int n)
{
  for (int i = 0; i < n; i++) {
    free (A[i]);
  }
  free(A);
}

// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
void
sub_cofactor_matrix (double **subcofA, double **A, int p, int q, int n)
{
  int i = 0, j = 0;
  for (int row = 0; row < n; row++)  {
    for (int col = 0; col < n; col++)  {
      if (row != p && col != q)  {
        //  Copying into temporary matrix only those element
        //  which are not in given row or column
        subcofA[i][j++] = A[row][col];
        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}

double
determinant (double **A, int n)
{
  int sign = 1;
  double det = 0, **subcofA;
  int i;
  //  Base case : if matrix contains single element
  if (n == 1) {
    return A[0][0];
  }

  subcofA = allocate_square_matrix (n - 1);
  // Iterate for each element of first row
  for (i = 0; i < n; i++)  {
    // Getting sub_cofactor of A[0][n]
    sub_cofactor_matrix (subcofA, A, 0, i, n);
    det += sign * A[0][i] * determinant (subcofA, n - 1);
    // terms are to be added with alternate sign
    sign *= -1;
  }

  destroy_square_matrix (subcofA, n - 1);

  return (det);
}

void
cofactor_matrix (double **cofA, double **A, int n)
{
  int   i, j;
  int   sign = 1;
  double **subcofA;

  if (n == 1)  {
    cofA[0][0] = 1;
    return;
  }

  subcofA = allocate_square_matrix (n - 1);
  for (i = 0; i < n; i++)  {
    for (j = 0; j < n; j++)  {
      sub_cofactor_matrix (subcofA, A, i, j, n);
      sign = int_pow (-1, i+j);
      cofA[i][j] = sign * determinant (subcofA, n - 1);
    }
  }
  destroy_square_matrix (subcofA, n - 1);
}

void
transpose_matrix (double **transA, double **A, int n)
{
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      transA[i][j] = A[j][i];
    }
  }

}

void
adjoint_matrix (double **adjA, double **A, int n)
{
  double **cofA;

  if (n == 1)  {
    adjA[0][0] = 1;
    return;
  }

  cofA = allocate_square_matrix (n);
  cofactor_matrix (cofA, A, n);
  transpose_matrix (adjA, cofA, n);

  destroy_square_matrix (cofA, n);
}

void
inverse_matrix (double **invA, double **A, int n)
{
  double det;

  if (n == 1) {
    invA[0][0] = 1.0 / A[0][0];
    return;
  }

  det = determinant(A, n);
  if (det == 0) {
    fprintf(stderr, "determinant of A matrix is 0!\n");
    return;
  }
  else {
    adjoint_matrix (invA, A, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        invA[i][j] /= det;
  }
}

void
matrix_multiplies_vector (double *b, double **A, double *x, int n)
{
  int i, j;
  for (i = 0; i < n; i++) {
    b[i] = .0;
    for (j = 0; j < n; j++) {
      b[i] += A[i][j] * x[j];
    }
  }
}

void
solve_x_from_Ax_b (double *x, double **A, double *b, int n)
{
  double **invA;

  invA = allocate_square_matrix (n);
  inverse_matrix (invA, A, n);
  matrix_multiplies_vector (x, invA, b, n);

  destroy_square_matrix (invA, n);
}
