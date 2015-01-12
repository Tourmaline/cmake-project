/* Ergo, version 3.4, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2014 Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

#include <stdlib.h>

#include "matrix_algebra.h"
#include "memorymanag.h"
#include "output.h"

#include "../matrix/mat_gblas.h"


#define USE_BLAS_MM


#if 0
#ifdef USE_BLAS_MM
#ifdef __cplusplus
extern "C" 
#endif
void dgemm_(const char *ta,const char *tb,
	    const int *n, const int *k, const int *l,
	    const double *alpha,const double *A,const int *lda,
	    const double *B, const int *ldb,
	    const double *beta,const double *C, const int *ldc);
#endif
#endif


void multiply_matrices_general_2(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB, bool initToZero) {
  if(An2 != Bn1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in multiply_matrices_general_2: (An2 != Bn1)");
    exit(0);
  }
#ifdef USE_BLAS_MM
  if(An1 == 0 || An2 == 0 || Bn1 == 0 || Bn2 == 0)
    return;
  /*  call gemm */
  ergo_real alpha = 1;
  ergo_real beta = 1;
  if(initToZero)
    beta = 0;
  mat::gemm("N", "N", &Bn2, &An1, &Bn1, &alpha, 
	    B, &Bn2, 
	    A, &An2, 
	    &beta, 
	    AB, &Bn2);
#else
  if(initToZero)
    memset(AB, 0, An1*Bn2*sizeof(ergo_real));
  for(int i = 0; i < An1; i++)
    for(int k = 0; k < An2; k++)
      for(int j = 0; j < Bn2; j++)
	AB[i*Bn2+j] += A[i*An2+k] * B[k*Bn2+j];
#endif
}



/*
  Standard matrix multiplication.
*/
void 
multiply_matrices_general(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB)
{
  int i, j;
  if(An2 != Bn1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, 
		"error in multiply_matrices_general: (An2 != Bn1)");
      exit(0);
    }
#ifdef USE_BLAS_MM
  if(An1 == 0 || An2 == 0 || Bn1 == 0 || Bn2 == 0)
    return;
  /*  call gemm */
  ergo_real alpha = 1;
  ergo_real beta = 0;
  ergo_real* ABtemp = (ergo_real*)ergo_malloc(An1*Bn2*sizeof(ergo_real));
  memset(ABtemp, 0, An1*Bn2*sizeof(ergo_real));
  mat::gemm("T", "T", &An1, &Bn2, &An2, &alpha, 
	    A, &An2, 
	    B, &Bn2, 
	    &beta, 
	    ABtemp, &An1);
  /*  do transpose of result */
  for(i = 0; i < An1; i++)
    for(j = 0; j < Bn2; j++)
      {
	AB[i*Bn2+j] = ABtemp[j*An1+i];
      }  
  ergo_free(ABtemp);
#else
  for(i = 0; i < An1; i++)
    for(j = 0; j < Bn2; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < An2; k++)
	  sum += A[i*An2+k] * B[k*Bn2+j];
	AB[i*Bn2+j] = sum;
      }
#endif
}


/*
  Matrix multiplication when the first matrix is transposed.
*/
void 
multiply_matrices_general_T_1(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB)
{
  int i, j;
  if(An1 != Bn1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in multiply_matrices_general_T_1: (An1 != Bn1)");
      exit(0);
    }
  ergo_real* At = (ergo_real*)ergo_malloc(An1*An2*sizeof(ergo_real));
  for(i = 0; i < An1; i++)
    for(j = 0; j < An2; j++)
      At[j*An1+i] = A[i*An2+j];
  multiply_matrices_general(An2, An1, Bn1, Bn2, At, B, AB);
  ergo_free(At);
}

/*
  Matrix multiplication when the second matrix is transposed.
*/
void 
multiply_matrices_general_T_2(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB)
{
  int i, j;
  if(An2 != Bn2)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in multiply_matrices_general_T_2: (An2 != Bn2)");
      exit(0);
    }
  ergo_real* Bt = (ergo_real*)ergo_malloc(Bn1*Bn2*sizeof(ergo_real));
  for(i = 0; i < Bn1; i++)
    for(j = 0; j < Bn2; j++)
      Bt[j*Bn1+i] = B[i*Bn2+j];
  multiply_matrices_general(An1, An2, Bn2, Bn1, A, Bt, AB);
  ergo_free(Bt);
}


void 
multiply2matrices(int n, ergo_real* A, ergo_real* B, ergo_real* AB)
{
#ifdef USE_BLAS_MM
  if(n == 0)
    return;
  /*  call dgemm */
  int i, j;
  int nn = n;
  ergo_real alpha = 1;
  ergo_real beta = 0;
  mat::gemm("T", "T", &nn, &nn, &nn, &alpha, 
	    A, &nn, 
	    B, &nn, 
	    &beta, 
	    AB, &nn);
  /*  do transpose of result */
  for(i = 0; i < n; i++)
    for(j = i; j < n; j++)
      {
	ergo_real temp = AB[i*n+j];
	AB[i*n+j] = AB[j*n+i];
	AB[j*n+i] = temp;
      }
#else
  int i, j, k;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(k = 0; k < n; k++)
	  sum += A[i*n+k] * B[k*n+j];
	AB[i*n+j] = sum;
      }
#endif
}


void 
multiply2matricesSymm(int n, ergo_real* A, ergo_real* B, ergo_real* AB)
{
  int i, j, k;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(k = 0; k < n; k++)
	  sum += A[i*n+k] * B[j*n+k];
	AB[i*n+j] = sum;
      }
}

void 
multiply2matricesSymmResult(int n, ergo_real* A, ergo_real* B, ergo_real* AB)
{
  int i, j, k;
  for(i = 0; i < n; i++)
    for(j = 0; j <= i; j++)
      {
	ergo_real sum = 0;
	for(k = 0; k < n; k++)
	  sum += A[i*n+k] * B[k*n+j];
	AB[i*n+j] = sum;
      }
  for(i = 0; i < n; i++)
    for(j = i+1; j < n; j++)
      AB[i*n+j] = AB[j*n+i];
}


void 
computeSquareOfSymmetricMatrix(int n, const ergo_real* Aa, const ergo_real* Ab, ergo_real* A2)
{
  int i, j;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	int k;
	for(k = 0; k < n; k++)
	  sum += Aa[i*n+k] * Ab[j*n+k];
	A2[i*n+j] = sum;
      }
  return;
}


void multiply3matrices(int n, ergo_real* A, ergo_real* B, ergo_real* C, ergo_real* ABC)
{
  ergo_real* AB = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  multiply2matrices(n, A, B, AB);
  multiply2matrices(n, AB, C, ABC);
  ergo_free(AB);
}



