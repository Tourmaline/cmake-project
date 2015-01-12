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
 
 /* This file belongs to the template_lapack part of the Ergo source 
  * code. The source files in the template_lapack directory are modified
  * versions of files originally distributed as CLAPACK, see the
  * Copyright/license notice in the file template_lapack/COPYING.
  */
 

#include <stdexcept>
#include <iostream>
#include <math.h>
#include <string.h>
#include <pthread.h>

#include "template_lapack_common.h"

typedef double realtype;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

#if 0
static void output_matrix(int n, const realtype* A, const char* name)
{
  std::cout << "output_matrix for matrix '" << name << "'\n";
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
	std::cout << "  " << A[i*n+j];
      std::cout << std::endl;
    }
}
#endif

static void*
thread_func(void* arg)
{
  for(int loop_k = 0; loop_k < 4444; loop_k++) {

    realtype A[4];
    A[0] = 1.3;  A[1] = 2.4;
    A[2] = 1.1;  A[3] = 2.7;
  
    realtype B[4];
    B[0] = 0.3;  B[1] = 0.4;
    B[2] = 0.1;  B[3] = 0.7;
  
    realtype C[4];
    C[0] = 0;  C[1] = 0;
    C[2] = 0;  C[3] = 0;
  
    int m = 2;
    int n = 2;
    int k = 2;
    realtype alpha = 1.0;
    realtype beta = 0.0;
    template_blas_gemm("N", "N", &m, &n, &k, &alpha, A, &m, B, &m, &beta, C, &m);

    //  std::cout << "template_blas_gemm finished.\n";

#if 0
    for(int i = 0; i < n; i++)
      {
	for(int j = 0; j < n; j++)
	  std::cout << "  " << C[i*n+j];
	std::cout << std::endl;
      }
#endif

    realtype D[4];
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
	{
	  realtype sum = 0;
	  for(int k = 0; k < n; k++)
	    sum += B[i*n+k] * A[k*n+j];
	  D[i*n+j] = sum;
	}
      
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
	{
	  if(fabs(C[i*n+j] - D[i*n+j]) > 1e-11)
	    {
	      throw std::runtime_error("ERROR! wrong result from gemm!");
	      return NULL;
	    }
	}




    realtype E[4];
    E[0] = 1.3;  E[1] = 0.4;
    E[2] = 0.4;  E[3] = 2.7;

    int info = 0;


    //output_matrix(n, E, "E before pptrf");
    template_lapack_pptrf("U", &n, E, &info);
    //output_matrix(n, E, "E after pptrf");

    int itype = 1;

    //  std::cout << "calling template_lapack_spgst\n";

    template_lapack_spgst(&itype, "U", &n, 
			  A, B, &info);

    //  std::cout << "calling template_lapack_tptri\n";

    template_lapack_tptri("U", "U", &n, 
			  A, &info);





    int ITYPE=1;
    int lwork = 10*n, i;
    double work[1000];
    double eigv[1000];

    //  std::cout << "calling template_lapack_sygv\n";

    template_lapack_sygv(&ITYPE, "V", "L", &n, A, &n, B, &n, 
			 eigv, work, &lwork, &i);

    //  std::cout << "calling template_lapack_trtri\n";

    char uplo[8];
    char diag[8];
    strcpy(uplo, "U");
    strcpy(diag, "U");
    template_lapack_trtri(uplo, diag, &n, A, &n, &info);








    //  std::cout << "testing template_lapack_gesv().." << std::endl;

    {
      integer n = 3;

      double A[n*n];

      A[0] = 1;
      A[1] = 2;
      A[2] = 0;

      A[3] = 0;
      A[4] = 3;
      A[5] = 3;

      A[6] = 1;
      A[7] = 1;
      A[8] = 0;

      double RHS[n];
      RHS[0] = 1;
      RHS[1] = 0;
      RHS[2] = 1;

      integer IPIV[n];
      integer NRHS = 1;

      integer info = -1;

      template_lapack_gesv(&n, &NRHS, A, &n, IPIV, RHS, &n, &info);
      if(info != 0)
	{
	  throw std::runtime_error("ERROR in template_lapack_gesv!");
	  return NULL;
	}
      //    std::cout << "OK"  << std::endl;

    double ref_solution[] = {-2, 0.33333333333333, 3};
    for(int i = 0; i < n; i++) {
      if(fabs(RHS[i] - ref_solution[i]) > 1e-11)
	throw std::runtime_error("ERROR: wrong solution from template_lapack_gesv!");
    }

#if 0
      std::cout << "solution:    ";
      for(int i = 0; i < n; i++)
	std::cout << RHS[i] << "    ";
      std::cout   << std::endl;
#endif
    }

  }

  pthread_mutex_lock(&mutex);
  std::cout << "Thread " << pthread_self() << " finished OK." << std::endl;
  pthread_mutex_unlock(&mutex);

  return NULL;
}



int main()
{
  const int nThreads = 8;
  pthread_t threads[nThreads];
  for(int i = 0; i < nThreads; i++) {
    if(pthread_create(&threads[i], NULL, thread_func, NULL) != 0)
      throw std::runtime_error("Error in pthread_create.");
  }
  for(int i = 0; i < nThreads; i++) {
    if(pthread_join(threads[i], NULL) != 0)
      throw std::runtime_error("Error in pthread_join.");
  }
  
  std::cout << "template_lapack_threaded test finished OK." << std::endl;
  
  return 0;
}

