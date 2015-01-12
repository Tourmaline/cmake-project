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

#include <memory.h>
#include "solve_lin_eq_syst.h"
#include "mat_gblas.h"
#include "output.h"


int 
solve_linear_equation_system(int n, 
			     const ergo_real* matrix, 
			     const ergo_real* RHS, 
			     ergo_real* resultVector)
{
  integer NRHS = 1;
  integer n2 = n;
  integer info;
  ergo_real* A = new ergo_real[n*n];
  integer* IPIV = new integer[n];
  
  memcpy(A, matrix, n*n*sizeof(ergo_real));
  memcpy(resultVector, RHS, n*sizeof(ergo_real));
  
  template_lapack_gesv(&n2, &NRHS, A, &n2, IPIV, resultVector, &n2, &info);
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "ERROR in dgesv_");
      return -1;
    }

  delete [] A;
  delete [] IPIV;

  return 0;
}





#if 0

int 
main()
{
  printf("hello\n");

  integer n = 3;

  ergo_real A[n*n];

  A[0] = 1;
  A[1] = 2;
  A[2] = 0;

  A[3] = 0;
  A[4] = 3;
  A[5] = 3;

  A[6] = 1;
  A[7] = 1;
  A[8] = 0;

  ergo_real RHS[n];
  RHS[0] = 1;
  RHS[1] = 0;
  RHS[2] = 1;

  integer IPIV[n];
  integer NRHS = 1;

  integer info = -1;

  dgesv_(&n, &NRHS, A, &n, IPIV, RHS, &n, &info);
  if(info != 0)
    {
      printf("ERROR!\n");
      return -1;
    }
  printf("OK\n");

  printf("solution:\n");
  int i;
  for(i = 0; i < n; i++)
    printf("%22.11f\n", (double)RHS[i]);
  
  return 0;
}

#endif
