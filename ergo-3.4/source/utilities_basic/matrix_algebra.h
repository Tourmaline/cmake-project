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

#include "realtype.h"

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

EXTERN_C void multiply2matrices(int n, ergo_real* A, ergo_real* B, ergo_real* AB);
EXTERN_C void multiply2matricesSymm(int n, ergo_real* A, ergo_real* B, ergo_real* AB);
EXTERN_C void multiply2matricesSymmResult(int n, ergo_real* A, ergo_real* B, ergo_real* AB);
EXTERN_C void computeSquareOfSymmetricMatrix(int n, 
				    const ergo_real* Aa, 
				    const ergo_real* Ab, 
				    ergo_real* A2);
EXTERN_C void multiply_matrices_general    (int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB);
EXTERN_C void multiply_matrices_general_2  (int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB, bool initToZero);
EXTERN_C void multiply_matrices_general_T_1(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB);
EXTERN_C void multiply_matrices_general_T_2(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB);
EXTERN_C void multiply3matrices(int n, ergo_real* A, ergo_real* B, ergo_real* C, ergo_real* ABC);

