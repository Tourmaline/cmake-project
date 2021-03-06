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

#ifndef DENSITYFITTING_HEADER
#define DENSITYFITTING_HEADER

#include <stdio.h>

#include "basisinfo.h"
#include "integrals_2el.h"


typedef struct {
  ergo_real *ptr;
  FILE *f;
  unsigned using_file:1;
} DensfitData;

DensfitData* densfit_init(const IntegralInfo* integralInfo,
			  const BasisInfoStruct & basisInfoDensFit);
void densfit_destroy(DensfitData *p);

int densfit_compute_alpha_beta_matrix_inverse(const IntegralInfo* integralInfo,
					      const BasisInfoStruct & basisInfoDensFit,
					      ergo_real* result_U_inverse);

int densfit_compute_gamma(const IntegralInfo* integralInfo,
			  const BasisInfoStruct & basisInfoMain,
			  const BasisInfoStruct & basisInfoDensFit,
			  ergo_real* densityMatrix,
			  ergo_real* result_gamma,
			  ergo_real threshold);

int densfit_compute_c_vector(const IntegralInfo* integralInfo,
			     const BasisInfoStruct & basisInfoDensFit,
			     DensfitData* U_inverse,
			     ergo_real* gamma,
			     ergo_real* result_c_vector);

int densfit_compute_J(const IntegralInfo* integralInfo,
		      const BasisInfoStruct & basisInfoMain,
		      const BasisInfoStruct & basisInfoDensFit,
		      ergo_real* c_vector,
		      ergo_real* result_J,
		      ergo_real threshold);

#endif /* DENSITYFITTING_HEADER */
