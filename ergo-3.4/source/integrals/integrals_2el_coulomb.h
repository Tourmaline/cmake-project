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

#ifndef INTEGRALS_2EL_COULOMB_HEADER
#define INTEGRALS_2EL_COULOMB_HEADER


#include "basisinfo.h"
#include "integrals_2el.h"
#include "basis_func_pair_list.h"


int
compute_J_by_boxes(const BasisInfoStruct & basisInfo,
		   const IntegralInfo & integralInfo,
		   const JK::Params& J_K_params,
		   ergo_real* J,
		   const ergo_real* dens);

int
compute_J_by_boxes_nosymm(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  const JK::Params& J_K_params,
			  ergo_real* J,
			  const ergo_real* dens);

int
compute_J_by_boxes_linear(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  const JK::Params& J_K_params,
			  const basis_func_index_pair_struct* basisFuncIndexPairList,
			  int basisFuncIndexPairCount,
			  const ergo_real* D_list,
			  ergo_real* result_J_list,
			  int noOfBasisFuncIndexPairs);


#endif
