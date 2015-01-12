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

#ifndef INTEGRALS_1EL_POTENTIAL
#define INTEGRALS_1EL_POTENTIAL

#include "basisinfo.h"
#include "basis_func_pair_list_1el.h"


ergo_real simplePrimVintegral_list(DistributionSpecStruct* list,
				   int nPrims,
				   const Atom* atom,
				   ergo_real threshold,
				   const IntegralInfo & integralInfo);

int compute_V_matrix_full(const BasisInfoStruct& basisInfo,
			  const IntegralInfo& integralInfo,
			  int nAtoms,
			  const Atom* atomList,
			  ergo_real threshold,
			  ergo_real* result);

int compute_V_and_gradient_linear(const BasisInfoStruct& basisInfo,
				  const IntegralInfo& integralInfo,
				  const Molecule& molecule,
				  ergo_real threshold,
				  ergo_real boxSize,
				  const basis_func_index_pair_struct_1el* basisFuncIndexPairList,
				  ergo_real* V_list,
				  int noOfBasisFuncIndexPairs,
				  bool compute_gradient_also,
				  const ergo_real* D_list,         // List of corresponding density matrix elemets; used for compute_gradient_also case, NULL otherwise
				  ergo_real* result_gradient_list  // list of result gradient values; used for compute_gradient_also case, NULL otherwise
				  );


#endif
