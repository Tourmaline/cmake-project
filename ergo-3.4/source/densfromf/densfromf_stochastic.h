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

#ifndef DENSFROMFSTOCHASTICHEADER
#define DENSFROMFSTOCHASTICHEADER

#include "realtype.h"
#include "matrix_typedefs.h"

int get_dens_from_fock_stochastic(int n, /**< System size. */
				  int noOfOccupiedOrbs, /**< Number of occupied orbitals. */
				  symmMatrix & resultDens, /**< (out) Resulting density matrix (D_S) in 'non-orthogonal basis'. 
							      ( D_S = Z*D_ort*ZT ) */ 
				  ergo_real factor, /**< Factor to scale the resulting density matrix. (for restricted vs unrestricted calc) */
				  symmMatrix const & Finput, /**< (in) Fock/Kohn-Sham matrix (F_S) in 'non-orthogonal basis'. 
								(written to file) */ 
				  triangMatrix const & invCholFactor, /**< (in) Inverse Cholesky factor of S. (written to file) */
				  mat::SizesAndBlocks const & matrixSizesAndBlocks, /**< Information about HML matrix block sizes etc. */
				  const std::vector< std::vector<ergo_real> > stochastic_orbitals /**< Vector of stochastic orbitals.  */
				  );

#endif
