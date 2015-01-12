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

#ifndef DENSITY_PROJECTION
#define DENSITY_PROJECTION

#include "basisinfo.h"
#include "matrix_typedefs.h"


/** load_density_and_project_full loads one or two density matrices (depending
    on value of noOfDensityMatrices) from the file specified by
    densityFileName. 
    @param densityFileName Name of file to load density matrices from
    @param noOfDensityMatrices Number of density matrices to load
    @param integralInfo static helper object for
    integral evaluation. 
    @param basisInfo the basis set that the density is
    to be expanded into (if the original density was expressed with
    help of other basis set, an apriopriate projection will be
    performed). 
    @param densityMatrixList must be already allocated and have
    proper dimension. 
    @param do_purification determines whether an additional
    purification is to be run after the projection. 
    @param noOfElectronsList is an one or two element array specyfying 
    the number of total electrons (one element) 
    or alpha and beta electrons (two elements).
    @param electronic_temperature Electronic temperature
*/
int load_density_and_project_full(const char *densityFileName,
				  int noOfDensityMatrices,
				  const IntegralInfo* integralInfo,
				  const BasisInfoStruct & basisInfo,
				  ergo_real** densityMatrixList,
				  int do_purification,
				  const int* noOfElectronsList,
				  ergo_real electronic_temperature);


/** load_density_and_project_sparse loads one or two density matrices (depending
    on value of noOfDensityMatrices) from the file specified by
    densityFileName. 

    The projection is done as follows:
    First, a matrix R is computed. 
    R is the overlap matrix between the two basis sets.
    Then RT * P * R is computed, where P is the starting guess density matrix
    read from file.
    To get a final projected density one could then multiply by S_inv from
    both sides, but to prepare for purification the matrix S*D*S is needed,
    so we skip multiplication by S_inv since it will anyway be cancelled out.

    @param densityFileName Name of file to load density matrices from
    @param noOfDensityMatrices Number of density matrices to load
    @param integralInfo static helper object for
    integral evaluation. 
    @param basisInfo the basis set that the density is
    to be expanded into (if the original density was expressed with
    help of other basis set, an apriopriate projection will be
    performed). 
    @param S_symm Overlap matrix
    @param densityMatrixList pointers to one or two empty matrices that will 
    contain the result.
    Purification is always run after the projection. 
    @param noOfElectronsList is an one or two element array specyfying 
    the number of total electrons (one element) 
    or alpha and beta electrons (two elements).
    @param matrix_size_block_info Information about HML matrix block sizes etc.
    @param matrixPermutationVec Permutation vector used when calling matrix lib.
    @param sparse_threshold Threshold used when truncating matrices.
    @param invCholFactor Inverse Cholesky factor of S.
    @param invCholFactor_euclnorm Euclidean norm of inverse Cholesky factor.
    @param gap_expected_lower_bound Expected lower bound for the band gap.
    @param purification_eigvalue_err_limit Requested accuracy in eigenvalues of D_ort. 
    @param purification_subspace_err_limit Requested accuracy in the occupied subspace.
    @param purification_truncation_norm Norm to be used for truncation in purification.
    @param purification_maxmul Maximum allowed number of matrix multiplications.
    @param purification_create_m_files Flag to create m-files with information about the purification process.
    @param use_diagonalization Flag to turn on diagonalization.
    @param use_diag_on_error Flag to fall back on diagonalization if purification fails.
    @param purification_ignore_failure Continue even if purification fails to converge.
    @param purification_use_rand_perturbation_for_alleigsint Apply a random perturbation to (try to) improve the convergence speed of Lanczos calculation of extremal eigenvalues.
    @param electronic_temperature Electronic temperature.
*/
int 
load_density_and_project_sparse(const char *densityFileName,
				int noOfDensityMatrices,
				const IntegralInfo* integralInfo,
				const BasisInfoStruct & basisInfo,
				symmMatrix & S_symm,
				symmMatrix** densityMatrixList,
				const int* noOfElectronsList,
				mat::SizesAndBlocks matrix_size_block_info,
				std::vector<int> const & matrixPermutationVec,
				ergo_real sparse_threshold,
				triangMatrix & invCholFactor,
				ergo_real invCholFactor_euclnorm,
				ergo_real gap_expected_lower_bound,
				ergo_real purification_eigvalue_err_limit,
				ergo_real purification_subspace_err_limit,
				mat::normType const purification_truncation_norm,
				int purification_maxmul,
				int purification_create_m_files,
				int use_diagonalization,
				int use_diag_on_error,
				int purification_ignore_failure,
				int purification_use_rand_perturbation_for_alleigsint,
				ergo_real electronic_temperature);




#endif /*  DENSITY_PROJECTION */
