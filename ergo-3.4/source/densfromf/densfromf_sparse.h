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

#ifndef DENSFROMFSPARSEHEADER
#define DENSFROMFSPARSEHEADER

#include "realtype.h"
#include "matrix_typedefs.h"

int get_dens_from_fock_sparse(int n, /**< System size. */
			      int noOfOccupiedOrbs, /**< Number of occupied orbitals. */
			      symmMatrix & resultDens, /**< (out) Resulting density matrix (D_S) in 'non-orthogonal basis'. 
							  ( D_S = Z*D_ort*ZT ) */ 
			      ergo_real factor, /**< Factor to scale the resulting density matrix. (for restricted vs unrestricted calc) */
			      symmMatrix const & Finput, /**< (in) Fock/Kohn-Sham matrix (F_S) in 'non-orthogonal basis'. 
							    (written to file) */ 
			      intervalType & homoInterval_Finput, /**< (out)
								     Output: Contains the homo eigenvalue of Finput. */
			      intervalType & lumoInterval_Finput, /**< (out)
								     Output: Contains the lumo eigenvalue of Finput. */
			      triangMatrix const & invCholFactor, /**< (in) Inverse Cholesky factor of S. (written to file) */
			      ergo_real invCholFactor_euclnorm, /**< Euclidean norm of inverse Cholesky factor. */
			      ergo_real gap_expected_lower_bound, /**< Expected lower bound for the gap to be used in early iterations. */
			      mat::SizesAndBlocks const & matrixSizesAndBlocks, /**< Information about HML matrix block sizes etc. */
			      symmMatrix & F_ort_prev, /**< (in/out) 
							  Input: Previous F matrix in orthogonal basis. (written to file)
							  Output: New F matrix in orthogonal basis ( ZT*Finput*Z ). (written to file) */
			      intervalType & homoInterval_F_ort_prev, /**< (in/out)
									 Input: Contains the homo eigenvalue of F_ort_prev.
									 Output: Contains the homo eigenvalue of F_ort_prev. */
			      intervalType & lumoInterval_F_ort_prev, /**< (in/out)
									 Input: Contains the lumo eigenvalue of F_ort_prev.
									 Output: Contains the lumo eigenvalue of F_ort_prev. */
			      ergo_real eigvalueErrorLimit, /**< (in) Requested accuracy in eigenvalues of D_ort. */
			      ergo_real subspaceErrorLimit, /**< (in) Requested accuracy in the occupied subspace of D_ort. */
			      mat::normType const truncationNormPurification, /**< Norm to be used for truncation in, before, 
										 and after purification. */
			      int maxMul, /**< Maximum allowed number of matrix multiplications. */
			      int create_m_files, /**< Flag to create m-files with information about the purification process. */
			      int ignore_purification_failure, /**< Continue even if purification fails to converge. */
			      int use_rand_perturbation_for_alleigsint, /**< Apply a random perturbation to (try to) improve the convergence speed of Lanczos calculation of extremal eigenvalues.  */
			      std::string stats_prefix, /**<  Prefix to be added to statistics files. */
			      std::map<std::string, double> & puri_stats, /**< Map to store stats for purification. */
			      int do_sparsity_investigation, /**< Flag to turn on sparsity investigation. */
			      int sparsity_plots_resolution_m, /**< Resolution in element magnitude histograms. */
			      int do_comparison_to_simple_purification, /**< Flag to turn on comparison to simple purification. */			      
			      int do_puri_mmul_tests, /**< Flag to turn on purification matrix-matrix multiplication tests. */
			      generalVector * eigVecLUMO = 0, /**<  LUMO eigenvector */
			      generalVector * eigVecHOMO = 0  /**<  HOMO eigenvector */
			      );

#endif
