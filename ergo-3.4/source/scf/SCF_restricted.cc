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

#include <sstream>
#include "SCF_restricted.h"
#include "output.h"
#include "scf_utils.h"
#include "utilities.h"
#include "diis_restricted.h"
#include "density_projection.h"
#include "densfromf_general.h"
#include "density_description_file.h"
#include "matrix_utilities.h"
#include "matrix_utilities_gao.h"
#include "machine_epsilon.h"
#include "units.h"
#include "atom_labels.h"
#include "integral_matrix_wrappers.h"


SCF_restricted::SCF_restricted(const Molecule& molecule_,
			       const Molecule& extraCharges_,
			       const BasisInfoStruct & basisInfo_, 
			       const BasisInfoStruct & basisInfoDensFit_,
			       const IntegralInfo& integralInfo_,
			       const char* guessDmatFileNamePtr,
			       const JK::Params& J_K_paramsPtr,
			       const Dft::GridParams& gridParams_,
			       const SCF::Options& scfopts,
			       const SCF::MatOptions& matOpts_,
			       ergo_real threshold_integrals_1el_input)
  :   SCF_general(molecule_,
		  extraCharges_,
		  basisInfo_, 
		  basisInfoDensFit_,
		  integralInfo_,
		  guessDmatFileNamePtr,
		  J_K_paramsPtr,
		  gridParams_,
		  scfopts,
                  matOpts_,
		  threshold_integrals_1el_input)
{
  DIIS = new DIISManagerRestricted;
}


SCF_restricted::~SCF_restricted()
{
  delete ((DIISManagerRestricted*)DIIS);
}


void SCF_restricted::get_Fock_matrix(symmMatrix & FockMatrix_)
{
  FockMatrix.readFromFile();
  FockMatrix_ = FockMatrix;
  FockMatrix.writeToFile();
}


void SCF_restricted::get_density_matrix(symmMatrix & densityMatrix_)
{
  densityMatrix.readFromFile();
  densityMatrix_ = densityMatrix;
  densityMatrix.writeToFile();
}


void SCF_restricted::initialize_matrices()
{
  densityMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);

  densityMatrix_core.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);

  twoel_matrix_core.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);

  FockMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  Fprev.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  Dprev.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  F_ort_prev.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  bestFockMatrixSoFar.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  bestFockMatrixSoFar2.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  ErrorMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
}


void SCF_restricted::check_params()
{
  if(noOfElectrons%2 != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error: odd number of electrons in restricted calculation");
      throw "error: odd number of electrons in restricted calculation";
    }
}


void SCF_restricted::get_starting_guess_density()
{
  // set up starting guess

  int n = basisInfo.noOfBasisFuncs;

  if(guessDmatFileName != NULL)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "getting starting guess density from file '%s'", guessDmatFileName);
      int noOfDensityMatrices = 1;
      symmMatrix* matrixList[2];
      matrixList[0] = &densityMatrix;

      if(load_density_and_project_sparse(guessDmatFileName,
					 noOfDensityMatrices,
					 &integralInfo, 
					 basisInfo,
					 S_symm,
					 matrixList,
					 &noOfElectrons,
					 matOpts.size_block_info,
					 matOpts.permutationHML,
					 matOpts.sparse_threshold,
					 invCholFactor,
					 invCholFactor_euclnorm,
					 scfopts.gap_expected_lower_bound,
					 scfopts.purification_eigvalue_err_limit * scfopts.puri_eig_acc_factor_for_guess,
					 scfopts.purification_subspace_err_limit,
					 scfopts.purification_truncation_norm,
					 scfopts.purification_maxmul,
					 scfopts.purification_create_m_files,
					 scfopts.use_diagonalization,
					 scfopts.use_diag_on_error_guess,
					 scfopts.purification_ignore_failure,
					 scfopts.purification_use_rand_perturbation_for_alleigsint,
					 scfopts.electronic_temperature) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in load_density_and_project_sparse");
	  throw "error in load_density_and_project_sparse";
	}
    }
  else
    {
      if(scfopts.use_simple_starting_guess == 1)
	{
	  if(get_simple_starting_guess_sparse(n, noOfElectrons, densityMatrix) != 0)
	    throw "error in get_simple_starting_guess_sparse";
	  densityMatrix.writeToFile();
	}
      else if(scfopts.use_diag_guess_from_file == 1)
	{
	  if(get_diag_matrix_from_file(n, densityMatrix, "diagdens.txt",
				       matOpts.permutationHML) != 0)
	    throw "error in get_diag_matrix_from_file";
	  densityMatrix.writeToFile();
	}
      else
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		    "calling get_dens_from_fock_general to diagonalize H_core for starting guess, n = %i, sparse_threshold = %g", 
		    n, (double)matOpts.sparse_threshold);
	  
	  symmMatrix F_ort_prev_dummy;
	  F_ort_prev_dummy.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
	  F_ort_prev_dummy.writeToFile();
	  intervalType homoInterval_dummy1(-1e22,1e22);
	  intervalType lumoInterval_dummy1(-1e22,1e22);
	  intervalType homoInterval_dummy2(-1e22,1e22);
	  intervalType lumoInterval_dummy2(-1e22,1e22);
	  densityMatrix.writeToFile();
	  std::map<std::string, double> puri_stats;
	  ergo_real entropyTermDummy = 0;
	  if(get_dens_from_fock_general(n,
					noOfElectrons / 2,
					scfopts.use_diagonalization,
					scfopts.use_diag_on_error_guess,
					scfopts.electronic_temperature,
					densityMatrix,
					2,
					entropyTermDummy,
					H_core_Matrix,
					homoInterval_dummy1,
					lumoInterval_dummy1,
					S_symm,
					invCholFactor,
					invCholFactor_euclnorm,
					scfopts.gap_expected_lower_bound,
					matOpts.size_block_info,
					F_ort_prev_dummy,
					homoInterval_dummy2,
					lumoInterval_dummy2,
					scfopts.purification_eigvalue_err_limit,
					scfopts.purification_subspace_err_limit,
					scfopts.purification_truncation_norm,
					scfopts.purification_maxmul,
					scfopts.purification_create_m_files,
					scfopts.purification_ignore_failure,
					scfopts.purification_use_rand_perturbation_for_alleigsint,
					scfopts.use_stochastic_orbs,
					stochastic_orbitals,
					"",
					puri_stats,
					scfopts.do_sparsity_investigation,
					scfopts.sparsity_plots_resolution_m,
					scfopts.do_comparison_to_simple_purification,
					scfopts.do_puri_mmul_tests) != 0)
	    {
	      throw "SCF_restricted::get_starting_guess_density: Error in get_dens_from_fock_general";
	    }
	} // END ELSE use H_core
    } // END ELSE no dmat given

  densityMatrix.readFromFile();
  output_sparsity_symm(n, densityMatrix, "starting guess density matrix");
  densityMatrix.writeToFile();

  densityMatrix_core.writeToFile(); // densityMatrix_core (if used) also needs to be written to file
}


void SCF_restricted::add_random_disturbance_to_starting_guess()
{
  if(scfopts.sg_disturb_specific_elements > SCF::DISTURB_ELEMENT_MAX_COUNT)
    throw "Error in SCF_restricted::add_random_disturbance_to_starting_guess: (scfopts.sg_disturb_specific_elements > SCF::DISTURB_ELEMENT_MAX_COUNT)";
  int n = basisInfo.noOfBasisFuncs;
  densityMatrix.readFromFile();
  add_disturbance_to_matrix(n, 
			    densityMatrix, 
			    scfopts.starting_guess_disturbance,
			    scfopts.sg_disturb_specific_elements,
			    scfopts.disturbedElementIndexVector,
			    matOpts.permutationHML);
  densityMatrix.writeToFile();
}


void SCF_restricted::initialize_homo_lumo_limits()
{
  intervalType hugeInterval(-1e22, 1e22);
  homoInterval_F_ort_prev = hugeInterval;
  lumoInterval_F_ort_prev = hugeInterval;
  homoInterval_Fprev = hugeInterval;
  lumoInterval_Fprev = hugeInterval;
}


void SCF_restricted::write_matrices_to_file()
{
  FockMatrix.writeToFile();
  Fprev.writeToFile();
  Dprev.writeToFile();
  bestFockMatrixSoFar.writeToFile();
  bestFockMatrixSoFar2.writeToFile();
  F_ort_prev.writeToFile();
}


static void output_diff_norm_values(symmMatrix const & F1, 
				    symmMatrix const & F2, 
				    ergo_real acc,
				    const char* name) {
  ergo_real E_norm_frob = symmMatrix::frob_diff(F1, F2);
  Util::TimeMeter timeMeterMixedDiff;
  ergo_real E_norm_mixed = symmMatrix::mixed_diff(F1, F2, acc);
  timeMeterMixedDiff.print(LOG_AREA_DENSFROMF, "symmMatrix::mixed_diff");
  Util::TimeMeter timeMeterEuclDiff;
  ergo_real E_norm_eucl = symmMatrix::eucl_diff(F1, F2, acc);
  timeMeterEuclDiff.print     (LOG_AREA_DENSFROMF, "symmMatrix::eucl_diff ");
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Nor of error matrix for '%s':", name);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "frob :  %22.15f", (double)E_norm_frob);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "mixed:  %22.15f", (double)E_norm_mixed);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "eucl:   %22.15f", (double)E_norm_eucl);    
}


void SCF_restricted::get_2e_part_and_energy()
{
  densityMatrix.readFromFile();
  densityMatrix_core.readFromFile();

  bool scan_do_invcholfactor_transf = (bool)scfopts.scan_do_invcholfactor_transf;

  if(scfopts.do_acc_scan_J)
    do_acc_scan_J(densityMatrix,
		  integralInfo,
		  basisInfo,
		  invCholFactor,
		  scan_do_invcholfactor_transf,
		  J_K_params,
		  matOpts.size_block_info,
		  matOpts.permutationHML,
		  scfopts.scan_no_of_steps,
		  scfopts.scan_start_thresh,
		  scfopts.scan_step_factor);
  
  if(scfopts.do_acc_scan_K)
    do_acc_scan_K(densityMatrix,
		  integralInfo,
		  basisInfo,
		  invCholFactor,
		  scan_do_invcholfactor_transf,
		  CAM_params,
		  J_K_params,
		  matOpts.size_block_info,
		  matOpts.permutationHML,
		  matOpts.inversePermutationHML,
		  scfopts.scan_no_of_steps,
		  scfopts.scan_start_thresh,
		  scfopts.scan_step_factor);
  
  if(scfopts.do_acc_scan_Vxc)
    do_acc_scan_Vxc(densityMatrix,
		    integralInfo,
		    basisInfo,
		    molecule,
		    gridParams,
		    noOfElectrons,
		    invCholFactor,
		    scan_do_invcholfactor_transf,
		    matOpts.size_block_info,
		    matOpts.permutationHML,
		    matOpts.inversePermutationHML,
		    scfopts.scan_no_of_steps,
		    scfopts.scan_start_thresh,
		    scfopts.scan_step_factor);

  symmMatrix G;
  G.resetSizesAndBlocks(matOpts.size_block_info,
			matOpts.size_block_info);
  
  if(get_2e_matrix_and_energy_sparse(basisInfo, 
				     basisInfoDensFit, 
				     molecule,
				     integralInfo, 
				     G,
				     densityMatrix, 
				     J_K_params,
				     CAM_params,
				     gridParams,
				     scfopts.use_dft,
				     &energy_2el, noOfElectrons,
				     densfit_data,
				     matOpts.size_block_info,
				     matOpts.permutationHML,
				     matOpts.inversePermutationHML,
				     scfopts.do_sparsity_investigation,
				     J_matrix,
				     K_matrix,
				     Fxc_matrix,
				     *curr_cycle_stats) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrix_and_energy_sparse");
      throw "error in get_2e_matrix_and_energy_sparse";
    }

  if(scfopts.compute_core_density == 1) {
    // Get energy_2el_valence and H_matrix_core, but only if densityMatrix_core is nonzero (it is zero the first time, then we just skip).
    if(densityMatrix_core.frob() != 0) {
      // Get energy_2el_valence
      {
	symmMatrix densityMatrix_valence(densityMatrix);
	densityMatrix_valence += (ergo_real)-1 * densityMatrix_core;
	symmMatrix G_valence;
	G_valence.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);  
	symmMatrix J_matrix_valence;
	symmMatrix K_matrix_valence;
	symmMatrix Fxc_matrix_valence;
	SCF_statistics curr_cycle_stats_dummy;
	int no_of_valence_electrons = noOfElectrons - scfopts.no_of_core_electrons;
	if(get_2e_matrix_and_energy_sparse(basisInfo, 
					   basisInfoDensFit, 
					   molecule,
					   integralInfo, 
					   G_valence,
					   densityMatrix_valence, 
					   J_K_params,
					   CAM_params,
					   gridParams,
					   scfopts.use_dft,
					   &energy_2el_valence, no_of_valence_electrons,
					   densfit_data,
					   matOpts.size_block_info,
					   matOpts.permutationHML,
					   matOpts.inversePermutationHML,
					   scfopts.do_sparsity_investigation,
					   J_matrix_valence,
					   K_matrix_valence,
					   Fxc_matrix_valence,
					   curr_cycle_stats_dummy) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "Error in get_2e_matrix_and_energy_sparse for densityMatrix_valence.");
	  throw "Error in get_2e_matrix_and_energy_sparse for densityMatrix_valence.";
	}
      }
      // Compute J_matrix_core
      {
	symmMatrix J_matrix_core;
	symmMatrix K_matrix_core;
	symmMatrix Fxc_matrix_core;
	SCF_statistics curr_cycle_stats_dummy;
	energy_2el_core = 0;
	if(get_2e_matrix_and_energy_sparse(basisInfo, 
					   basisInfoDensFit, 
					   molecule,
					   integralInfo, 
					   twoel_matrix_core,
					   densityMatrix_core, 
					   J_K_params,
					   CAM_params,
					   gridParams,
					   scfopts.use_dft,
					   &energy_2el_core, scfopts.no_of_core_electrons,
					   densfit_data,
					   matOpts.size_block_info,
					   matOpts.permutationHML,
					   matOpts.inversePermutationHML,
					   scfopts.do_sparsity_investigation,
					   J_matrix_core,
					   K_matrix_core,
					   Fxc_matrix_core,
					   curr_cycle_stats_dummy) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "Error in get_2e_matrix_and_energy_sparse for densityMatrix_core.");
	  throw "Error in get_2e_matrix_and_energy_sparse for densityMatrix_core.";
	}
      }
    }
  }

  densityMatrix.writeToFile();
  densityMatrix_core.writeToFile();

  // Check that G matrix is free from "inf", "nan" etc.
  if(check_if_matrix_contains_strange_elements(G, matOpts.inversePermutationHML))
    throw std::runtime_error("error in SCF_restricted::get_2e_part_and_energy(): G matrix contains inf or nan.");
  
  // calculate Fock matrix F = H_core + G
  H_core_Matrix.readFromFile();
  FockMatrix.readFromFile();
  FockMatrix = H_core_Matrix + G;
  H_core_Matrix.writeToFile();
  G.clear();

  // Save a copy of FockMatrix before truncation if verification requested.
  symmMatrix FockMatrixBeforeTruncation;
  if(scfopts.do_f_thresh_verification == 1)
    FockMatrixBeforeTruncation = FockMatrix;

  if ( scfopts.do_sparsity_investigation == 1 )
    // Create m-file for dist vs mag
    output_distance_vs_magnitude( basisInfo, FockMatrix, 
				  matOpts.inversePermutationHML,
				  "dist_vs_mag_F_S",
				  scfopts.sparsity_plots_resolution_r,  
				  scfopts.sparsity_plots_resolution_m);
  
  // Do truncation of FockMatrix taking into account gap of Fprev and norm of Z.
  ergo_real gapOfFprevMin = lumoInterval_Fprev.low() - homoInterval_Fprev.upp();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "About to truncate FockMatrix, gap of Fprev >= %22.11f", (double)gapOfFprevMin);
  {
    // Compare FockMatrix to Fprev to check how far eigenvalues can have moved.
    Fprev.readFromFile();
    ergo_real maxEigValMovement_frob = symmMatrix::frob_diff(FockMatrix, Fprev);
    ergo_real acc = std::sqrt(get_machine_epsilon());
    Util::TimeMeter timeMeterMixedDiff;
    ergo_real maxEigValMovement_mixed = symmMatrix::mixed_diff(FockMatrix, Fprev, acc) + acc;
    timeMeterMixedDiff.print(LOG_AREA_SCF, "symmMatrix::mixed_diff for F vs Fprev maxEigValMovement_mixed");
    Util::TimeMeter timeMeterEuclDiff;
    ergo_real maxEigValMovement_eucl = symmMatrix::eucl_diff(FockMatrix, Fprev, acc) + acc;
    timeMeterEuclDiff.print(LOG_AREA_SCF,  "symmMatrix::eucl_diff  for F vs Fprev maxEigValMovement_eucl ");
    Fprev.writeToFile();
    // Increase HOMO/LUMO intervals so that they for sure contain the HOMO and LUMO eigenvalues of F_ort
    intervalType homoInterval = homoInterval_Fprev; 
    intervalType lumoInterval = lumoInterval_Fprev; 
    homoInterval.increase(maxEigValMovement_eucl);
    lumoInterval.increase(maxEigValMovement_eucl);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxEigValMovement_frob  = %22.11f", (double)maxEigValMovement_frob);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxEigValMovement_mixed = %22.11f", (double)maxEigValMovement_mixed);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxEigValMovement_eucl  = %22.11f", (double)maxEigValMovement_eucl);
    // Now we have homoInterval and lumoInterval valid for FockMatrix.
    ergo_real gapMin = lumoInterval.low() - homoInterval.upp();
    ergo_real gapMax = lumoInterval.upp() - homoInterval.low();
    ergo_real threshold_1;
    // Choose subspace error as for purification.
    ergo_real subspaceThr_1     = 0.1 * scfopts.purification_subspace_err_limit;
    // We consider the gap to be accurately known if the uncertainty is at most 10 %
    if ( gapMin > 0 && (gapMax - gapMin) / gapMin < 0.1 ) 
      // Gap is accurately known: we use gapMin
      threshold_1 = subspaceThr_1 * gapMin / (1+subspaceThr_1);
    else 
      // Gap is not accurately known. To avoid choosing a very tight
      // threshold value due to a small lower bound for the gap, we
      // use the largest of 'gap_expected_lower_bound' and calculated
      // 'gapMin':
      threshold_1 = gapMin > scfopts.gap_expected_lower_bound ?
	subspaceThr_1 * gapMin / (1+subspaceThr_1) :
	subspaceThr_1 * scfopts.gap_expected_lower_bound / (1+subspaceThr_1);
    /* Truncate matrix taking into account that we are in
       'non-orthogonal basis', passing invCholFactor to thresh */
    invCholFactor.readFromFile();
    Util::TimeMeter timeMeterEuclThresh;
    int n = basisInfo.noOfBasisFuncs;
    double nnzF_S_before_trunc_pc  = (double)FockMatrix.nnz()  * 100 / ((double)n*n); 
    ergo_real truncError_1 = FockMatrix.eucl_thresh( threshold_1, &invCholFactor );
    double nnzF_S_after_trunc_pc  = (double)FockMatrix.nnz()  * 100 / ((double)n*n); 
    invCholFactor.writeToFile();
    timeMeterEuclThresh.print(LOG_AREA_SCF,  "FockMatrix.eucl_thresh() (with Z)");
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "Truncated FockMatrix (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
	      threshold_1, truncError_1, nnzF_S_before_trunc_pc, nnzF_S_after_trunc_pc);
    this->curr_cycle_stats->add_value("investigation_nnz_percentage_F_S", nnzF_S_after_trunc_pc);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Computing Euclidean norm of FockMatrix...");
    Util::TimeMeter timeMeterEuclF;
    ergo_real FockMatrix_euclnorm = FockMatrix.eucl( std::sqrt(std::numeric_limits<ergo_real>::epsilon()) );
    timeMeterEuclF.print(LOG_AREA_SCF, "FockMatrix.eucl()");
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Euclidean norm of FockMatrix = %22.11f", (double)FockMatrix_euclnorm);

    if ( scfopts.do_sparsity_investigation == 1 ) {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_sparsity_investigation requested, testing truncation of H_core, J, K, and F_xc ");
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "and outputting m-files for dist vs magnitude figs for H_core, J, K, and F_xc.");
      invCholFactor.readFromFile();
      if ( !H_core_Matrix.is_empty() ) {
	symmMatrix H_core_copy(H_core_Matrix);
	H_core_copy.readFromFile();
	// Create m-file for dist vs mag
	output_distance_vs_magnitude( basisInfo, H_core_copy, 
				      matOpts.inversePermutationHML,
				      "dist_vs_mag_Hcore",
				      scfopts.sparsity_plots_resolution_r,  
				      scfopts.sparsity_plots_resolution_m);
	// M-file ready, continue with truncation test
	double nnz_before_trunc_pc  = (double)H_core_copy.nnz()  * 100 / ((double)n*n); 
	ergo_real truncError = H_core_copy.eucl_thresh( threshold_1, &invCholFactor );
	double nnz_after_trunc_pc  = (double)H_core_copy.nnz()  * 100 / ((double)n*n); 
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		  "Truncated H_core (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
		  threshold_1, truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);	
	this->curr_cycle_stats->add_value("investigation_nnz_percentage_H_core", nnz_after_trunc_pc);
      } // end if H_core ...
      else {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "No H_core matrix stored, skipping H_core truncation test.");      
      }
      if ( !J_matrix.is_empty() ) {
	// Create m-file for dist vs mag
	output_distance_vs_magnitude( basisInfo, J_matrix, 
				      matOpts.inversePermutationHML,
				      "dist_vs_mag_J",
				      scfopts.sparsity_plots_resolution_r,  
				      scfopts.sparsity_plots_resolution_m);
	// M-file ready, continue with truncation test
	double nnz_before_trunc_pc  = (double)J_matrix.nnz()  * 100 / ((double)n*n); 
	ergo_real truncError = J_matrix.eucl_thresh( threshold_1, &invCholFactor );
	double nnz_after_trunc_pc  = (double)J_matrix.nnz()  * 100 / ((double)n*n); 
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		  "Truncated J      (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
		  threshold_1, truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);	
	this->curr_cycle_stats->add_value("investigation_nnz_percentage_J", nnz_after_trunc_pc);
      } // end if J ...
      else {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "No J matrix stored, skipping J truncation test.");      
      }
      if ( !K_matrix.is_empty() ) {
	// Create m-file for dist vs mag
	output_distance_vs_magnitude( basisInfo, K_matrix, 
				      matOpts.inversePermutationHML,
				      "dist_vs_mag_K",
				      scfopts.sparsity_plots_resolution_r,  
				      scfopts.sparsity_plots_resolution_m);
	// M-file ready, continue with truncation test
	double nnz_before_trunc_pc  = (double)K_matrix.nnz()  * 100 / ((double)n*n); 
	ergo_real truncError = K_matrix.eucl_thresh( threshold_1, &invCholFactor );
	double nnz_after_trunc_pc  = (double)K_matrix.nnz()  * 100 / ((double)n*n); 
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		  "Truncated K      (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
		  threshold_1, truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);	
	this->curr_cycle_stats->add_value("investigation_nnz_percentage_K", nnz_after_trunc_pc);
      } // end if K ...
      else {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "No K matrix stored, skipping K truncation test.");      
      }
      if ( !Fxc_matrix.is_empty() ) {
	// Create m-file for dist vs mag
	output_distance_vs_magnitude( basisInfo, Fxc_matrix, 
				      matOpts.inversePermutationHML,
				      "dist_vs_mag_Fxc",
				      scfopts.sparsity_plots_resolution_r,  
				      scfopts.sparsity_plots_resolution_m);
	// M-file ready, continue with truncation test
	double nnz_before_trunc_pc  = (double)Fxc_matrix.nnz()  * 100 / ((double)n*n); 
	ergo_real truncError = Fxc_matrix.eucl_thresh( threshold_1, &invCholFactor );
	double nnz_after_trunc_pc  = (double)Fxc_matrix.nnz()  * 100 / ((double)n*n); 
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		  "Truncated Fxc    (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
		  threshold_1, truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);	
	this->curr_cycle_stats->add_value("investigation_nnz_percentage_Fxc", nnz_after_trunc_pc);
      } // end if Fxc ...
      else {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "No Fxc matrix stored, skipping Fxc truncation test.");      
      }
      invCholFactor.writeToFile();
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_sparsity_investigation, testing truncation of H_core, J, K, and F_xc finished.");      
    } // end if ( scfopts.do_sparsity_investigation == 1 )
  }
  //  FockMatrix.frob_thresh(matOpts.sparse_threshold);

  if(scfopts.do_f_thresh_verification == 1) {
    Util::TimeMeter timeMeterThreshVerification;
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_f_thresh_verification requested, computing error matrices...");
    symmMatrix F1(FockMatrixBeforeTruncation);
    symmMatrix F2(FockMatrix);
    ergo_real acc = std::sqrt(get_machine_epsilon());
    // First get norm of error matrix in non-orthogonal basis.
    output_diff_norm_values(F1, F2, acc, "F in in non-orthogonal basis");
    // Now get norm of error matrix in non-orthogonal basis.
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_f_thresh_verification requested, doing Z multiplications...");
    invCholFactor.readFromFile();
    Util::TimeMeter timeMeterZFZ1;
    F1 = transpose(invCholFactor) * F1 * invCholFactor;
    timeMeterZFZ1.print(LOG_AREA_DENSFROMF, "transpose(invCholFactor) * F1 * invCholFactor");
    Util::TimeMeter timeMeterZFZ2;
    F2 = transpose(invCholFactor) * F2 * invCholFactor;
    timeMeterZFZ2.print(LOG_AREA_DENSFROMF, "transpose(invCholFactor) * F2 * invCholFactor");
    invCholFactor.writeToFile();
    // Now we have the matrix in orthogonal basis, before and after truncation.
    output_diff_norm_values(F1, F2, acc, "F in in orthogonal basis");
    timeMeterThreshVerification.print(LOG_AREA_SCF,  "do_f_thresh_verification stuff");
  }

  FockMatrix.writeToFile();
}



void SCF_restricted::output_sparsity_S_F_D(SCF_statistics & stats)
{
  int n = basisInfo.noOfBasisFuncs;
  S_symm.readFromFile();
  output_sparsity_symm(n, S_symm, "S");
  stats.add_value("nnz_S", S_symm.nnz());
  S_symm.writeToFile();
  FockMatrix.readFromFile();
  output_sparsity_symm(n, FockMatrix, "F");
  stats.add_value("nnz_F", FockMatrix.nnz());
  FockMatrix.writeToFile();
  densityMatrix.readFromFile();
  output_sparsity_symm(n, densityMatrix, "D");
  stats.add_value("nnz_D", densityMatrix.nnz());
  densityMatrix.writeToFile();
}


void SCF_restricted::calculate_energy()
{
  // calculate energy
  H_core_Matrix.readFromFile();
  densityMatrix.readFromFile();
  energy = symmMatrix::trace_ab(densityMatrix, H_core_Matrix) + energy_2el;
  if(scfopts.compute_core_density == 1) {
    densityMatrix_core.readFromFile();
    symmMatrix densityMatrix_valence(densityMatrix);
    densityMatrix_valence += (ergo_real)-1 * densityMatrix_core;
    energy_of_valence = symmMatrix::trace_ab(densityMatrix_valence, H_core_Matrix) + symmMatrix::trace_ab(densityMatrix_valence, twoel_matrix_core) + energy_2el_valence;
    energy_reference = symmMatrix::trace_ab(densityMatrix_core, H_core_Matrix) + energy_2el_core + nuclearEnergy;
    densityMatrix_core.writeToFile();
  }
  densityMatrix.writeToFile();
  H_core_Matrix.writeToFile();
  energy += nuclearEnergy;
}


void SCF_restricted::get_FDSminusSDF()
{
  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling compute_FDSminusSDF_sparse, n = %i", n);
  densityMatrix.readFromFile();
  S_symm.readFromFile();
  FockMatrix.readFromFile();
  compute_FDSminusSDF_sparse(n, FockMatrix, densityMatrix, S_symm,
			     ErrorMatrix, matOpts.sparse_threshold);
  S_symm.writeToFile();
  FockMatrix.writeToFile();
  densityMatrix.writeToFile();
  // write to file and read back again to reduce memory fragmentation.
  ErrorMatrix.writeToFile();
  ErrorMatrix.readFromFile();
  output_sparsity(n, ErrorMatrix, "FDS-SDF");
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::get_FDSminusSDF finished.");
}


void SCF_restricted::get_error_measure()
{
  ergo_real error_maxabs = compute_maxabs_sparse(ErrorMatrix);  
  ergo_real error_frob = ErrorMatrix.frob();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxabs FDS-SDF is %8.3g", (double)error_maxabs);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "frob   FDS-SDF is %8.3g", (double)error_frob);
  errorMeasure = error_maxabs;
}


void SCF_restricted::add_to_DIIS_list()
{
  FockMatrix.readFromFile();
  if(((DIISManagerRestricted*)DIIS)->AddIterationToList(FockMatrix, ErrorMatrix) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in DIIS AddIterationToList");
      throw "error in DIIS AddIterationToList";
    }
  FockMatrix.writeToFile();
}


void SCF_restricted::update_best_fock_so_far()
{
  Fprev.readFromFile();
  bestFockMatrixSoFar.readFromFile();
  bestFockMatrixSoFar = Fprev;
  Fprev.writeToFile();
  bestFockMatrixSoFar.writeToFile();
  FockMatrix.readFromFile();
  bestFockMatrixSoFar2.readFromFile();
  bestFockMatrixSoFar2 = FockMatrix;
  FockMatrix.writeToFile();
  bestFockMatrixSoFar2.writeToFile();
}


void SCF_restricted::combine_old_fock_matrices(ergo_real stepLength)
{
  bestFockMatrixSoFar.readFromFile();
  bestFockMatrixSoFar2.readFromFile();
  FockMatrix.readFromFile();
  FockMatrix = 0;
  FockMatrix += stepLength * bestFockMatrixSoFar2;
  FockMatrix += (1 - stepLength) * bestFockMatrixSoFar;
  FockMatrix.writeToFile();
  bestFockMatrixSoFar.writeToFile();
  bestFockMatrixSoFar2.writeToFile();
}


void SCF_restricted::use_diis_to_get_new_fock_matrix()
{
  symmMatrix newFsymm;
  if(((DIISManagerRestricted*)DIIS)->GetCombinedFockMatrix(newFsymm) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in DIIS.GetCombinedFockMatrix");
      throw "error in DIIS.GetCombinedFockMatrix";
    }
  FockMatrix.readFromFile();
  FockMatrix = newFsymm;
  FockMatrix.writeToFile();
}


void SCF_restricted::clear_diis_list()
{
  ((DIISManagerRestricted*)DIIS)->ClearList();
}


void SCF_restricted::clear_error_matrices()
{
  ErrorMatrix.clear();
}


void SCF_restricted::save_current_fock_as_fprev()
{
  // save current Fock matrix as Fprev
  FockMatrix.readFromFile();
  Fprev.readFromFile();
  Fprev = FockMatrix;
  FockMatrix.writeToFile();
  Fprev.writeToFile();
}


void SCF_restricted::get_new_density_matrix()
{
  int n = basisInfo.noOfBasisFuncs;

  if ( scfopts.do_sparsity_investigation_reppuri == 1 ) {
    if(scfopts.shift_using_prev_density_matrix != 0)
      throw "SCF_restricted::get_new_density_matrix: error: do_sparsity_investigation_reppuri cannot be combined with shift_using_prev_density_matrix.";
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "do_sparsity_investigation_reppuri, running purification three times using different norms.");
    mat::normType ntype[3] = {mat::frobNorm, mat::mixedNorm, mat::euclNorm};
    std::map<std::string, double> puri_stats;
    for (int i = 0; i<3; i++) {      
      // Copy input parameters that will be overwritten by purification
      symmMatrix F_ort_prev_copy(F_ort_prev);
      intervalType homoInterval_F_ort_prev_copy(homoInterval_F_ort_prev);
      intervalType lumoInterval_F_ort_prev_copy(lumoInterval_F_ort_prev);
      if(get_dens_from_fock_general(n,
				    noOfElectrons / 2,
				    scfopts.use_diagonalization,
				    scfopts.use_diag_on_error,
				    scfopts.electronic_temperature,
				    densityMatrix,
				    2,
				    electronicEntropyTerm,
				    FockMatrix, // should be same as Fprev now
				    homoInterval_Fprev, // output
				    lumoInterval_Fprev, // output
				    S_symm,
				    invCholFactor,
				    invCholFactor_euclnorm,
				    scfopts.gap_expected_lower_bound,
				    matOpts.size_block_info,
				    F_ort_prev_copy,
				    homoInterval_F_ort_prev_copy,
				    lumoInterval_F_ort_prev_copy,
				    scfopts.purification_eigvalue_err_limit,
				    scfopts.purification_subspace_err_limit,
				    ntype[i],
				    scfopts.purification_maxmul,
				    scfopts.purification_create_m_files,
				    scfopts.purification_ignore_failure,
				    scfopts.purification_use_rand_perturbation_for_alleigsint,
				    scfopts.use_stochastic_orbs,
				    stochastic_orbitals,
				    "investigation_"+mat::getNormTypeString(ntype[i])+"_",
				    puri_stats,
				    scfopts.do_sparsity_investigation,
				    scfopts.sparsity_plots_resolution_m,
				    0,
				    scfopts.do_puri_mmul_tests) != 0)
	{
	  throw "SCF_restricted::get_new_density_matrix: Error in get_dens_from_fock_general";
	}
      
    } // end for       
    // Create m-file for dist vs mag
    densityMatrix.readFromFile();
    output_distance_vs_magnitude( basisInfo, densityMatrix, 
				  matOpts.inversePermutationHML,
				  "dist_vs_mag_D_S",
				  scfopts.sparsity_plots_resolution_r,  
				  scfopts.sparsity_plots_resolution_m);
    densityMatrix.writeToFile();      
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_sparsity_investigation for purification finished, continuing as usual.");    
    this->curr_cycle_stats->add_values(puri_stats);
  } // end if investigation


  // As input to the density matrix construction routine, the default
  // is to use FockMatrix. However, if shift_using_prev_density_matrix
  // we use a modified matrix instead.
  symmMatrix* F_effective = &FockMatrix;
  symmMatrix F_modified;
  if(scfopts.shift_using_prev_density_matrix != 0) {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Using shift_using_prev_density_matrix, shifting by %9.5f a.u. = %9.5f eV", 
	      (double)scfopts.shift_using_prev_density_matrix, (double)scfopts.shift_using_prev_density_matrix / UNIT_one_eV);
    F_modified = FockMatrix;
    F_modified.readFromFile();
    // Get matrix SDS
    symmMatrix SDS(densityMatrix);
    SDS.readFromFile();
    transform_with_S(SDS);
    // Use factor 0.5 since this is restricted case.
    F_modified += (ergo_real)-0.5 * scfopts.shift_using_prev_density_matrix * SDS;
    F_modified.writeToFile();
    F_effective = &F_modified;
  }

  std::map<std::string, double> puri_stats;
  if(get_dens_from_fock_general(n,
				noOfElectrons / 2,
				scfopts.use_diagonalization,
				scfopts.use_diag_on_error,
				scfopts.electronic_temperature,
				densityMatrix,
				2,
				electronicEntropyTerm,
				*F_effective, // should normally be same as Fprev now
				homoInterval_Fprev, // output
				lumoInterval_Fprev, // output
				S_symm,
				invCholFactor,
				invCholFactor_euclnorm,
				scfopts.gap_expected_lower_bound,
				matOpts.size_block_info,
				F_ort_prev,
				homoInterval_F_ort_prev,
				lumoInterval_F_ort_prev,
				scfopts.purification_eigvalue_err_limit,
				scfopts.purification_subspace_err_limit,
				scfopts.purification_truncation_norm,
				scfopts.purification_maxmul,
				scfopts.purification_create_m_files,
				scfopts.purification_ignore_failure,
				scfopts.purification_use_rand_perturbation_for_alleigsint,
				scfopts.use_stochastic_orbs,
				stochastic_orbitals,
				"",
				puri_stats,
				scfopts.do_sparsity_investigation,
				scfopts.sparsity_plots_resolution_m,
				scfopts.do_comparison_to_simple_purification,
				scfopts.do_puri_mmul_tests,
				&eigVecLUMO,
				&eigVecHOMO) != 0)
    {
      throw "SCF_restricted::get_new_density_matrix: Error in get_dens_from_fock_general";
    }
  this->curr_cycle_stats->add_values(puri_stats);

  if(scfopts.compute_core_density == 1) {
    std::map<std::string, double> puri_stats_unused;
    intervalType homoInterval_dummy1(-1e22,1e22);
    intervalType lumoInterval_dummy1(-1e22,1e22);
    intervalType homoInterval_dummy2(-1e22,1e22);
    intervalType lumoInterval_dummy2(-1e22,1e22);
    if(get_dens_from_fock_general(n,
				  scfopts.no_of_core_electrons / 2,
				  scfopts.use_diagonalization,
				  scfopts.use_diag_on_error,
				  scfopts.electronic_temperature,
				  densityMatrix_core,
				  2,
				  electronicEntropyTerm,
				  *F_effective, // should normally be same as Fprev now
				  homoInterval_dummy1, // output
				  lumoInterval_dummy1, // output
				  S_symm,
				  invCholFactor,
				  invCholFactor_euclnorm,
				  scfopts.gap_expected_lower_bound,
				  matOpts.size_block_info,
				  F_ort_prev,
				  homoInterval_dummy2,
				  lumoInterval_dummy2,
				  scfopts.purification_eigvalue_err_limit,
				  scfopts.purification_subspace_err_limit,
				  scfopts.purification_truncation_norm,
				  scfopts.purification_maxmul,
				  scfopts.purification_create_m_files,
				  scfopts.purification_ignore_failure,
				  scfopts.purification_use_rand_perturbation_for_alleigsint,
				  scfopts.use_stochastic_orbs,
				  stochastic_orbitals,
				  "",
				  puri_stats_unused,
				  0,
				  0,
				  0,
				  0,
				  NULL,
				  NULL) != 0)
      {
	throw "SCF_restricted::get_new_density_matrix: Error in get_dens_from_fock_general for core density matrix.";
      }
  }

  // Report sparsity of D and trace(DS)
  S_symm.readFromFile();
  densityMatrix.readFromFile();
  output_sparsity_symm(n, densityMatrix, "new density matrix");
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Tr( D * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix, S_symm));
  densityMatrix.writeToFile();
  S_symm.writeToFile();
}


static int write_matrix_to_file(symmMatrix & M, const std::vector<int> & inversePermutationHML, const BasisInfoStruct & basisInfo, const char* fileName) {
  M.readFromFile();
  matrix_description_struct matrixList[2];
  int nvalues = M.nvalues();
  std::vector<int> rowind; 
  rowind.reserve(nvalues);
  std::vector<int> colind; 
  colind.reserve(nvalues);
  std::vector<ergo_real> values; 
  values.reserve(nvalues);
  M.get_all_values(rowind, colind, values, inversePermutationHML, inversePermutationHML);
  M.writeToFile();
  matrixList[0].nvalues = nvalues;
  matrixList[0].rowind = &rowind[0];
  matrixList[0].colind = &colind[0];
  matrixList[0].values = &values[0];
  if(ddf_writeShellListAndDensityMatricesToFile_sparse(basisInfo, 1, matrixList, fileName) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "Error in ddf_writeShellListAndDensityMatricesToFile_sparse.");
    return -1;
  }
  return 0;
}

void SCF_restricted::write_density_to_file()
{
  if(write_matrix_to_file(densityMatrix, matOpts.inversePermutationHML, basisInfo, "density.bin") != 0)
    throw "error in SCF_restricted::write_density_to_file(): write_matrix_to_file failed for densityMatrix.";
  if(scfopts.compute_core_density == 1) {
    if(write_matrix_to_file(densityMatrix_core, matOpts.inversePermutationHML, basisInfo, "density_core.bin") != 0)
      throw "error in SCF_restricted::write_density_to_file(): write_matrix_to_file failed for densityMatrix_core.";
  }
}


void SCF_restricted::save_final_potential()
{
  if(save_symmetric_matrix(FockMatrix, basisInfo, "potential.bin",
			   matOpts.inversePermutationHML) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF,
              "error in ddf_writeShellListAndDensityMatricesToFile");
    throw "error in ddf_writeShellListAndDensityMatricesToFile";
  }
}


void SCF_restricted::save_full_matrices_for_matlab()
{
  int n = basisInfo.noOfBasisFuncs;

  FockMatrix.readFromFile();
  write_full_matrix(n, FockMatrix, "matrix_F",
		    matOpts.inversePermutationHML);
  FockMatrix.writeToFile();

  S_symm.readFromFile();
  write_full_matrix(n, S_symm, "matrix_S",
		    matOpts.inversePermutationHML);
  S_symm.writeToFile();

  densityMatrix.readFromFile();
  write_full_matrix(n, densityMatrix, "matrix_D",
		    matOpts.inversePermutationHML);
  densityMatrix.writeToFile();
}


void SCF_restricted::output_density_images()
{
  Util::TimeMeter timeMeter;

  int n = basisInfo.noOfBasisFuncs;

  ergo_real* densityMatrixFull_tot  = new ergo_real[n*n];
  ergo_real* densityMatrixFull_spin = new ergo_real[n*n];

  // Get full matrix version of density matrix, and empty spin density matrix.
  {
    std::vector<ergo_real> densityMatrixFull(n*n);
    
    densityMatrix.readFromFile();
    densityMatrix.fullMatrix(densityMatrixFull,
			     matOpts.inversePermutationHML,
			     matOpts.inversePermutationHML);
    densityMatrix.writeToFile();
        
    for(int i = 0; i < n*n; i++)
      {
	densityMatrixFull_tot [i] = densityMatrixFull[i];
	densityMatrixFull_spin[i] = 0;
      }
  }

  do_density_images(basisInfo,
		    molecule,
		    densityMatrixFull_tot, 
		    densityMatrixFull_spin,
		    scfopts.output_density_images_boxwidth);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::output_density_images finished OK.");
  timeMeter.print(LOG_AREA_SCF, "SCF_restricted::output_density_images");
}


static ergo_real rand_minus1_to_1(int rand_param) {
  if(rand_param == 0) {
    // ELIAS NOTE 2013-11-14: this is the old bad implementation that gives a slightly "unfair" distribution in the [-1,1] interval.
    ergo_real randomNumber = (ergo_real)rand() / RAND_MAX;
    // Now randomNumber is between 0 and 1
    randomNumber *= 2;
    // Now randomNumber is between 0 and 2
    randomNumber -= 1;
    // Now randomNumber is between -1 and 1
    return randomNumber;
  }
  else if(rand_param == 1) {
    ergo_real h = (ergo_real)1.0 / ((ergo_real)RAND_MAX+1);
    ergo_real randomNumber = (ergo_real)rand() / ((ergo_real)RAND_MAX+1);
    // Now randomNumber is between 0 and 1-h, and the spacing between possible values is h.
    randomNumber *= 2;
    // Now randomNumber is between 0 and 2-2*h, and the spacing between possible values is 2*h.
    randomNumber -= 1;
    // Now randomNumber is between -1 and 1-2*h, and the spacing between possible values is 2*h.
    randomNumber += h;
    // Now randomNumber is between -1+h and 1-h, and the spacing between possible values is 2*h.
    return randomNumber;
  }
  else if(rand_param == 2) {
    ergo_real randomNumber_a = ((ergo_real)RAND_MAX+1)*(ergo_real)rand() + (ergo_real)rand();
    ergo_real randomNumber_a_max = ((ergo_real)RAND_MAX+1)*((ergo_real)RAND_MAX+1) - 1;
    // Now randomNumber_a is between 0 and randomNumber_a_max, with spacing 1 between possible values.
    ergo_real h = (ergo_real)1.0 / (randomNumber_a_max+1);
    ergo_real randomNumber = randomNumber_a / (randomNumber_a_max+1);
    // Now randomNumber is between 0 and 1-h, and the spacing between possible values is h.
    randomNumber *= 2;
    // Now randomNumber is between 0 and 2-2*h, and the spacing between possible values is 2*h.
    randomNumber -= 1;
    // Now randomNumber is between -1 and 1-2*h, and the spacing between possible values is 2*h.
    randomNumber += h;
    // Now randomNumber is between -1+h and 1-h, and the spacing between possible values is 2*h.
    return randomNumber;
  }
  else
    throw "Error in rand_minus1_to_1: unsupported rand_param value.";
}


static ergo_real get_random_vector_with_norm_leq_one_recursive(std::vector<ergo_real> & v, int n, int rand_param) {
  if(n == 1) {
    v.resize(1);
    v[0] = rand_minus1_to_1(rand_param);
    ergo_real norm = std::sqrt(v[0]*v[0]);
    return norm;
  }
  while(1) {
    // Consider two parts.
    int n1 = n - n/2;
    int n2 = n - n1;
    std::vector<ergo_real> v1;
    std::vector<ergo_real> v2;
    ergo_real norm1 = get_random_vector_with_norm_leq_one_recursive(v1, n1, rand_param);
    ergo_real norm2 = get_random_vector_with_norm_leq_one_recursive(v2, n2, rand_param);
    ergo_real norm = std::sqrt(norm1*norm1+norm2*norm2);
    if(norm <= 1) {
      v.resize(n);
      for(int i = 0; i < n1; i++)
	v[i] = v1[i];
      for(int i = 0; i < n2; i++)
	v[n1+i] = v2[i];
      // Done!
      return norm;
    }
  }
}

static void get_normalized_random_vector_2(std::vector<ergo_real> & v, int n, int rand_param) {
  ergo_real norm = get_random_vector_with_norm_leq_one_recursive(v, n, rand_param);
  for(int i = 0; i < n; i++)
    v[i] /= norm;
}

// ELIAS NOTE 2013-11-13: this routine seems wrong, does not seem to give a "fair" distribution of vectors.
static void get_normalized_random_vector_recursive(std::vector<ergo_real> & v, int n, int rand_param) {
  if(n == 1) {
    v.resize(1);
    if(rand() % 2 == 0)
      v[0] = 1;
    else
      v[0] = -1;
    return;
  }
  // Consider two parts.
  int n1 = n - n/2;
  int n2 = n - n1;
  std::vector<ergo_real> v1;
  std::vector<ergo_real> v2;
  get_normalized_random_vector_recursive(v1, n1, rand_param);
  get_normalized_random_vector_recursive(v2, n2, rand_param);
  ergo_real x1 = 0;
  ergo_real x2 = 0;
  while(1) {
    x1 = rand_minus1_to_1(rand_param);
    x2 = rand_minus1_to_1(rand_param);
    if(x1*x1+x2*x2 <= 1)
      break;
  }
  // OK, now we have x1 and x2 selected randomly in a "fair" way.
  ergo_real norm = std::sqrt(x1*x1+x2*x2);
  x1 /= norm;
  x2 /= norm;
  // OK, now the coordinate pair (x1,x2) points in a random direction in the 2D space.
  v.resize(n);
  for(int i = 0; i < n1; i++)
    v[i] = x1 * v1[i];
  for(int i = 0; i < n2; i++)
    v[n1+i] = x2 * v2[i];
  // Done!
}

/* ELIAS NOTE 2014-07-13: Commented out this unused routine to silence compiler warning.
static void get_normalized_random_vector_naive(std::vector<ergo_real> & v, int n, int rand_param) {
  v.resize(n);
  ergo_real sqSum = 0;
  while(1) {
    sqSum = 0;
    for(int k = 0; k < n; k++) {
      ergo_real x = rand_minus1_to_1(rand_param);
      v[k] = x;
      sqSum += x*x;
    }
    if(sqSum <= 1)
      break;
  } // end while
  // Normalize vector
  ergo_real norm = std::sqrt(sqSum);
  for(int k = 0; k < n; k++)
    v[k] /= norm;
}
*/

/* ELIAS NOTE 2014-07-13: Commented out this unused routine to silence compiler warning.
static void do_random_vector_tests(int rand_param) {
  int N = 200000;
  int n = 3;
  ergo_real limit = 0.95;
  ergo_real probability_naive_min = 1;
  ergo_real probability_naive_max = 0;
  ergo_real probability_recursive_min = 1;
  ergo_real probability_recursive_max = 0;
  ergo_real probability_new_min = 1;
  ergo_real probability_new_max = 0;
  for(int outerLoop = 0; outerLoop < 20; outerLoop++) {
    int first_component_larger_than_limit_count_naive = 0;
    int first_component_larger_than_limit_count_recursive = 0;
    int first_component_larger_than_limit_count_new = 0;
    for(int i = 0; i < N; i++) {
      std::vector<ergo_real> v_naive;
      get_normalized_random_vector_naive(v_naive, n, rand_param);
      if(v_naive[0] > limit)
	first_component_larger_than_limit_count_naive++;
      std::vector<ergo_real> v_recursive;
      get_normalized_random_vector_recursive(v_recursive, n, rand_param);
      if(v_recursive[0] > limit)
	first_component_larger_than_limit_count_recursive++;
      std::vector<ergo_real> v_new;
      get_normalized_random_vector_2(v_new, n, rand_param);
      if(v_new[0] > limit)
	first_component_larger_than_limit_count_new++;
    }
    ergo_real probability_naive = (ergo_real)first_component_larger_than_limit_count_naive / N;
    ergo_real probability_recursive = (ergo_real)first_component_larger_than_limit_count_recursive / N;
    ergo_real probability_new = (ergo_real)first_component_larger_than_limit_count_new / N;
    if(probability_naive > probability_naive_max)
      probability_naive_max = probability_naive;
    if(probability_naive < probability_naive_min)
      probability_naive_min = probability_naive;
    if(probability_recursive > probability_recursive_max)
      probability_recursive_max = probability_recursive;
    if(probability_recursive < probability_recursive_min)
      probability_recursive_min = probability_recursive;
    if(probability_new > probability_new_max)
      probability_new_max = probability_new;
    if(probability_new < probability_new_min)
      probability_new_min = probability_new;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_random_vector_tests(), probability_naive     min max = %12.6f %12.6f.", probability_naive_min,     probability_naive_max    );
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_random_vector_tests(), probability_recursive min max = %12.6f %12.6f.", probability_recursive_min, probability_recursive_max);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_random_vector_tests(), probability_new       min max = %12.6f %12.6f.", probability_new_min,       probability_new_max      );
}
*/

void SCF_restricted::prepare_stochastic_orbitals() {
  //  do_random_vector_tests();
  Util::TimeMeter timeMeter;
  int n = basisInfo.noOfBasisFuncs;
  int nOrbs = scfopts.stochastic_orbs_no_of_vectors;
  stochastic_orbitals.resize(nOrbs);
  if(scfopts.stochastic_orbs_use_unit_vectors == 1) {
    // In this case we simply set the "stochastic orbitals" to unit vectors, and continue in a round-robin way if nOrbs>n.
    for(int i = 0; i < nOrbs; i++) {
      stochastic_orbitals[i].resize(n);
      for(int k = 0; k < n; k++)
	stochastic_orbitals[i][k] = 0;
      stochastic_orbitals[i][i%n] = 1;
    }
  }
  else {
    for(int i = 0; i < nOrbs; i++) {
#if 1
      get_normalized_random_vector_2(stochastic_orbitals[i], n, scfopts.stochastic_orbs_rand_param);
#else
      get_normalized_random_vector_naive(stochastic_orbitals[i], n);
#endif
    }
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::prepare_stochastic_orbitals finished OK.");
  timeMeter.print(LOG_AREA_SCF, "SCF_restricted::prepare_stochastic_orbitals");
}


void SCF_restricted::output_csr_matrices_for_gao()
{
  int n = basisInfo.noOfBasisFuncs;

  FockMatrix.readFromFile();
  write_CSR_matrix(n,
		   FockMatrix,
		   "CSR_Fock_matrix.bin",
		   1e-5,
		   matOpts.inversePermutationHML);
  FockMatrix.writeToFile();

  S_symm.readFromFile();
  write_CSR_matrix(n,
		   S_symm,
		   "CSR_overlap_matrix.bin",
		   1e-5,
		   matOpts.inversePermutationHML);
  S_symm.writeToFile();
}


void SCF_restricted::write_diag_dens_to_file()
{
  int n = basisInfo.noOfBasisFuncs;

  densityMatrix.readFromFile();
  write_diag_elements_to_file(n, densityMatrix, "diagdens.txt",
			      matOpts.permutationHML);
  densityMatrix.writeToFile();
}

void SCF_restricted::report_final_results()
{

}


void SCF_restricted::do_spin_flip(int atomCount)
{
  throw "error: SCF_restricted::do_spin_flip does not make sense, should not have been called.";
}


void SCF_restricted::save_density_as_prevdens()
{
  densityMatrix.readFromFile();
  Dprev.readFromFile();
  Dprev = densityMatrix;
  densityMatrix.writeToFile();
  Dprev.writeToFile();
}


void SCF_restricted::report_density_difference()
{
  if(scfopts.do_report_density_diff == 0) {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "SCF_restricted::report_density_difference() skipping: (scfopts.do_report_density_diff == 0).");
    return;
  }
  Util::TimeMeter tm;
  densityMatrix.readFromFile();
  Dprev.readFromFile();
  symmMatrix diff(densityMatrix);
  diff += (ergo_real)-1.0 * Dprev;
  ergo_real diff_eucl = GetEuclideanNormOfMatrix(diff);
  densityMatrix.writeToFile();
  Dprev.writeToFile();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "SCF_restricted::report_density_difference, diff_eucl = %22.11f", (double)diff_eucl);
  tm.print(LOG_AREA_SCF, "SCF_restricted::report_density_difference");
}

void SCF_restricted::compute_dipole_moment()
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::compute_dipole_moment");
  densityMatrix.readFromFile();  
  get_dipole_moment(densityMatrix, basisInfo, matOpts.size_block_info, matOpts.permutationHML, molecule);
  densityMatrix.writeToFile();
}

void SCF_restricted::do_mulliken_pop_stuff()
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::do_mulliken_pop_stuff");
  densityMatrix.readFromFile();  
  S_symm.readFromFile();  
  do_mulliken_atomic_charges(densityMatrix,
			     S_symm,
			     basisInfo,
			     matOpts.size_block_info,
			     matOpts.permutationHML,
			     matOpts.inversePermutationHML,
			     molecule);
  densityMatrix.writeToFile();
  S_symm.writeToFile();
}

void SCF_restricted::create_mtx_files_F(int const scfIter) 
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx file for Fock matrix");
  std::stringstream ss_fileName;
  ss_fileName << "F_matrix_" << scfIter;
  std::stringstream ss_id;
  ss_id << scfopts.calculation_identifier << " - effective Hamiltonian matrix, SCF cycle " << scfIter;  
  FockMatrix.readFromFile();
  write_matrix_in_matrix_market_format( FockMatrix, matOpts.inversePermutationHML, ss_fileName.str(), 
					ss_id.str(), scfopts.method_and_basis_set );
  FockMatrix.writeToFile();
}
void SCF_restricted::create_mtx_files_D(int const scfIter) 
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx file for density matrix");
  std::stringstream ss_fileName;
  ss_fileName << "D_matrix_" << scfIter;
  std::stringstream ss_id;
  ss_id << scfopts.calculation_identifier << " - density matrix, SCF cycle " << scfIter;
  densityMatrix.readFromFile();
  write_matrix_in_matrix_market_format( densityMatrix, matOpts.inversePermutationHML, ss_fileName.str(), 
					ss_id.str(), scfopts.method_and_basis_set );
  densityMatrix.writeToFile();
}

void SCF_restricted::create_homo_eigvec_file() const
{
  if (eigVecHOMO.is_empty()) {
    do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Failed to output HOMO eigenvector to file. No HOMO eigenvector stored.");
    return;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Storing HOMO eigenvector to file homo_coefficient_vec.txt.");
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> homo_vec_perm(n);
  eigVecHOMO.fullvector(homo_vec_perm); 
  // now we have the permuted vector
  std::vector<ergo_real> homo_vec(n);
  for (int ind = 0; ind < n; ind++) 
    homo_vec[ind] = homo_vec_perm[matOpts.inversePermutationHML[ind]];
  char ffname[888];
  sprintf(ffname, "homo_coefficient_vec.txt");
  std::ofstream ff(ffname);
  for (int ind = 0; ind < n; ind++) 
    ff << homo_vec[ind] << std::endl;
  ff.close();  
}

void SCF_restricted::create_lumo_eigvec_file() const
{
  if (eigVecLUMO.is_empty()) {
    do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Failed to output LUMO eigenvector to file. No LUMO eigenvector stored.");
    return;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Storing LUMO eigenvector to file lumo_coefficient_vec.txt.");
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> lumo_vec_perm(n);
  eigVecLUMO.fullvector(lumo_vec_perm); 
  // now we have the permuted vector
  std::vector<ergo_real> lumo_vec(n);
  for (int ind = 0; ind < n; ind++) 
    lumo_vec[ind] = lumo_vec_perm[matOpts.inversePermutationHML[ind]];
  char ffname[888];
  sprintf(ffname, "lumo_coefficient_vec.txt");
  std::ofstream ff(ffname);
  for (int ind = 0; ind < n; ind++) 
    ff << lumo_vec[ind] << std::endl;
  ff.close();  
}


static void 
output_orbital_coeffs_in_gabedit_order(const BasisInfoStruct & basisInfo,
				       std::vector<int> const & shellIdxList,
				       std::ofstream & ff,
				       std::vector<ergo_real> const & orbital_vec) {
  int count = 0;
  for(int i = 0; i < basisInfo.noOfShells; i++) {
    int k = shellIdxList[i];
    int startIdx = basisInfo.shellList[k].startIndexInMatrix;
    switch(basisInfo.shellList[k].shellType) {
    case 0: // s-type shell
      ff << count+1 << "   " << orbital_vec[startIdx] << std::endl; count++;
      break;
    case 1: // p-type shell
      ff << count+1 << "   " << orbital_vec[startIdx+2] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+0] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+1] << std::endl; count++;
      break;
    case 2: // d-type shell
      ff << count+1 << "   " << orbital_vec[startIdx+2] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+3] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+1] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+4] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+0] << std::endl; count++;
      break;
    case 3: // f-type shell
      ff << count+1 << "   " << orbital_vec[startIdx+3] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+4] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+2] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+5] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+1] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+6] << std::endl; count++;
      ff << count+1 << "   " << orbital_vec[startIdx+0] << std::endl; count++;
      break;
    default: throw "error in output_orbital_coeffs_in_gabedit_order: shell types beyond f not implemented!";
    }
  }
  if(count != basisInfo.noOfBasisFuncs)
    throw "error in output_orbital_coeffs_in_gabedit_order: (count != basisInfo.noOfBasisFuncs)";
}

void SCF_restricted::create_gabedit_file() const
{
  if (eigVecHOMO.is_empty() || eigVecLUMO.is_empty()) {
    do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Failed to output HOMO/LUMO eigenvectors to gabedit file; no HOMO/LUMO info available.");
    return;
  }
  if(basisInfo.use_6_d_funcs == 1) {
    do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Failed to output HOMO/LUMO eigenvectors to gabedit file; not implemented for use_6_d_funcs case.");
    return;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating Gabedit file with HOMO/LUMO eigenvector info.");
  int n = basisInfo.noOfBasisFuncs;
  // Get HOMO
  std::vector<ergo_real> homo_vec_perm(n);
  eigVecHOMO.fullvector(homo_vec_perm); 
  // now we have the permuted vector
  std::vector<ergo_real> homo_vec(n);
  for (int ind = 0; ind < n; ind++) 
    homo_vec[ind] = homo_vec_perm[matOpts.permutationHML[ind]];
  //homo_vec[ind] = homo_vec_perm[matOpts.inversePermutationHML[ind]];
  // Get LUMO
  std::vector<ergo_real> lumo_vec_perm(n);
  eigVecLUMO.fullvector(lumo_vec_perm); 
  // now we have the permuted vector
  std::vector<ergo_real> lumo_vec(n);
  for (int ind = 0; ind < n; ind++) 
    lumo_vec[ind] = lumo_vec_perm[matOpts.permutationHML[ind]];
  //lumo_vec[ind] = lumo_vec_perm[matOpts.inversePermutationHML[ind]];
  // Create Gabedit file.
  const char fileName [] = "gabeditfile.gab";
  std::ofstream ff(fileName);
  /* FIXME: check if we should use "Cart" or "Sphe" here. That is,
     should we use use_6_d_funcs? */
  int use_6_d_funcs = 0;
  ff << "[Gabedit Format] Sphe" << std::endl;
  ff << "[Atoms] Angs" << std::endl;
  for(int i = 0; i < molecule.getNoOfAtoms(); i++) {
    char atomLabelString[4];
    get_atom_label_from_charge_int(molecule.getAtom(i).charge, atomLabelString, 4);
    ff << atomLabelString << " " << i+1 << " " << molecule.getAtom(i).charge
       << "   " << molecule.getAtom(i).coords[0] / UNIT_one_Angstrom
       << "   " << molecule.getAtom(i).coords[1] / UNIT_one_Angstrom
       << "   " << molecule.getAtom(i).coords[2] / UNIT_one_Angstrom 
       << std::endl;
  }
  std::vector<int> shellIdxList(basisInfo.noOfShells);
  int shellIdxCounter = 0;
  SquareFuncIntegrator sfi;
  ff << "[Basis]" << std::endl;
  for(int i = 0; i < molecule.getNoOfAtoms(); i++) {
    ff << i+1 << " 0" << std::endl;
    // Now output info about shells for this atom.
    for(int k = 0; k < basisInfo.noOfShells; k++) {
      // Check if this shell belongs to the current atom.
      ergo_real absdx = std::abs(basisInfo.shellList[k].centerCoords[0] - molecule.getAtom(i).coords[0]);
      ergo_real absdy = std::abs(basisInfo.shellList[k].centerCoords[1] - molecule.getAtom(i).coords[1]);
      ergo_real absdz = std::abs(basisInfo.shellList[k].centerCoords[2] - molecule.getAtom(i).coords[2]);
      ergo_real distlimit = 0.01;
      if(absdx > distlimit || absdy > distlimit || absdz > distlimit)
	continue;
      // OK, now we know this shell is at least very near the current atom.
      shellIdxList[shellIdxCounter] = k;
      shellIdxCounter++;
      char shellChar = 'x';
      int shellType = basisInfo.shellList[k].shellType;
      switch(shellType) {
      case 0: shellChar = 's'; break;
      case 1: shellChar = 'p'; break;
      case 2: shellChar = 'd'; break;
      case 3: shellChar = 'f'; break;
      default: throw "SCF_restricted::create_gabedit_file error: shell types beyond f not implemented!";
      }
      ff << shellChar << " " << basisInfo.shellList[k].noOfContr << " 1.00" << std::endl;
      for(int contridx = 0; contridx < basisInfo.shellList[k].noOfContr; contridx++) {
	ergo_real exponent = basisInfo.shellList[k].exponentList[contridx];
	ergo_real shellFactor = sfi.getShellFactor(integralInfo, exponent, shellType, use_6_d_funcs);
	ergo_real coeff = basisInfo.shellList[k].coeffList[contridx] / shellFactor;
	ff << exponent << "   " << coeff << std::endl;
      }
    }
    // Blank line before shells for next atom.
    ff << std::endl;
  }
  if(shellIdxCounter != basisInfo.noOfShells)
    throw "Error: (shellIdxCounter != basisInfo.noOfShells)";
  // MO section.
  ff << "[MO]" << std::endl;
  // HOMO
  ff << "Spin=Alpha" << std::endl;
  ff << "Occup=   2.000000" << std::endl;
  output_orbital_coeffs_in_gabedit_order(basisInfo, shellIdxList, ff, homo_vec);
  // LUMO
  ff << "Spin=Alpha" << std::endl;
  ff << "Occup=   0.000000" << std::endl;
  output_orbital_coeffs_in_gabedit_order(basisInfo, shellIdxList, ff, lumo_vec);
  // Blank line before end of file.
  ff << std::endl;
  // Close file.
  ff.close();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Gabedit file '%s' with HOMO/LUMO eigenvector info created OK.", fileName);
}


void SCF_restricted::update_subspace_diff()
{
  densityMatrix.readFromFile();
  Dprev.readFromFile();
  ergo_real acc = std::sqrt(get_machine_epsilon());

  symmMatrix diff(densityMatrix);
  diff += (ergo_real)-1.0 * Dprev;

  transform_with_S(diff);
  transform_with_invChol(diff);

  // Compensate for factor 2 (restricted case)
  diff *= (ergo_real)0.5;

  ergo_real diff_eucl = diff.eucl(acc);

  densityMatrix.writeToFile();
  Dprev.writeToFile();

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::update_subspace_diff, diff_eucl = %22.11f", (double)diff_eucl);
  curr_subspace_diff = diff_eucl;

}



struct RandomNumber {
  ergo_real accumulate(ergo_real& a, int const dummy1, int const dummy2) {
    a = rand() / (ergo_real)RAND_MAX;
    return 0;
  }
};


/** Transform matrix A to S*A*S */
void SCF_restricted::transform_with_S(symmMatrix & A)
{
  S_symm.readFromFile();

  normalMatrix S_norm(S_symm);
  normalMatrix A_norm(A);

  normalMatrix SA(S_symm);
  SA = (ergo_real)1.0 * S_norm * A_norm;
  normalMatrix SAS(S_symm);
  SAS = (ergo_real)1.0 * SA * S_norm;

  A = SAS;

  S_symm.writeToFile();
}


/** Transform matrix A to invCholT*A*invChol */
void SCF_restricted::transform_with_invChol(symmMatrix & A)
{
  invCholFactor.readFromFile();
  A = transpose(invCholFactor) * A * invCholFactor;
  invCholFactor.writeToFile();
}


void SCF_restricted::get_non_ort_err_mat_normalized_in_ort_basis(symmMatrix & randomMatrix, int transform_with_S_also)
{
  symmMatrix randomMatrix1;
  randomMatrix1.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  symmMatrix randomMatrix2;
  randomMatrix2.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  randomMatrix1.random();
  randomMatrix2.random();
  randomMatrix = 0;
  randomMatrix += (ergo_real) 1.0 * randomMatrix1;
  randomMatrix += (ergo_real)-1.0 * randomMatrix2;
  
  symmMatrix randomMatrix_ort(randomMatrix);
  
  if(transform_with_S_also)
    transform_with_S(randomMatrix_ort);

  transform_with_invChol(randomMatrix_ort);
  
  ergo_real acc = std::sqrt(get_machine_epsilon());
  ergo_real randomMatrix_Norm = randomMatrix.eucl(acc);
  ergo_real randomMatrix_ort_Norm = randomMatrix_ort.eucl(acc);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "norms of randomMatrix and randomMatrix_ort : %22.11f %22.11f", 
	    (double)randomMatrix_Norm, (double)randomMatrix_ort_Norm);

  // Normalize randomMatrix so that randomMatrix_ort would have norm 1.
  randomMatrix *= (ergo_real)(1.0 / randomMatrix_ort_Norm);

  ergo_real randomMatrixNormAfterNormalization = randomMatrix.eucl(acc);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "randomMatrixNormAfterNormalization = %22.11f", (double)randomMatrixNormAfterNormalization);
};



void SCF_restricted::disturb_dens_matrix(ergo_real subspaceError)
{
  ergo_real gap = 2;
  ergo_real desiredErrorNorm = gap * subspaceError / (1 + subspaceError);

  symmMatrix randomMatrix;
  randomMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  get_non_ort_err_mat_normalized_in_ort_basis(randomMatrix, 1);

  densityMatrix.readFromFile();

  densityMatrix += desiredErrorNorm * randomMatrix;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Disturbed density matrix with desiredErrorNorm = %22.11f  (gap = %8.3f)", 
	    (double)desiredErrorNorm, (double)gap);

  densityMatrix.writeToFile();  
}



void SCF_restricted::disturb_dens_matrix_exact_try(const symmMatrix & randomMatrix,
						   const symmMatrix & orgDensMatrix,
						   ergo_real disturbanceFactor,
						   ergo_real & resultSinTheta,
						   symmMatrix & resultDensMatrix)
{
  symmMatrix D(orgDensMatrix);
  D += (ergo_real)1.0 * disturbanceFactor * randomMatrix;

  symmMatrix SDS_symm(D);
  transform_with_S(SDS_symm);
  SDS_symm *= (ergo_real)-1.0;

  symmMatrix F_ort_prev_dummy;
  F_ort_prev_dummy.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  F_ort_prev_dummy.writeToFile();
  intervalType homoInterval_tmp1(-1e22,1e22);
  intervalType lumoInterval_tmp1(-1e22,1e22);
  intervalType homoInterval_tmp2(-1e22,1e22);
  intervalType lumoInterval_tmp2(-1e22,1e22);

  resultDensMatrix.writeToFile();
  SDS_symm.writeToFile();
  int use_diag = 0;
  int use_diag_on_error = 0;
  int n = basisInfo.noOfBasisFuncs;
  std::map<std::string, double> puri_stats;
  ergo_real electronicEntropyTermDummy = 0;
  if(get_dens_from_fock_general(n,
				noOfElectrons / 2,
				use_diag,
				use_diag_on_error,
				scfopts.electronic_temperature,
				resultDensMatrix,
				2,
				electronicEntropyTermDummy,
				SDS_symm,
				homoInterval_tmp1,
				lumoInterval_tmp1,
				S_symm,
				invCholFactor,
				invCholFactor_euclnorm,
				scfopts.gap_expected_lower_bound,
				matOpts.size_block_info,
				F_ort_prev_dummy,
				homoInterval_tmp2,
				lumoInterval_tmp2,
				scfopts.purification_eigvalue_err_limit,
				scfopts.purification_subspace_err_limit,
				scfopts.purification_truncation_norm,
				scfopts.purification_maxmul,
				scfopts.purification_create_m_files,
				scfopts.purification_ignore_failure,
				scfopts.purification_use_rand_perturbation_for_alleigsint,
				scfopts.use_stochastic_orbs,
				stochastic_orbitals,
				"",
				puri_stats,
				scfopts.do_sparsity_investigation,
				scfopts.sparsity_plots_resolution_m,
				scfopts.do_comparison_to_simple_purification,
				scfopts.do_puri_mmul_tests) != 0)
    {
      throw "SCF_restricted::get_starting_guess_density: Error in get_dens_from_fock_general";
    }

  // OK, now we have computed D_Pure which is the purified version of SDS_symm.
  // But D_Pure is not in orthogonal basis.

  resultDensMatrix.readFromFile();

  symmMatrix diff(resultDensMatrix);
  diff += (ergo_real)-1.0 * orgDensMatrix;
  
  transform_with_S(diff);
  transform_with_invChol(diff);

  // Compensate for factor 2 (restricted case)
  diff *= (ergo_real)0.5;

  ergo_real acc = std::sqrt(get_machine_epsilon());

  resultSinTheta = diff.eucl(acc);
}


void SCF_restricted::disturb_dens_matrix_exact(ergo_real subspaceError)
{
  //ergo_real gap = 2;
  //ergo_real desiredErrorNorm = gap * subspaceError / (1 + subspaceError);

  symmMatrix randomMatrix;
  randomMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  get_non_ort_err_mat_normalized_in_ort_basis(randomMatrix, 1);

  symmMatrix newDensMatrix;
  newDensMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);

  densityMatrix.readFromFile();

  ergo_real currSinTheta;
  ergo_real disturbanceFactor_min = 0;
  ergo_real disturbanceFactor_max = 5;

  int iterCount = 0;
  do {
    iterCount++;
    if(iterCount > 44)
      throw "error in SCF_restricted::disturb_dens_matrix_exact, iterCount esceeded limit.";
    ergo_real disturbanceFactor = 0.5 * (disturbanceFactor_min + disturbanceFactor_max);
    disturb_dens_matrix_exact_try(randomMatrix,
				  densityMatrix,
				  disturbanceFactor,
				  currSinTheta,
				  newDensMatrix);
    if(currSinTheta < subspaceError)
      disturbanceFactor_min = disturbanceFactor;
    else
      disturbanceFactor_max = disturbanceFactor;
  } while(std::fabs(currSinTheta - subspaceError) > 0.001*subspaceError);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::disturb_dens_matrix_exact done, iterCount = %2i", iterCount);

  densityMatrix = newDensMatrix;

  densityMatrix.writeToFile();  
}



void SCF_restricted::disturb_fock_matrix(ergo_real subspaceError)
{

  symmMatrix F_ort_prev_dummy;
  F_ort_prev_dummy.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  F_ort_prev_dummy.writeToFile();
  intervalType homoInterval_tmp1(-1e22,1e22);
  intervalType lumoInterval_tmp1(-1e22,1e22);
  intervalType homoInterval_tmp2(-1e22,1e22);
  intervalType lumoInterval_tmp2(-1e22,1e22);

  symmMatrix densityMatrix_dummy;
  densityMatrix_dummy.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  densityMatrix_dummy.writeToFile();
  int use_diag = 0;
  int use_diag_on_error = 0;
  int n = basisInfo.noOfBasisFuncs;

  std::map<std::string, double> puri_stats;
  if(get_dens_from_fock_general(n,
				noOfElectrons / 2,
				use_diag,
				use_diag_on_error,
				scfopts.electronic_temperature,
				densityMatrix_dummy,
				2,
				electronicEntropyTerm,
				FockMatrix,
				homoInterval_tmp1,
				lumoInterval_tmp1,
				S_symm,
				invCholFactor,
				invCholFactor_euclnorm,
				scfopts.gap_expected_lower_bound,
				matOpts.size_block_info,
				F_ort_prev_dummy,
				homoInterval_tmp2,
				lumoInterval_tmp2,
				scfopts.purification_eigvalue_err_limit,
				scfopts.purification_subspace_err_limit,
				scfopts.purification_truncation_norm,
				scfopts.purification_maxmul,
				scfopts.purification_create_m_files,
				scfopts.purification_ignore_failure,
				scfopts.purification_use_rand_perturbation_for_alleigsint,
				scfopts.use_stochastic_orbs,
				stochastic_orbitals,
				"",
				puri_stats,
				scfopts.do_sparsity_investigation,
				scfopts.sparsity_plots_resolution_m,
				scfopts.do_comparison_to_simple_purification,
				scfopts.do_puri_mmul_tests) != 0)
    {
      throw "SCF_restricted::get_starting_guess_density: Error in get_dens_from_fock_general";
    }

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_restricted::disturb_fock_matrix, interval sizes: %22.11f %22.11f",
	    (double)(homoInterval_tmp2.upp() - homoInterval_tmp2.low()), (double)(lumoInterval_tmp2.upp() - lumoInterval_tmp2.low()));
  
  ergo_real gap = lumoInterval_tmp2.low() - homoInterval_tmp2.upp();

  ergo_real desiredErrorNorm = gap * subspaceError / (1 + subspaceError);

  symmMatrix randomMatrix;
  randomMatrix.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  get_non_ort_err_mat_normalized_in_ort_basis(randomMatrix, 0);

  FockMatrix.readFromFile();

  FockMatrix += desiredErrorNorm * randomMatrix;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Disturbed Fock matrix with desiredErrorNorm = %22.11f  (gap = %8.3f)", 
	    (double)desiredErrorNorm, (double)gap);

  FockMatrix.writeToFile();
}


static ergo_real get_nucl_energy_for_given_mol_and_dens(const IntegralInfo& integralInfo,
							const Molecule& molecule,
							const BasisInfoStruct& basisInfo,
							const symmMatrix & D,
							ergo_real threshold_integrals_1el,
							mat::SizesAndBlocks const & matrix_size_block_info,
							std::vector<int> const & permutationHML) {
  ergo_real nuclearRepulsionEnergy = molecule.getNuclearRepulsionEnergy();
  ergo_real elecNuclEnergy = get_electron_nuclear_attraction_energy(integralInfo,
								    molecule,
								    basisInfo,
								    D,
								    threshold_integrals_1el,
								    matrix_size_block_info,
								    permutationHML);
  return nuclearRepulsionEnergy + elecNuclEnergy;
}


/* Compute gradient of energy with respect to nuclear positions, for
   fixed electron density. */
void
SCF_restricted::compute_gradient_fixeddens() {
  // Since we here regard the electron density as fixed, there are
  // only two terms in the energy which give nonzero contributions:
  // the nuclear-electron interaction term, and the nuclear-nuclear
  // interaction term.
  densityMatrix.readFromFile();
  int nAtoms = molecule.getNoOfAtoms();
  std::vector<ergo_real> gradient(nAtoms*3);
  get_gradient_for_given_mol_and_dens(integralInfo,
				      molecule,
				      basisInfo,
				      densityMatrix,
				      threshold_integrals_1el,
				      matOpts.size_block_info,
				      matOpts.permutationHML,
				      &gradient[0]);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Gradient of energy with respect to nuclear positions, for fixed electron density:");
  for(int i = 0; i < nAtoms; i++) {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Atom %6d: %22.11f %22.11f %22.11f", 
	      i,
	      (double)gradient[i*3+0],
	      (double)gradient[i*3+1],
	      (double)gradient[i*3+2]);
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "(End of gradient)");

  if(scfopts.verify_gradient_fixeddens == 1) {
    std::vector<ergo_real> gradient_for_verification(nAtoms*3);
    for(int i = 0; i < nAtoms; i++) {
      for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
	const ergo_real h = 1e-3;
	Molecule moleculeTmp = molecule;
	Atom atomTmp = molecule.getAtom(i);
	atomTmp.coords[coordIdx] += h;
	moleculeTmp.replaceAtom(i, atomTmp);
	ergo_real E1 = get_nucl_energy_for_given_mol_and_dens(integralInfo,
							      moleculeTmp,
							      basisInfo,
							      densityMatrix,
							      threshold_integrals_1el,
							      matOpts.size_block_info,
							      matOpts.permutationHML);
	moleculeTmp = molecule;
	atomTmp = molecule.getAtom(i);
	atomTmp.coords[coordIdx] -= h;
	moleculeTmp.replaceAtom(i, atomTmp);
	ergo_real E2 = get_nucl_energy_for_given_mol_and_dens(integralInfo,
							      moleculeTmp,
							      basisInfo,
							      densityMatrix,
							      threshold_integrals_1el,
							      matOpts.size_block_info,
							      matOpts.permutationHML);
	ergo_real gradientComponent = (E1 - E2) / (2 * h);
	gradient_for_verification[i*3+coordIdx] = gradientComponent;
      } // END FOR coordIdx
    } // END FOR i
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Gradient of energy with respect to nuclear positions, for fixed electron density (for verification, computed using finite differences):");
    for(int i = 0; i < nAtoms; i++) {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Atom %6d: %22.11f %22.11f %22.11f", 
		i,
		(double)gradient_for_verification[i*3+0],
		(double)gradient_for_verification[i*3+1],
		(double)gradient_for_verification[i*3+2]);
    }
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "(End of gradient)");
    ergo_real maxAbsDiff = -1;
    int atomIdx_saved = -1;
    int coordIdx_saved = -1;
    for(int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
      for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
	ergo_real absdiff = std::fabs(gradient[atomIdx*3+coordIdx] - gradient_for_verification[atomIdx*3+coordIdx]);
	if(absdiff > maxAbsDiff) {
	  maxAbsDiff = absdiff;
	  atomIdx_saved = atomIdx;
	  coordIdx_saved = coordIdx;
	}
      }
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxAbsDiff %22.11f = %9.4g, for atomIdx %d and coordIdx %d", maxAbsDiff, maxAbsDiff, atomIdx_saved, coordIdx_saved);
  }

  densityMatrix.writeToFile();
}



