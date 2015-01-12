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
#include "SCF_unrestricted.h"
#include "output.h"
#include "scf_utils.h"
#include "utilities.h"
#include "diis_unrestricted.h"
#include "density_projection.h"
#include "densfromf_sparse.h"
#include "densfromf_full.h"
#include "densfromf_general.h"
#include "density_description_file.h"
#include "densitymanager.h"
#include "matrix_utilities.h"
#include "units.h"


SCF_unrestricted::SCF_unrestricted(const Molecule& molecule_,
				   const Molecule& extraCharges_,
				   const BasisInfoStruct & basisInfo_, 
				   const BasisInfoStruct & basisInfoDensFit_,
				   const IntegralInfo& integralInfo_,
				   const char* guessDmatFileName_,
				   const JK::Params& J_K_params_,
				   const Dft::GridParams& gridParams_,
				   const SCF::Options& scfopts,
				   const SCF::MatOptions& matOpts,
				   ergo_real threshold_integrals_1el_input,
				   int alpha_beta_diff_input)
  :   SCF_general(molecule_,
		  extraCharges_,
		  basisInfo_, 
		  basisInfoDensFit_,
		  integralInfo_,
		  guessDmatFileName_,
		  J_K_params_,
		  gridParams_,
		  scfopts,
		  matOpts,
		  threshold_integrals_1el_input),
      alpha_beta_diff(alpha_beta_diff_input)
{
  DIIS = new DIISManagerUnrestricted;
  if(determine_number_of_electrons_unrestricted(noOfElectrons, alpha_beta_diff, &noOfElectrons_alpha, &noOfElectrons_beta) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in determine_number_of_electrons_unrestricted");
      throw "error in determine_number_of_electrons_unrestricted";
    }
}

SCF_unrestricted::~SCF_unrestricted()
{
  delete ((DIISManagerUnrestricted*)DIIS);
}


void SCF_unrestricted::get_Fock_matrices(symmMatrix & FockMatrix_a, symmMatrix & FockMatrix_b)
{
  FockMatrix_alpha.readFromFile();
  FockMatrix_a = FockMatrix_alpha;
  FockMatrix_alpha.writeToFile();
  FockMatrix_beta.readFromFile();
  FockMatrix_b = FockMatrix_beta;
  FockMatrix_beta.writeToFile();
}


void SCF_unrestricted::get_no_of_electrons(int & noOfElectrons_a, int & noOfElectrons_b)
{
  noOfElectrons_a = noOfElectrons_alpha;
  noOfElectrons_b = noOfElectrons_beta;
}


void SCF_unrestricted::initialize_matrices()
{
  densityMatrix_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  densityMatrix_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  FockMatrix_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  FockMatrix_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  Fprev_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  Fprev_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  Dprev_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  Dprev_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  F_ort_prev_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  F_ort_prev_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  bestFockMatrixSoFar_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  bestFockMatrixSoFar_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  bestFockMatrixSoFar2_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  bestFockMatrixSoFar2_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  ErrorMatrix_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  ErrorMatrix_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  G_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  G_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
}



void SCF_unrestricted::check_params()
{
}




void SCF_unrestricted::get_starting_guess_density()
{
  // set up starting guess

  int n = basisInfo.noOfBasisFuncs;

  if(guessDmatFileName != NULL)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "getting starting guess density from file '%s'", guessDmatFileName);
      int noOfDensityMatrices = 2;


      symmMatrix* matrixList[2];
      matrixList[0] = &densityMatrix_alpha;
      matrixList[1] = &densityMatrix_beta;

      int noOfElectronsList[2];
      noOfElectronsList[0] = noOfElectrons_alpha + scfopts.starting_guess_spin_diff;
      noOfElectronsList[1] = noOfElectrons_beta  - scfopts.starting_guess_spin_diff;

      if(load_density_and_project_sparse(guessDmatFileName,
					 noOfDensityMatrices,
					 &integralInfo,
					 basisInfo,
					 S_symm,
					 matrixList,
					 noOfElectronsList,
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
	  if(get_simple_starting_guess_sparse(n, noOfElectrons_alpha, densityMatrix_alpha) != 0)
	    throw "error in get_simple_starting_guess_sparse";
	  densityMatrix_alpha.writeToFile();
	  if(get_simple_starting_guess_sparse(n, noOfElectrons_beta, densityMatrix_beta) != 0)
	    throw "error in get_simple_starting_guess_sparse";
	  densityMatrix_beta.writeToFile();
	}
      else if(scfopts.use_diag_guess_from_file == 1)
	{
	  if(get_diag_matrix_from_file(n, densityMatrix_alpha, "diagdens_alpha.txt",matOpts.permutationHML) != 0)
	    throw "error in get_diag_matrix_from_file";
	  if(get_diag_matrix_from_file(n, densityMatrix_beta , "diagdens_beta.txt",matOpts.permutationHML) != 0)
	    throw "error in get_diag_matrix_from_file";
	}
      else
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		    "calling get_dens_from_fock_general for alpha and beta to diagonalize H_core for starting guesses, n = %i, matOpts.sparse_threshold = %g", 
		    n, (double)matOpts.sparse_threshold);

	  // ALPHA
	  symmMatrix F_ort_prev_dummy_alpha;
	  F_ort_prev_dummy_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
	  F_ort_prev_dummy_alpha.writeToFile();
	  intervalType homoInterval_dummy1_alpha(-1e22,1e22);
	  intervalType lumoInterval_dummy1_alpha(-1e22,1e22);
	  intervalType homoInterval_dummy2_alpha(-1e22,1e22);
	  intervalType lumoInterval_dummy2_alpha(-1e22,1e22);
	  densityMatrix_alpha.writeToFile();
	  std::map<std::string, double> puri_stats_alpha;
	  ergo_real electronicEntropyTermDummy = 0;
	  std::vector< std::vector<ergo_real> > stochasticOrbsDummy;
	  if(get_dens_from_fock_general(n,
					noOfElectrons_alpha,
					scfopts.use_diagonalization,
					scfopts.use_diag_on_error_guess,
					scfopts.electronic_temperature,
					densityMatrix_alpha,
					1,
					electronicEntropyTermDummy,
					H_core_Matrix, 
					homoInterval_dummy1_alpha,
					lumoInterval_dummy1_alpha,
					S_symm,
					invCholFactor,
					invCholFactor_euclnorm,
					scfopts.gap_expected_lower_bound,
					matOpts.size_block_info,
					F_ort_prev_dummy_alpha,
					homoInterval_dummy2_alpha,
					lumoInterval_dummy2_alpha,
					scfopts.purification_eigvalue_err_limit,
					scfopts.purification_subspace_err_limit,
					scfopts.purification_truncation_norm,
					scfopts.purification_maxmul,
					scfopts.purification_create_m_files,
					scfopts.purification_ignore_failure,
					scfopts.purification_use_rand_perturbation_for_alleigsint,
					0,
					stochasticOrbsDummy,
					"",
					puri_stats_alpha,
					scfopts.do_sparsity_investigation,
					scfopts.sparsity_plots_resolution_m,
					scfopts.do_comparison_to_simple_purification,
					scfopts.do_puri_mmul_tests) != 0)
	    {
	      throw "SCF_unrestricted::get_starting_guess_density: Error in get_dens_from_fock_general for alpha.";
	    }

	  // BETA
	  symmMatrix F_ort_prev_dummy_beta;
	  F_ort_prev_dummy_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
	  F_ort_prev_dummy_beta.writeToFile();
	  intervalType homoInterval_dummy1_beta(-1e22,1e22);
	  intervalType lumoInterval_dummy1_beta(-1e22,1e22);
	  intervalType homoInterval_dummy2_beta(-1e22,1e22);
	  intervalType lumoInterval_dummy2_beta(-1e22,1e22);
	  densityMatrix_beta.writeToFile();
	  std::map<std::string, double> puri_stats_beta;
	  if(get_dens_from_fock_general(n,
					noOfElectrons_beta,
					scfopts.use_diagonalization,
					scfopts.use_diag_on_error_guess,
					scfopts.electronic_temperature,
					densityMatrix_beta,
					1,
					electronicEntropyTermDummy,
					H_core_Matrix, 
					homoInterval_dummy1_beta,
					lumoInterval_dummy1_beta,
					S_symm,
					invCholFactor,
					invCholFactor_euclnorm,
					scfopts.gap_expected_lower_bound,
					matOpts.size_block_info,
					F_ort_prev_dummy_beta,
					homoInterval_dummy2_beta,
					lumoInterval_dummy2_beta,
					scfopts.purification_eigvalue_err_limit,
					scfopts.purification_subspace_err_limit,
					scfopts.purification_truncation_norm,
					scfopts.purification_maxmul,
					scfopts.purification_create_m_files,
					scfopts.purification_ignore_failure,
					scfopts.purification_use_rand_perturbation_for_alleigsint,
					0,
					stochasticOrbsDummy,
					"",
					puri_stats_beta,
					scfopts.do_sparsity_investigation,
					scfopts.sparsity_plots_resolution_m,
					scfopts.do_comparison_to_simple_purification,
					scfopts.do_puri_mmul_tests) != 0)
	    {
	      throw "SCF_unrestricted::get_starting_guess_density: Error in get_dens_from_fock_general for beta.";
	    }

	} // END ELSE use H_core to get starting guesses
    } // END ELSE no dmat given

  densityMatrix_alpha.readFromFile();
  densityMatrix_beta.readFromFile();
  output_sparsity_symm(n, densityMatrix_alpha, "starting guess density matrix (alpha)");
  output_sparsity_symm(n, densityMatrix_beta , "starting guess density matrix (beta )");
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.writeToFile();
  
}



void SCF_unrestricted::add_random_disturbance_to_starting_guess()
{
  if(scfopts.sg_disturb_specific_elements > SCF::DISTURB_ELEMENT_MAX_COUNT)
    throw "Error in SCF_unrestricted::add_random_disturbance_to_starting_guess: (scfopts.sg_disturb_specific_elements > SCF::DISTURB_ELEMENT_MAX_COUNT)";
  int n = basisInfo.noOfBasisFuncs;
  densityMatrix_alpha.readFromFile();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "SCF_unrestricted::add_random_disturbance_to_starting_guess, scfopts.starting_guess_disturbance = %7.3f",
	    scfopts.starting_guess_disturbance);
  add_disturbance_to_matrix(n, 
			    densityMatrix_alpha, 
			    scfopts.starting_guess_disturbance,
			    scfopts.sg_disturb_specific_elements,
			    scfopts.disturbedElementIndexVector,
			    matOpts.permutationHML);
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.readFromFile();
  add_disturbance_to_matrix(n, 
			    densityMatrix_beta , 
			    -1 * scfopts.starting_guess_disturbance,
			    scfopts.sg_disturb_specific_elements,
			    scfopts.disturbedElementIndexVector,
			    matOpts.permutationHML);
  densityMatrix_beta.writeToFile();
}



void SCF_unrestricted::initialize_homo_lumo_limits()
{
  intervalType hugeInterval(-1e22, 1e22);
  homoInterval_F_ort_prev_alpha = hugeInterval;
  lumoInterval_F_ort_prev_alpha = hugeInterval;
  homoInterval_F_ort_prev_beta = hugeInterval;
  lumoInterval_F_ort_prev_beta = hugeInterval;
  homoInterval_Fprev_alpha = hugeInterval;
  lumoInterval_Fprev_alpha = hugeInterval;
  homoInterval_Fprev_beta = hugeInterval;
  lumoInterval_Fprev_beta = hugeInterval;
}



void SCF_unrestricted::write_matrices_to_file()
{
  FockMatrix_alpha.writeToFile();
  FockMatrix_beta.writeToFile();

  Fprev_alpha.writeToFile();
  Fprev_beta.writeToFile();

  Dprev_alpha.writeToFile();
  Dprev_beta.writeToFile();

  bestFockMatrixSoFar_alpha.writeToFile();
  bestFockMatrixSoFar_beta.writeToFile();

  bestFockMatrixSoFar2_alpha.writeToFile();
  bestFockMatrixSoFar2_beta.writeToFile();

  F_ort_prev_alpha.writeToFile();
  F_ort_prev_beta.writeToFile();
}



template<typename T>
static void
printMat(const T& m, const char *msg, int nbast)
{
#if 0
  ergo_real* m_full = new ergo_real[nbast*nbast];
  m.fullmatrix(m_full, nbast, nbast);
  puts(msg);
  for(int row=0; row<nbast; row++) {
    for(int col=0; col<nbast; col++)
      printf(" %10.6f", m_full[row + col*nbast]);
    puts("");
  }
  delete []m_full;
#endif
}

void SCF_unrestricted::get_2e_part_and_energy()
{
  // calculate Fock matrices 
  if(!scfopts.force_restricted) {
    // F_alpha = H_core + G_alpha
    // F_beta  = H_core + G_beta
    densityMatrix_alpha.readFromFile();
    densityMatrix_beta.readFromFile();
    if(get_2e_matrices_and_energy_sparse_unrestricted(basisInfo,
						      basisInfoDensFit, 
						      molecule,
						      integralInfo, 
						      CAM_params,
						      G_alpha,
						      G_beta,
						      densityMatrix_alpha, 
						      densityMatrix_beta, 
						      J_K_params,
						      gridParams,
						      scfopts.use_dft,
						      &energy_2el,
						      noOfElectrons,
						      densfit_data,
						      matOpts.size_block_info,
						      matOpts.permutationHML,
						      matOpts.inversePermutationHML) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrices_and_energy_sparse_unrestricted");
	throw "error in get_2e_matrices_and_energy_sparse_unrestricted";
      }
    densityMatrix_alpha.writeToFile();
    densityMatrix_beta.writeToFile();

    H_core_Matrix.readFromFile();
    FockMatrix_alpha.readFromFile();
    FockMatrix_alpha = H_core_Matrix + G_alpha;
    FockMatrix_alpha.frob_thresh(matOpts.sparse_threshold);
    printMat(FockMatrix_alpha, "FockMatrix_alpha", basisInfo.noOfBasisFuncs);
    FockMatrix_alpha.writeToFile();
    FockMatrix_beta.readFromFile();
    FockMatrix_beta = H_core_Matrix + G_beta;
    FockMatrix_beta.frob_thresh(matOpts.sparse_threshold);
    FockMatrix_beta.writeToFile();
    H_core_Matrix.writeToFile();
  } else {
    // F_alpha = F_beta = Frestricted
    densityMatrix_alpha.readFromFile();
    densityMatrix_beta.readFromFile();
    if(get_2e_matrices_and_energy_restricted_open(basisInfo,
						  basisInfoDensFit, 
						  molecule,
						  integralInfo, 
						  CAM_params,
						  G_alpha, /*J_a+J_b+X_a+X_b*/
						  G_beta,  /*J_a+J_b+X_a*/
						  densityMatrix_alpha, 
						  densityMatrix_beta, 
						  J_K_params,
						  gridParams,
						  scfopts.use_dft,
						  &energy_2el,
						  noOfElectrons,
						  densfit_data,
						  matOpts.size_block_info,
						  matOpts.permutationHML,
						  matOpts.inversePermutationHML) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrices_and_energy_restricted_open");
	throw "error in get_2e_matrices_and_energy_restricted_open";
      }

    printMat(densityMatrix_alpha, "densityMatrix_alpha",
	     basisInfo.noOfBasisFuncs);
    printMat(densityMatrix_beta, "densityMatrix_beta",
	     basisInfo.noOfBasisFuncs);

    symmMatrix S, Fc, Fo, FcMinusFo, fockConstrained, syTmp;
    normalMatrix Sc, So, tmp, tmp1, tmp2;
    Sc.resetSizesAndBlocks(matOpts.size_block_info,
				  matOpts.size_block_info);
    So.resetSizesAndBlocks(matOpts.size_block_info,
				  matOpts.size_block_info);
    tmp.resetSizesAndBlocks(matOpts.size_block_info,
				   matOpts.size_block_info);
    tmp1.resetSizesAndBlocks(matOpts.size_block_info,
				    matOpts.size_block_info);
    tmp2.resetSizesAndBlocks(matOpts.size_block_info,
				    matOpts.size_block_info);
    fockConstrained.resetSizesAndBlocks(matOpts.size_block_info,
					       matOpts.size_block_info);

    H_core_Matrix.readFromFile();
    FockMatrix_alpha.readFromFile();
    FockMatrix_beta.readFromFile();

    get_overlap_matrix(S);
    printMat(S, "overlap",  basisInfo.noOfBasisFuncs);

    syTmp = densityMatrix_alpha - densityMatrix_beta;
    So = syTmp*S;

    Sc = densityMatrix_beta*S;
    printMat(Sc, "Sc", basisInfo.noOfBasisFuncs);

    Fc = H_core_Matrix + G_alpha;
    Fo = H_core_Matrix + G_beta;  Fo *= 0.5;

    printMat(Fc, "Fc", basisInfo.noOfBasisFuncs);
    printMat(Fo, "Fo", basisInfo.noOfBasisFuncs);
    tmp = 1;
    tmp -= So;
    //    tmp = (-1.0)*So + one;
    tmp1 = Fc*tmp;
    tmp2 = transpose(tmp1)*tmp;
    fockConstrained = tmp2;

    tmp = 1;
    tmp -= Sc;
    tmp1 = Fo * tmp;
    tmp2 = transpose(tmp1) * tmp;
    syTmp = tmp2;
    printMat(tmp2, "a2", basisInfo.noOfBasisFuncs);
    fockConstrained += syTmp;

    printMat(fockConstrained, "a1+a2", basisInfo.noOfBasisFuncs);

    FcMinusFo = Fc - Fo;
    
    tmp2 = FcMinusFo*Sc;
    tmp1 = transpose(tmp2)*So;
    tmp = transpose(tmp1);
    tmp1 += tmp;
    syTmp = tmp1;
    fockConstrained += syTmp;
    printMat(syTmp, "a3", basisInfo.noOfBasisFuncs);

    printMat(fockConstrained, "a1+a2+a3", basisInfo.noOfBasisFuncs);

    FockMatrix_alpha = fockConstrained;
    FockMatrix_beta  = fockConstrained;
    FockMatrix_alpha.writeToFile();
    FockMatrix_beta.writeToFile();
    H_core_Matrix.writeToFile();
    densityMatrix_alpha.writeToFile();
    densityMatrix_beta.writeToFile();
  }
}




void SCF_unrestricted::output_sparsity_S_F_D(SCF_statistics & stats)
{
  int n = basisInfo.noOfBasisFuncs;
  S_symm.readFromFile();
  output_sparsity_symm(n, S_symm, "S");
  S_symm.writeToFile();
  FockMatrix_alpha.readFromFile();
  output_sparsity_symm(n, FockMatrix_alpha, "F_alpha");
  FockMatrix_alpha.writeToFile();
  FockMatrix_beta.readFromFile();
  output_sparsity_symm(n, FockMatrix_beta, "F_beta ");
  FockMatrix_beta.writeToFile();
  densityMatrix_alpha.readFromFile();
  output_sparsity_symm(n, densityMatrix_alpha, "D_alpha");
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.readFromFile();
  output_sparsity_symm(n, densityMatrix_beta, "D_beta ");
  densityMatrix_beta.writeToFile();
}




void SCF_unrestricted::calculate_energy()
{
  // calculate energy
  H_core_Matrix.readFromFile();
  densityMatrix_alpha.readFromFile();
  densityMatrix_beta.readFromFile();
  energy = 0;
  energy += symmMatrix::trace_ab(densityMatrix_alpha, H_core_Matrix);
  energy += symmMatrix::trace_ab(densityMatrix_beta , H_core_Matrix);
  energy += energy_2el;
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.writeToFile();
  H_core_Matrix.writeToFile();
  energy += nuclearEnergy;
}




void SCF_unrestricted::get_FDSminusSDF()
{
  int n = basisInfo.noOfBasisFuncs;

  // ALPHA
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling compute_FDSminusSDF_sparse for ALPHA, n = %i", n);
  densityMatrix_alpha.readFromFile();
  S_symm.readFromFile();
  FockMatrix_alpha.readFromFile();
  compute_FDSminusSDF_sparse(n, FockMatrix_alpha, densityMatrix_alpha, S_symm,
			     ErrorMatrix_alpha, matOpts.sparse_threshold);
  S_symm.writeToFile();
  FockMatrix_alpha.writeToFile();
  densityMatrix_alpha.writeToFile();
  // write to file and read back again to reduce memory fragmentation.
  ErrorMatrix_alpha.writeToFile();
  ErrorMatrix_alpha.readFromFile();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_FDSminusSDF_sparse for ALPHA finished.");

  // BETA
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling compute_FDSminusSDF_sparse for BETA, n = %i", n);
  densityMatrix_beta.readFromFile();
  S_symm.readFromFile();
  FockMatrix_beta.readFromFile();
  compute_FDSminusSDF_sparse(n, FockMatrix_beta, densityMatrix_beta, S_symm,
			     ErrorMatrix_beta, matOpts.sparse_threshold);
  S_symm.writeToFile();
  FockMatrix_beta.writeToFile();
  densityMatrix_beta.writeToFile();
  // write to file and read back again to reduce memory fragmentation.
  ErrorMatrix_beta.writeToFile();
  ErrorMatrix_beta.readFromFile();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_FDSminusSDF_sparse for BETA finished.");

  output_sparsity(n, ErrorMatrix_alpha, "FDS-SDF (alpha)");
  output_sparsity(n, ErrorMatrix_beta , "FDS-SDF (beta )");
}




void SCF_unrestricted::get_error_measure()
{
  ergo_real error_maxabs_alpha = compute_maxabs_sparse(ErrorMatrix_alpha);
  ergo_real error_maxabs_beta  = compute_maxabs_sparse(ErrorMatrix_beta );
  
  ergo_real error_frob_alpha = ErrorMatrix_alpha.frob();
  ergo_real error_frob_beta  = ErrorMatrix_beta .frob();
  
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxabs FDS-SDF (alpha) is %8.3g", (double)error_maxabs_alpha);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxabs FDS-SDF (beta ) is %8.3g", (double)error_maxabs_beta);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "frob   FDS-SDF (alpha) is %8.3g", (double)error_frob_alpha);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "frob   FDS-SDF (beta ) is %8.3g", (double)error_frob_beta);

  errorMeasure = 0.5 * (error_maxabs_alpha + error_maxabs_beta);
}





void SCF_unrestricted::add_to_DIIS_list()
{
  FockMatrix_alpha.readFromFile();
  FockMatrix_beta.readFromFile();
  if(((DIISManagerUnrestricted*)DIIS)->AddIterationToList(FockMatrix_alpha, 
							  FockMatrix_beta,
							  ErrorMatrix_alpha,
							  ErrorMatrix_beta) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in DIIS AddIterationToList");
      throw "error in DIIS AddIterationToList";
    }
  FockMatrix_alpha.writeToFile();
  FockMatrix_beta.writeToFile();
}





void SCF_unrestricted::update_best_fock_so_far()
{
  // ALPHA
  Fprev_alpha.readFromFile();
  bestFockMatrixSoFar_alpha.readFromFile();
  bestFockMatrixSoFar_alpha = Fprev_alpha;
  Fprev_alpha.writeToFile();
  bestFockMatrixSoFar_alpha.writeToFile();
  FockMatrix_alpha.readFromFile();
  bestFockMatrixSoFar2_alpha.readFromFile();
  bestFockMatrixSoFar2_alpha = FockMatrix_alpha;
  FockMatrix_alpha.writeToFile();
  bestFockMatrixSoFar2_alpha.writeToFile();
  // BETA
  Fprev_beta.readFromFile();
  bestFockMatrixSoFar_beta.readFromFile();
  bestFockMatrixSoFar_beta = Fprev_beta;
  Fprev_beta.writeToFile();
  bestFockMatrixSoFar_beta.writeToFile();
  FockMatrix_beta.readFromFile();
  bestFockMatrixSoFar2_beta.readFromFile();
  bestFockMatrixSoFar2_beta = FockMatrix_beta;
  FockMatrix_beta.writeToFile();
  bestFockMatrixSoFar2_beta.writeToFile();
}




void SCF_unrestricted::combine_old_fock_matrices(ergo_real stepLength)
{
  // ALPHA
  bestFockMatrixSoFar_alpha.readFromFile();
  bestFockMatrixSoFar2_alpha.readFromFile();
  FockMatrix_alpha.readFromFile();
  FockMatrix_alpha = 0;
  FockMatrix_alpha += stepLength * bestFockMatrixSoFar2_alpha;
  FockMatrix_alpha += (1 - stepLength) * bestFockMatrixSoFar_alpha;
  FockMatrix_alpha.writeToFile();
  bestFockMatrixSoFar_alpha.writeToFile();
  bestFockMatrixSoFar2_alpha.writeToFile();
  // BETA
  bestFockMatrixSoFar_beta.readFromFile();
  bestFockMatrixSoFar2_beta.readFromFile();
  FockMatrix_beta.readFromFile();
  FockMatrix_beta = 0;
  FockMatrix_beta += stepLength * bestFockMatrixSoFar2_beta;
  FockMatrix_beta += (1 - stepLength) * bestFockMatrixSoFar_beta;
  FockMatrix_beta.writeToFile();
  bestFockMatrixSoFar_beta.writeToFile();
  bestFockMatrixSoFar2_beta.writeToFile();
}




void SCF_unrestricted::use_diis_to_get_new_fock_matrix()
{
  symmMatrix newFsymm_alpha;
  symmMatrix newFsymm_beta;
  if(((DIISManagerUnrestricted*)DIIS)->GetCombinedFockMatrices(newFsymm_alpha, newFsymm_beta) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in DIIS GetCombinedFockMatrices");
      throw "error in DIIS GetCombinedFockMatrices";
    }
  // ALPHA
  FockMatrix_alpha.readFromFile();
  FockMatrix_alpha = newFsymm_alpha;
  FockMatrix_alpha.frob_thresh(matOpts.sparse_threshold);
  FockMatrix_alpha.writeToFile();
  // BETA
  FockMatrix_beta.readFromFile();
  FockMatrix_beta = newFsymm_beta;
  FockMatrix_beta.frob_thresh(matOpts.sparse_threshold);
  FockMatrix_beta.writeToFile();
}





void SCF_unrestricted::clear_diis_list()
{
  ((DIISManagerUnrestricted*)DIIS)->ClearList();
}




void SCF_unrestricted::clear_error_matrices()
{
  ErrorMatrix_alpha.clear();
  ErrorMatrix_beta.clear();
}




void SCF_unrestricted::save_current_fock_as_fprev()
{
  // ALPHA
  FockMatrix_alpha.readFromFile();
  Fprev_alpha.readFromFile();
  Fprev_alpha = FockMatrix_alpha;
  FockMatrix_alpha.writeToFile();
  Fprev_alpha.writeToFile();
  // BETA
  FockMatrix_beta.readFromFile();
  Fprev_beta.readFromFile();
  Fprev_beta = FockMatrix_beta;
  FockMatrix_beta.writeToFile();
  Fprev_beta.writeToFile();
}




void SCF_unrestricted::get_new_density_matrix()
{
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "SCF_unrestricted::get_new_density_matrix, noOfElectrons_alpha noOfElectrons_beta = %5i %5i",
	    noOfElectrons_alpha, noOfElectrons_beta);

  // ALPHA
  std::map<std::string, double> puri_stats_alpha;
  ergo_real electronicEntropyTerm_alpha = 0;
  std::vector< std::vector<ergo_real> > stochasticOrbsDummy;
  if(get_dens_from_fock_general(n,
				noOfElectrons_alpha, 
				scfopts.use_diagonalization,
				scfopts.use_diag_on_error,
				scfopts.electronic_temperature,
				densityMatrix_alpha,
				1,
				electronicEntropyTerm_alpha,
				FockMatrix_alpha, 
				homoInterval_Fprev_alpha,
				lumoInterval_Fprev_alpha,
				S_symm,
				invCholFactor,
				invCholFactor_euclnorm,
				scfopts.gap_expected_lower_bound,
				matOpts.size_block_info,
				F_ort_prev_alpha,
				homoInterval_F_ort_prev_alpha,
				lumoInterval_F_ort_prev_alpha,
				scfopts.purification_eigvalue_err_limit,
				scfopts.purification_subspace_err_limit,
				scfopts.purification_truncation_norm,
				scfopts.purification_maxmul,
				scfopts.purification_create_m_files,
				scfopts.purification_ignore_failure,
				scfopts.purification_use_rand_perturbation_for_alleigsint,
				0,
				stochasticOrbsDummy,
				"",
				puri_stats_alpha,
				scfopts.do_sparsity_investigation,
				scfopts.sparsity_plots_resolution_m,
				scfopts.do_comparison_to_simple_purification,
				scfopts.do_puri_mmul_tests) != 0)
    {
      throw "SCF_unrestricted::get_new_density_matrix: Error in get_dens_from_fock_general for alpha.";
    }

  // BETA
  std::map<std::string, double> puri_stats_beta;
  ergo_real electronicEntropyTerm_beta = 0;
  if(get_dens_from_fock_general(n,
				noOfElectrons_beta,
				scfopts.use_diagonalization,
				scfopts.use_diag_on_error,
				scfopts.electronic_temperature,
				densityMatrix_beta,
				1,
				electronicEntropyTerm_beta,
				FockMatrix_beta, 
				homoInterval_Fprev_beta,
				lumoInterval_Fprev_beta,
				S_symm,
				invCholFactor,
				invCholFactor_euclnorm,
				scfopts.gap_expected_lower_bound,
				matOpts.size_block_info,
				F_ort_prev_beta,
				homoInterval_F_ort_prev_beta,
				lumoInterval_F_ort_prev_beta,
				scfopts.purification_eigvalue_err_limit,
				scfopts.purification_subspace_err_limit,
				scfopts.purification_truncation_norm,
				scfopts.purification_maxmul,
				scfopts.purification_create_m_files,
				scfopts.purification_ignore_failure,
				scfopts.purification_use_rand_perturbation_for_alleigsint,
				0,
				stochasticOrbsDummy,
				"",
				puri_stats_beta,
				scfopts.do_sparsity_investigation,
				scfopts.sparsity_plots_resolution_m,
				scfopts.do_comparison_to_simple_purification,
				scfopts.do_puri_mmul_tests) != 0)
    {
      throw "SCF_unrestricted::get_new_density_matrix: Error in get_dens_from_fock_general for beta.";
    }

  electronicEntropyTerm = electronicEntropyTerm_alpha + electronicEntropyTerm_beta;

  // Report trace(DS) and spin info.
  S_symm.readFromFile();
  densityMatrix_alpha.readFromFile();
  densityMatrix_beta.readFromFile();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Tr( D_alpha * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix_alpha, S_symm));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Tr( D_beta  * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix_beta , S_symm));
  symmMatrix spinDensityMatrix(densityMatrix_alpha);
  spinDensityMatrix += (ergo_real)(-1.0) * densityMatrix_beta;
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.writeToFile();
  S_symm.writeToFile();
  normalMatrix spinDensityMatrix_normal(spinDensityMatrix);
  spinDensityMatrix.clear();
  ergo_real maxabs = compute_maxabs_sparse(spinDensityMatrix_normal);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "maxabs element in spin density matrix : %22.11f", (double)maxabs);

  ergo_real S2_exact, S2;
  get_S2(S2_exact, S2);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "<S2>_exact = %7.3f, <S2> = %22.11f", (double)S2_exact, (double)S2);
}



void SCF_unrestricted::write_density_to_file()
{
  matrix_description_struct matrixList[2];
  // ALPHA
  densityMatrix_alpha.readFromFile();
  int nvalues_alpha = densityMatrix_alpha.nvalues();
  std::vector<int> rowind_alpha;
  rowind_alpha.reserve(nvalues_alpha);
  std::vector<int> colind_alpha;
  colind_alpha.reserve(nvalues_alpha);
  std::vector<ergo_real> values_alpha;
  values_alpha.reserve(nvalues_alpha);
  densityMatrix_alpha.get_all_values(rowind_alpha,
				     colind_alpha,
				     values_alpha,
				     matOpts.inversePermutationHML,
				     matOpts.inversePermutationHML);
  densityMatrix_alpha.writeToFile();
  matrixList[0].nvalues = nvalues_alpha;
  matrixList[0].rowind = &rowind_alpha[0];
  matrixList[0].colind = &colind_alpha[0];
  matrixList[0].values = &values_alpha[0];
  // BETA
  densityMatrix_beta.readFromFile();
  int nvalues_beta = densityMatrix_beta.nvalues();
  std::vector<int> rowind_beta;
  rowind_beta.reserve(nvalues_beta);
  std::vector<int> colind_beta;
  colind_beta.reserve(nvalues_beta);
  std::vector<ergo_real> values_beta;
  values_beta.reserve(nvalues_beta);
  densityMatrix_beta.get_all_values(rowind_beta,
				    colind_beta,
				    values_beta,
				    matOpts.inversePermutationHML,
				    matOpts.inversePermutationHML);
  densityMatrix_beta.writeToFile();
  matrixList[1].nvalues = nvalues_beta;
  matrixList[1].rowind = &rowind_beta[0];
  matrixList[1].colind = &colind_beta[0];
  matrixList[1].values = &values_beta[0];

  if(ddf_writeShellListAndDensityMatricesToFile_sparse(basisInfo,
						       2,
						       matrixList,
						       "density.bin") != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ddf_writeShellListAndDensityMatricesToFile_sparse");
      throw "error in ddf_writeShellListAndDensityMatricesToFile_sparse";
    }

}


void SCF_unrestricted::save_final_potential()
{
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, 
		"error: save_final_potential not implemented for unrestricted case.");
      throw "error: save_final_potential not implemented for unrestricted case.";
}


void SCF_unrestricted::output_density_images()
{
  Util::TimeMeter timeMeter;

  int n = basisInfo.noOfBasisFuncs;

  ergo_real* densityMatrixFull_tot  = new ergo_real[n*n];
  ergo_real* densityMatrixFull_spin = new ergo_real[n*n];

  // Get full matrix versions of density matrices
  {
    std::vector<ergo_real> densityMatrixFull_a(n*n);
    std::vector<ergo_real> densityMatrixFull_b(n*n);
    
    densityMatrix_alpha.readFromFile();
    densityMatrix_alpha.fullMatrix(densityMatrixFull_a,
				   matOpts.inversePermutationHML,
				   matOpts.inversePermutationHML);
    densityMatrix_alpha.writeToFile();
    
    densityMatrix_beta.readFromFile();
    densityMatrix_beta.fullMatrix(densityMatrixFull_b,
				  matOpts.inversePermutationHML,
				  matOpts.inversePermutationHML);
    densityMatrix_beta.writeToFile();
    
    for(int i = 0; i < n*n; i++)
      {
	densityMatrixFull_tot [i] = densityMatrixFull_a[i] + densityMatrixFull_b[i];
	densityMatrixFull_spin[i] = densityMatrixFull_a[i] - densityMatrixFull_b[i];
      }
  }

  do_density_images(basisInfo,
		    molecule,
		    densityMatrixFull_tot, 
		    densityMatrixFull_spin,
		    scfopts.output_density_images_boxwidth);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_unrestricted::output_density_images finished OK.");
  timeMeter.print(LOG_AREA_SCF, "SCF_unrestricted::output_density_images");
}


void SCF_unrestricted::prepare_stochastic_orbitals() {
  throw "Error: SCF_unrestricted::prepare_stochastic_orbitals() not implemented.";
}


void SCF_unrestricted::do_spin_flip(int atomCount)
{
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_unrestricted::do_spin_flip, n = %6i", n);

  ergo_real* densityMatrixFull_tot  = new ergo_real[n*n];
  ergo_real* densityMatrixFull_spin = new ergo_real[n*n];

  // Get full matrix versions of density matrices
  std::vector<ergo_real> densityMatrixFull_a(n*n);
  std::vector<ergo_real> densityMatrixFull_b(n*n);
  
  densityMatrix_alpha.readFromFile();
  densityMatrix_beta.readFromFile();
  S_symm.readFromFile();

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Before spinflip: Tr( D_alpha * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix_alpha, S_symm));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Before spinflip: Tr( D_beta  * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix_beta , S_symm));

  densityMatrix_alpha.fullMatrix(densityMatrixFull_a,
				 matOpts.inversePermutationHML,
				 matOpts.inversePermutationHML);  
  densityMatrix_beta.fullMatrix(densityMatrixFull_b,
				matOpts.inversePermutationHML,
				matOpts.inversePermutationHML);
  
  for(int i = 0; i < n*n; i++)
    {
      densityMatrixFull_tot [i] = densityMatrixFull_a[i] + densityMatrixFull_b[i];
      densityMatrixFull_spin[i] = densityMatrixFull_a[i] - densityMatrixFull_b[i];
    }

  // Now modify spin density matrix
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	// get atom index for basis funcs i and j
	int atomIndex_i = -1;
	int atomIndex_j = -1;
	for(int k = 0; k < molecule.getNoOfAtoms(); k++)
	  {
	    int iok = 1;
	    int jok = 1;
	    for(int coord = 0; coord < 3; coord++)
	      {
		if(std::fabs(molecule.getAtom(k).coords[coord] - basisInfo.basisFuncList[i].centerCoords[coord]) > 1e-5)
		  iok = 0;
		if(std::fabs(molecule.getAtom(k).coords[coord] - basisInfo.basisFuncList[j].centerCoords[coord]) > 1e-5)
		  jok = 0;
	      }
	    if(iok == 1)
	      atomIndex_i = k;
	    if(jok == 1)
	      atomIndex_j = k;
	  } // END FOR k
	
	// Change spin density matrix element if both basis functions belong to a spin-flip atom.
	if(atomIndex_i < atomCount && atomIndex_j < atomCount)
	  densityMatrixFull_spin[i*n+j] *= -1;
	
      } // END FOR i j

  // Now the spin density matrix has been modified. Recreate new alpha- and beta- density matrices.
  for(int i = 0; i < n*n; i++)
    {
      densityMatrixFull_a[i] = 0.5 * (densityMatrixFull_tot[i] + densityMatrixFull_spin[i]);
      densityMatrixFull_b[i] = 0.5 * (densityMatrixFull_tot[i] - densityMatrixFull_spin[i]);
    }

  densityMatrix_alpha.assignFromFull(densityMatrixFull_a,
				     matOpts.permutationHML,
				     matOpts.permutationHML);
  densityMatrix_beta.assignFromFull(densityMatrixFull_b,
				    matOpts.permutationHML,
				    matOpts.permutationHML);
  
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "After spinflip: Tr( D_alpha * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix_alpha, S_symm));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "After spinflip: Tr( D_beta  * S ) = %22.11f", (double)symmMatrix::trace_ab(densityMatrix_beta , S_symm));

  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.writeToFile();
  S_symm.writeToFile();

  delete [] densityMatrixFull_tot;
  delete [] densityMatrixFull_spin;
}




void SCF_unrestricted::output_csr_matrices_for_gao()
{
  throw "error: SCF_unrestricted::output_csr_matrices_for_gao not implemented.";
}


void SCF_unrestricted::write_diag_dens_to_file()
{
  int n = basisInfo.noOfBasisFuncs;

  densityMatrix_alpha.readFromFile();
  write_diag_elements_to_file(n, densityMatrix_alpha, "diagdens_alpha.txt",
			      matOpts.permutationHML);
  densityMatrix_alpha.writeToFile();
  
  densityMatrix_beta.readFromFile();
  write_diag_elements_to_file(n, densityMatrix_beta , "diagdens_beta.txt",
			      matOpts.permutationHML);
  densityMatrix_beta.writeToFile();
}


void SCF_unrestricted::save_full_matrices_for_matlab()
{
  throw "error: SCF_unrestricted::save_full_matrices_for_matlab not implemented.";
}


void SCF_unrestricted::report_final_results()
{
  ergo_real S2_exact, S2;
  get_S2(S2_exact, S2);
  do_output(LOG_CAT_RESULTS, LOG_AREA_SCF, "FINAL <S2> = %22.11f", (double)S2);
}


void SCF_unrestricted::get_S2(ergo_real & S2_exact, ergo_real & S2)
{
  // FIXME: the <S2> stuff should work also if N_alpha < N_beta.
  if(noOfElectrons_alpha >= noOfElectrons_beta)
    {
      // Compute <S2>
      S_symm.readFromFile();
      densityMatrix_alpha.readFromFile();
      densityMatrix_beta.readFromFile();
      
      normalMatrix B_alpha;
      B_alpha.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
      B_alpha = (ergo_real)1.0 * densityMatrix_alpha * S_symm;
      
      normalMatrix B_beta;
      B_beta.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
      B_beta = (ergo_real)1.0 * densityMatrix_beta * S_symm;
      
      ergo_real x = normalMatrix::trace_ab(B_alpha, B_beta);
      S2_exact = ( (ergo_real)(noOfElectrons_alpha - noOfElectrons_beta) / 2) * ( ( (ergo_real)(noOfElectrons_alpha - noOfElectrons_beta) / 2) + 1);
      S2 = S2_exact + noOfElectrons_beta - x;
      
      densityMatrix_alpha.writeToFile();
      densityMatrix_beta.writeToFile();
      S_symm.writeToFile();
    }
  else
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "Error: <S2> computation not impl for na < nb");
      S2_exact = 0;
      S2 = 0;
    }
}


void SCF_unrestricted::save_density_as_prevdens()
{
  // alpha
  densityMatrix_alpha.readFromFile();
  Dprev_alpha.readFromFile();
  Dprev_alpha = densityMatrix_alpha;
  densityMatrix_alpha.writeToFile();
  Dprev_alpha.writeToFile();
  // beta
  densityMatrix_beta.readFromFile();
  Dprev_beta.readFromFile();
  Dprev_beta = densityMatrix_beta;
  densityMatrix_beta.writeToFile();
  Dprev_beta.writeToFile();
}



void SCF_unrestricted::report_density_difference()
{
  // alpha
  densityMatrix_alpha.readFromFile();
  Dprev_alpha.readFromFile();
  symmMatrix diff_alpha(densityMatrix_alpha);
  diff_alpha += (ergo_real)-1.0 * Dprev_alpha;
  ergo_real diff_eucl_alpha = GetEuclideanNormOfMatrix(diff_alpha);
  densityMatrix_alpha.writeToFile();
  Dprev_alpha.writeToFile();

  // beta
  densityMatrix_beta.readFromFile();
  Dprev_beta.readFromFile();
  symmMatrix diff_beta(densityMatrix_beta);
  diff_beta += (ergo_real)-1.0 * Dprev_beta;
  ergo_real diff_eucl_beta = GetEuclideanNormOfMatrix(diff_beta);
  densityMatrix_beta.writeToFile();
  Dprev_beta.writeToFile();

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_unrestricted::report_density_difference, diff_eucl_alpha = %22.11f", (double)diff_eucl_alpha);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_unrestricted::report_density_difference, diff_eucl_beta  = %22.11f", (double)diff_eucl_beta );
}

void SCF_unrestricted::compute_dipole_moment()
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_unrestricted::compute_dipole_moment");
  densityMatrix_alpha.readFromFile();
  densityMatrix_beta.readFromFile();
  symmMatrix densityMatrix_tot(densityMatrix_alpha);
  densityMatrix_tot += (ergo_real)1.0 * densityMatrix_beta;
  get_dipole_moment(densityMatrix_tot, basisInfo, matOpts.size_block_info, matOpts.permutationHML, molecule);  
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.writeToFile();
}

void SCF_unrestricted::do_mulliken_pop_stuff()
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_unrestricted::do_mulliken_pop_stuff");
  densityMatrix_alpha.readFromFile();
  densityMatrix_beta.readFromFile();
  S_symm.readFromFile();
  symmMatrix densityMatrix_tot(densityMatrix_alpha);
  densityMatrix_tot += (ergo_real)1.0 * densityMatrix_beta;
  do_mulliken_atomic_charges(densityMatrix_tot,
			     S_symm,
			     basisInfo,
			     matOpts.size_block_info,
			     matOpts.permutationHML,
			     matOpts.inversePermutationHML,
			     molecule);
  densityMatrix_tot.clear(); /* This is to reduce memory usage. */
  symmMatrix spinDensityMatrix(densityMatrix_alpha);
  spinDensityMatrix += (ergo_real)-1.0 * densityMatrix_beta;
  do_mulliken_spin_densities(spinDensityMatrix,
			     S_symm,
			     basisInfo,
			     matOpts.size_block_info,
			     matOpts.permutationHML,
			     matOpts.inversePermutationHML,
			     molecule);
  densityMatrix_alpha.writeToFile();
  densityMatrix_beta.writeToFile();
  S_symm.writeToFile();
}

void SCF_unrestricted::create_mtx_files_F(int const scfIter) {
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx files for Fock matrices");
  { // alpha
    std::stringstream ss_fileName;
    ss_fileName << "F_matrix_alpha_" << scfIter;
    std::stringstream ss_id;
    ss_id << scfopts.calculation_identifier << " - effective Hamiltonian matrix (alpha), SCF cycle " << scfIter;  
    FockMatrix_alpha.readFromFile();
    write_matrix_in_matrix_market_format( FockMatrix_alpha, matOpts.inversePermutationHML, ss_fileName.str(), 
					  ss_id.str(), scfopts.method_and_basis_set );
    FockMatrix_alpha.writeToFile();
  }
  { // beta
    std::stringstream ss_fileName;
    ss_fileName << "F_matrix_beta_" << scfIter;
    std::stringstream ss_id;
    ss_id << scfopts.calculation_identifier << " - effective Hamiltonian matrix (beta), SCF cycle " << scfIter;  
    FockMatrix_beta.readFromFile();
    write_matrix_in_matrix_market_format( FockMatrix_beta, matOpts.inversePermutationHML, ss_fileName.str(), 
					  ss_id.str(), scfopts.method_and_basis_set );
    FockMatrix_beta.writeToFile();
  }  
}
void SCF_unrestricted::create_mtx_files_D(int const scfIter) {
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx files for density matrices");
  { // alpha
    std::stringstream ss_fileName;
    ss_fileName << "D_matrix_alpha_" << scfIter;
    std::stringstream ss_id;
    ss_id << scfopts.calculation_identifier << " - density matrix (alpha), SCF cycle " << scfIter;  
    densityMatrix_alpha.readFromFile();
    write_matrix_in_matrix_market_format( densityMatrix_alpha, matOpts.inversePermutationHML, ss_fileName.str(), 
					  ss_id.str(), scfopts.method_and_basis_set );
    densityMatrix_alpha.writeToFile();
  }
  { // beta
    std::stringstream ss_fileName;
    ss_fileName << "D_matrix_beta_" << scfIter;
    std::stringstream ss_id;
    ss_id << scfopts.calculation_identifier << " - density matrix (beta), SCF cycle " << scfIter;  
    densityMatrix_beta.readFromFile();
    write_matrix_in_matrix_market_format( densityMatrix_beta, matOpts.inversePermutationHML, ss_fileName.str(), 
					  ss_id.str(), scfopts.method_and_basis_set );
    densityMatrix_beta.writeToFile();
  }  
}

void SCF_unrestricted::create_homo_eigvec_file() const
{
  do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Output of HOMO eigenvector not implemented for unrestricted calculations.");
  return;
}
void SCF_unrestricted::create_lumo_eigvec_file() const
{
  do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Output of LUMO eigenvector not implemented for unrestricted calculations.");
  return;
}
void SCF_unrestricted::create_gabedit_file() const
{
  do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Output of gabedit file not implemented for unrestricted calculations.");
  return;
}


void SCF_unrestricted::update_subspace_diff()
{
  throw "Error: SCF_unrestricted::update_subspace_diff() not implemented";
}


void SCF_unrestricted::disturb_fock_matrix(ergo_real subspaceError)
{
  throw "SCF_unrestricted::disturb_fock_matrix() not implemented";
}

void SCF_unrestricted::disturb_dens_matrix(ergo_real subspaceError)
{
  throw "SCF_unrestricted::disturb_dens_matrix() not implemented";
}

void SCF_unrestricted::disturb_dens_matrix_exact(ergo_real subspaceError)
{
  throw "SCF_unrestricted::disturb_dens_matrix_exact() not implemented";
}

void SCF_unrestricted::compute_gradient_fixeddens()
{
  throw "SCF_unrestricted::compute_gradient_fixeddens() not implemented";  
}

