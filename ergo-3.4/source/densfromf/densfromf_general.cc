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

/** @file densfromf_general.cc

    \brief Routine get_dens_from_fock_general() for getting density matrix from a given Fock matrix. 
    This routine calls either get_dens_from_fock_sparse() or get_dens_from_fock_full().

    @author: Elias Rudberg <em>responsible</em>. 
*/
#include "densfromf_general.h"
#include "densfromf_sparse.h"
#include "densfromf_full.h"
#include "densfromf_stochastic.h"
#include "output.h"
#include "utilities.h"


int get_dens_from_fock_general(int n,
			       int noOfOccupiedOrbs,
			       int use_diagonalization,
			       int use_diag_on_error,
			       ergo_real electronicTemperature,
			       symmMatrix & resultDens,
			       ergo_real factor,
			       ergo_real & resultEntropyTerm,
			       symmMatrix & Finput, // written to file
			       intervalType & homoInterval_Finput,
			       intervalType & lumoInterval_Finput,
			       symmMatrix & overlapMatrix,
			       triangMatrix const & invCholFactor, // written to file
			       ergo_real invCholFactor_euclnorm,
			       ergo_real gap_expected_lower_bound, 
			       mat::SizesAndBlocks const & matrixSizesAndBlocks,
			       symmMatrix & F_ort_prev, // written to file
			       intervalType & homoInterval_F_ort_prev,
			       intervalType & lumoInterval_F_ort_prev,
			       ergo_real eigvalueErrorLimit,
			       ergo_real subspaceErrorLimit,
			       mat::normType const truncationNormPurification,
			       int maxMul,
			       int create_m_files,
			       int ignore_purification_failure,
			       int use_rand_perturbation_for_alleigsint,
			       int use_stochastic_orbitals,
			       const std::vector< std::vector<ergo_real> > stochastic_orbitals,
			       std::string stats_prefix,
			       std::map<std::string, double> & puri_stats,
			       int do_sparsity_investigation,
			       int sparsity_plots_resolution_m,
			       int do_comparison_to_simple_purification,
			       int do_puri_mmul_tests,
			       generalVector * eigVecLUMO, 
			       generalVector * eigVecHOMO
			       )
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_dens_from_fock_general, n = %i, use_diagonalization = %i, use_diag_on_error = %i, use_stochastic_orbitals = %i", 
	    n, use_diagonalization, use_diag_on_error, use_stochastic_orbitals);
  resultEntropyTerm = 0; // In nonzero temperature case, this will be set to nonzero value later.
  if(use_stochastic_orbitals == 1) {
    resultDens.readFromFile();
    resultDens.clear();
    if(get_dens_from_fock_stochastic(n, 
				     noOfOccupiedOrbs, 
				     resultDens,
				     factor,
				     Finput,
				     invCholFactor,
				     matrixSizesAndBlocks,
				     stochastic_orbitals) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "get_dens_from_fock_general: error in get_dens_from_fock_stochastic; aborting.");
      return -1;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_dens_from_fock_stochastic finished");
    resultDens.writeToFile();
    return 0;
  }
  int use_diag = 0;
  if(use_diagonalization == 1)
    use_diag = 1;
  else
    {
      // Try purification
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		"calling get_dens_from_fock_sparse, n = %6i, subspaceErrorLimit = %g, eigvalueErrorLimit = %g", 
		n, (double)subspaceErrorLimit, (double)eigvalueErrorLimit);
      if(electronicTemperature != 0)
	throw "Error: (electronicTemperature != 0) not implemented for sparse case.";
      resultDens.readFromFile();
      resultDens.clear();
      if(get_dens_from_fock_sparse(n, 
				   noOfOccupiedOrbs, 
				   resultDens,
				   factor,
				   Finput,
				   homoInterval_Finput,
				   lumoInterval_Finput,
				   invCholFactor,
				   invCholFactor_euclnorm,
				   gap_expected_lower_bound,
				   matrixSizesAndBlocks,
				   F_ort_prev,
				   homoInterval_F_ort_prev,
				   lumoInterval_F_ort_prev,
				   eigvalueErrorLimit,
				   subspaceErrorLimit,
				   truncationNormPurification,
				   maxMul,
				   create_m_files,
				   ignore_purification_failure,
				   use_rand_perturbation_for_alleigsint, 
				   stats_prefix,
				   puri_stats,
				   do_sparsity_investigation,
				   sparsity_plots_resolution_m,
				   do_comparison_to_simple_purification,
				   do_puri_mmul_tests,
				   eigVecLUMO,
				   eigVecHOMO) != 0)
	{
	  if(use_diag_on_error)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "get_dens_from_fock_general: error in get_dens_from_fock_sparse; trying with diagonalization instead.");
	      use_diag = 1;
	    }
	  else
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "get_dens_from_fock_general: error in get_dens_from_fock_sparse; aborting.");
	      return -1;
	    }
	}
      else
	{
	  // Purification success!
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "get_dens_from_fock_sparse finished OK.");
	}
      resultDens.writeToFile();
    }
  if(use_diag == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling get_dens_from_fock, n = %i", n);
      
      std::vector<ergo_real> F_full(n*n);
      std::vector<ergo_real> S_full(n*n);      

      {
	// Create full matrix versions of F and S
	normalMatrix* tmpMat;
	Finput.readFromFile();
	tmpMat = new normalMatrix(Finput);
	Finput.writeToFile();
	tmpMat->fullMatrix(F_full);
	delete tmpMat;
	overlapMatrix.readFromFile();
	tmpMat = new normalMatrix(overlapMatrix);
	overlapMatrix.writeToFile();
	tmpMat->fullMatrix(S_full);
	delete tmpMat;
      }
      
      std::vector<ergo_real> densityMatrixFull(n*n);
      std::vector<ergo_real> eigVecLUMO_tmp(n);
      std::vector<ergo_real> eigVecHOMO_tmp(n);


      if(get_dens_from_fock_full(n, 
				 noOfOccupiedOrbs, 
				 &densityMatrixFull[0], 
				 &F_full[0], 
				 &S_full[0], 
				 factor,
				 electronicTemperature,
				 resultEntropyTerm,
				 &eigVecLUMO_tmp[0],
				 &eigVecHOMO_tmp[0]) != 0)
	{
	  throw "error in get_dens_from_fock_full";
	}
      
      resultDens.readFromFile();
      resultDens.assignFromFull(densityMatrixFull);
      resultDens.writeToFile();
      if (eigVecLUMO) 
	eigVecLUMO->assign_from_full(eigVecLUMO_tmp, matrixSizesAndBlocks);
      if (eigVecHOMO)
	eigVecHOMO->assign_from_full(eigVecHOMO_tmp, matrixSizesAndBlocks);
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_dens_from_fock finished");
    }

  return 0;
}

