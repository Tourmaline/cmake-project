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

/* 
   Projection of electron density from one basis set to another.
*/

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>

#include "density_description_file.h"
#include "density_projection.h"
#include "densfromf_full.h"
#include "densfromf_sparse.h"
#include "densfromf_general.h"
#include "integrals_general.h"
#include "operator_matrix.h"
#include "matrix_algebra.h"
#include "memorymanag.h"
#include "output.h"
#include "utilities.h"
#include "matrix_utilities.h"


#if 0
static void
printmatr(ergo_real* M, int n1, int n2, const char* s)
{
  printf("matrix '%s':\n", s);
  int i, j;
  for(i = 0; i < n1; i++)
    {
      for(j = 0; j < n2; j++)
	printf("%8.3f ", M[i*n2+j]);
      printf("\n");
    }
}

static void
printmatr2(ergo_real* M, int n1, int n2, const char* s)
{
  printf("matrix '%s':\n", s);
  int i, j;
  for(i = 0; i < n1; i++)
    {
      for(j = 0; j < n2; j++)
	printf("%8.3f ", M[j*n1+i]);
      printf("\n");
    }
}

static void
printmatrsparse(int n, int m, const normalMatrix & M_sparse, const char* s)
{
  ergo_real* M = new ergo_real[n*m];
  M_sparse.fullmatrix(M, n, m);
  printmatr2(M, n, m, s);
  delete M;
}
#endif


int
load_density_and_project_full(const char *densityFileName,
			      int noOfDensityMatrices,
			      const IntegralInfo* integralInfo,
			      const BasisInfoStruct & basisInfo,
			      ergo_real** densityMatrixList,
			      int do_purification,
			      const int* noOfElectronsList,
			      ergo_real electronic_temperature)
{
  int n = basisInfo.noOfBasisFuncs;
  ergo_real* densityMatrixListForStartingGuess[2];
  BasisInfoStruct* basisInfoStartingGuess = NULL;
    
  if(ddf_load_density(densityFileName, noOfDensityMatrices,
                      *integralInfo, &basisInfoStartingGuess,
                      densityMatrixListForStartingGuess) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density");
      return -1;
    }

  /* Projection part, could be a separate routine, too. */

  // get matrix R : overlap of main basis and startingguess basis
  int n_sg = basisInfoStartingGuess->noOfBasisFuncs;
  ergo_real* R = ergo_new(n*n_sg, ergo_real);
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "calling compute_overlap_matrix for R");
  if(compute_overlap_matrix(*basisInfoStartingGuess, basisInfo, R) != 0)
  {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_overlap_matrix for matrix R");
      return -1;
  }

  // get matrix S : main overlap matrix
  ergo_real* S = ergo_new(n*n, ergo_real);
  if(compute_overlap_matrix(basisInfo, basisInfo, S) != 0)
  {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_overlap_matrix for matrix S");
      return -1;
  }

  // get matrix Sinv : inverse of main overlap matrix
  //  ergo_real* Sinv = ergo_new(n*n, ergo_real);

  
  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: get_inverse_of_posdef_symm_matrix not implemented.");
  return -1;
  /*
  if(get_inverse_of_posdef_symm_matrix(n, S, Sinv) != 0)
  {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_inverse_of_posdef_symm_matrix");
      return -1;
  }
  */
  /* ELIAS NOTE 2011-09-15: 
     Remainder of this routine removed since it could anyway not be used since get_inverse_of_posdef_symm_matrix was not implemented. */
}




static int
compute_R_matrix_sparse(const BasisInfoStruct & basisInfo_A,
			const BasisInfoStruct & basisInfo_B,
			normalMatrix & result_R,
			ergo_real sparse_threshold,
                        std::vector<int> const & matrixPermutationVec,
                        std::vector<int> const & matrixPermutationVec_sg)
{
  int n_A = basisInfo_A.noOfBasisFuncs;
  int n_B = basisInfo_B.noOfBasisFuncs;

  std::vector<int> nvaluesList(n_A);
  std::vector< std::vector<int> > colindList(n_A);
  std::vector< std::vector<ergo_real> > valuesList(n_A);
  
  if(compute_operator_matrix_sparse(basisInfo_A, 
				    basisInfo_B,
				    0,
				    0,
				    0,
				    n_A,
				    n_B,
				    nvaluesList,
				    colindList,
				    valuesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_operator_matrix_sparse");
      return -1;
    }
  
  // Now convert result to three vectors to prepare for HML storage.
  long nvalues = 0;
  for(int i = 0; i < n_A; i++)
    nvalues += nvaluesList[i];
  // allocate vectors
  std::vector<int> rowind(nvalues);
  std::vector<int> colind(nvalues);
  std::vector<ergo_real> values(nvalues);
  // populate vectors
  long count = 0;
  for(long i = 0; i < n_A; i++)
    {
      for(long j = 0; j < nvaluesList[i]; j++)
	{
	  rowind[count] = i;
	  colind[count] = colindList[i][j];
	  values[count] = valuesList[i][j];
	  count++;
	} // END FOR j
    } // END FOR i

  // Now the information is in rowind colind values.

  /* FIXME: there is confusion about the meaning of "rows" and "columns".
     Here it is reversed, otherwise it did not work.
     We should use the same convention for rows and columns everywhere. */
  result_R.assign_from_sparse(colind,
			      rowind,
			      values,
			      matrixPermutationVec_sg,
                              matrixPermutationVec);
  
  result_R.eucl_thresh(sparse_threshold);

  // Write to file and read again to reduce memory fragmentation.
  result_R.writeToFile();
  result_R.readFromFile();

  return 0;
}




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
				ergo_real electronic_temperature)
{
  Util::TimeMeter timeMeter;
  int n = basisInfo.noOfBasisFuncs;
  BasisInfoStruct* basisInfoStartingGuess = NULL;
  
  long nvaluesList[2];
  int* rowindList[2];
  int* colindList[2];
  ergo_real* valuesList[2];
  // Call ddf_load_density_sparse. It will allocate rowindList, colindList, valuesList,
  // which we will have to delete afterwards.
  int noOfDensitiesRead = 0;
  if(ddf_load_density_sparse(densityFileName,
			     *integralInfo, &basisInfoStartingGuess,
			     &noOfDensitiesRead,
			     rowindList,
			     colindList,
			     valuesList,
			     nvaluesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density_sparse");
      return -1;
    }

  int n_sg = basisInfoStartingGuess->noOfBasisFuncs;

  // Get permutation object for starting guess basis set.
  int sparse_block_size = 20;
  int blockSizeFactor = 8;
  mat::SizesAndBlocks matrix_size_block_info_sg 
    = prepareMatrixSizesAndBlocks(n_sg,
                                  sparse_block_size,
                                  blockSizeFactor, 
                                  blockSizeFactor, 
                                  blockSizeFactor);
  std::vector<int> matrixPermutationVec_sg;
  getMatrixPermutation(*basisInfoStartingGuess,
                       sparse_block_size,
                       blockSizeFactor, blockSizeFactor, blockSizeFactor,
                       matrixPermutationVec_sg);
  
  // Now we have one or two matrices in memory stored as vectors. Convert to symmMatrix form.
  symmMatrix dens_list_1_sparse[2];
  int maxDensities = noOfDensityMatrices > noOfDensitiesRead
    ? noOfDensityMatrices : noOfDensitiesRead;
  
  for(int i = 0; i < maxDensities; i++)
    dens_list_1_sparse[i].resetSizesAndBlocks
      (matrix_size_block_info_sg, matrix_size_block_info_sg);

  for(int i = 0; i < noOfDensitiesRead; i++)
    {
      std::vector<int> rowIndTmp(rowindList[i], rowindList[i] + nvaluesList[i]);
      std::vector<int> colIndTmp(colindList[i], colindList[i] + nvaluesList[i]);
      std::vector<ergo_real> valuesTmp(valuesList[i], valuesList[i] + nvaluesList[i]);
      dens_list_1_sparse[i].assign_from_sparse(rowIndTmp, colIndTmp, valuesTmp,
                                               matrixPermutationVec_sg,
                                               matrixPermutationVec_sg);
      dens_list_1_sparse[i].eucl_thresh(sparse_threshold);
      delete []rowindList[i];
      delete []colindList[i];
      delete []valuesList[i];
    }  

  /* Handle conversions between restricted and unrestricted. */
  if(noOfDensitiesRead == 1 && noOfDensityMatrices == 2) {
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "converting restricted density into unrestricted");
    dens_list_1_sparse[1] = dens_list_1_sparse[0];
  } else if(noOfDensitiesRead == 2 && noOfDensityMatrices == 1) {
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "converting unrestricted density into restricted");
    dens_list_1_sparse[0] += dens_list_1_sparse[1];
  }

  // Reduce memory fragmentation by writing to file and reading in again.
  for(int i = 0; i < noOfDensityMatrices; i++)
    dens_list_1_sparse[i].writeToFile();
  for(int i = 0; i < noOfDensityMatrices; i++)
    dens_list_1_sparse[i].readFromFile();


  /* Projection part, could be a separate routine, too. */

  // get matrix R : overlap of main basis and startingguess basis
  output_current_memory_usage(LOG_AREA_SCF, "Before getting compute_R_matrix_sparse");

  normalMatrix R_sparse;
  R_sparse.resetSizesAndBlocks
    (matrix_size_block_info_sg, matrix_size_block_info);
  
  if(compute_R_matrix_sparse(basisInfo,
			     *basisInfoStartingGuess,
			     R_sparse,
			     sparse_threshold,
                             matrixPermutationVec,
                             matrixPermutationVec_sg) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_R_matrix_sparse");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_SCF, "After compute_R_matrix_sparse");

  // Now calculate RT * P * R for each density matrix
  // FIXME: implement BT * P * B in matrix library. B is normalMatrix, P is symmMatrix. Result is symm.
  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      normalMatrix PR_sparse;
      PR_sparse.resetSizesAndBlocks
        (matrix_size_block_info_sg, matrix_size_block_info);

      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "trying to create PR_sparse");

      Util::TimeMeter timeMeterPR;
      PR_sparse = (ergo_real)1.0 * dens_list_1_sparse[i] * R_sparse;
      timeMeterPR.print(LOG_AREA_MAIN, "PR_sparse = P * R multiplication");

      PR_sparse.eucl_thresh(sparse_threshold);

      output_current_memory_usage(LOG_AREA_MAIN, "After creating PR_sparse");

      dens_list_1_sparse[i].clear();

      normalMatrix RT;
      RT.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info_sg);
      RT = transpose(R_sparse);
      
      normalMatrix RT_P_R;
      RT_P_R.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info);
      
      RT_P_R = (ergo_real)1.0 * RT * PR_sparse;

      output_current_memory_usage(LOG_AREA_MAIN, "After creating RT_P_R");
      RT.clear();
      PR_sparse.clear();

      RT_P_R.eucl_thresh(sparse_threshold);
      *densityMatrixList[i] = RT_P_R;
      densityMatrixList[i]->eucl_thresh(sparse_threshold);
    }

  // Now we do not need R_sparse any more.
  R_sparse.clear();

  // Projection done. Now do purification to force idempotency and correct trace.
  // densityMatrixList now contains RT * P * R for each density matrix.

  output_current_memory_usage(2, "While doing projection of starting guess (2).");

  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      symmMatrix SDS;
      SDS.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info);
      SDS = *densityMatrixList[i];
      SDS *= -1.0;

      output_current_memory_usage(LOG_AREA_MAIN, "After creating matrix -SDS");
	  
      ergo_real factor = 2;
      if(noOfDensityMatrices == 2)
	factor = 1;
      if(noOfElectronsList == NULL)
	return -1;

      symmMatrix F_ort_prev_dummy;
      F_ort_prev_dummy.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info);

      intervalType homoInterval_dummy1(-1e22,1e22);
      intervalType lumoInterval_dummy1(-1e22,1e22);      
      intervalType homoInterval_dummy2(-1e22,1e22);
      intervalType lumoInterval_dummy2(-1e22,1e22);      

      int noOfOccupiedOrbs = noOfElectronsList[i];
      if(noOfDensityMatrices == 1)
	noOfOccupiedOrbs /= 2;
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
		"calling get_dens_from_fock_general for SDS to force idempotency and correct trace of starting guess");

      // The result of purification will be placed in densityMatrixList[i], but it is not used as input,
      // so we can clear it now to free up some memory.
      densityMatrixList[i]->clear();
      densityMatrixList[i]->writeToFile();

      F_ort_prev_dummy.writeToFile();
      SDS.writeToFile();

      std::map<std::string, double> puri_stats;
      ergo_real electronicEntropyTerm_dummy = 0;
      std::vector< std::vector<ergo_real> > stochasticOrbsDummy;
      if(get_dens_from_fock_general(n,
				    noOfOccupiedOrbs,
				    use_diagonalization,
				    use_diag_on_error,
				    electronic_temperature,
				    *densityMatrixList[i],
				    factor,
				    electronicEntropyTerm_dummy,
				    SDS, 
				    homoInterval_dummy1,
				    lumoInterval_dummy1,
				    S_symm,
				    invCholFactor,
				    invCholFactor_euclnorm,
				    gap_expected_lower_bound,
				    matrix_size_block_info,
				    F_ort_prev_dummy,
				    homoInterval_dummy2,
				    lumoInterval_dummy2,
				    purification_eigvalue_err_limit,
				    purification_subspace_err_limit,
				    purification_truncation_norm,
				    purification_maxmul,
				    purification_create_m_files,
				    purification_ignore_failure,
				    purification_use_rand_perturbation_for_alleigsint,
				    0,
				    stochasticOrbsDummy,
				    "",
				    puri_stats,
				    0, // no sparsity investigation
				    1, // dummy 1 here (no investigation)
				    0, // no comparison to simple purification
				    0  // purification mmul tests
				    ) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_dens_from_fock_general.");
	  return -1;
	}

      output_current_memory_usage(LOG_AREA_MAIN, "After get_dens_from_fock_general");
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "get_dens_from_fock_general finished OK.");
    } // END FOR i
  delete basisInfoStartingGuess;

  timeMeter.print(LOG_AREA_MAIN, "load_density_and_project_sparse");

  return 0;
}


