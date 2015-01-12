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

/** @file densfromf_sparse.cc

    \brief Routine get_dens_from_fock_sparse() for getting density matrix from a given Fock matrix using purification.

    @author: Elias Rudberg <em>responsible</em>. 
*/
#include "densfromf_sparse.h"
#include "output.h"
#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include "utilities.h"
#include "matrix_utilities.h"
#include "TC2.h"
#include "PuriInfo.h"
#include "Purification.h"
#include "Purification_scaled.h"
#include "units.h"
#include "machine_epsilon.h"
#include "simple_purification.h"
#include "puri_mmul_tests.h"
#include "AllocatorManager.h"

#define RUN_NEW_PURIFICATION 0

typedef mat::DebugLevelLow debugPolicy;

typedef mat::PuriInfo<ergo_real, generalVector, debugPolicy> puriInfoType;

typedef mat::Purification<ergo_real, symmMatrix, debugPolicy> purificationType;

typedef mat::Interval<ergo_real> intervalType;


static int globalCounter = 0; // For naming purification statistics matlab files. 
                              // Now also for naming m-files from comparison with simple purification.


intervalType getAllEigsInterval(int n,
				symmMatrix & F, // not const because we want to write+read it
				mat::SizesAndBlocks const & matrixSizesAndBlocks,
				int use_rand_perturbation)
{
  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	    "getAllEigsInterval, use_rand_perturbation = %i", 
	    use_rand_perturbation);

  ergo_real lambdaMinGers, lambdaMaxGers;
  output_current_memory_usage(LOG_AREA_DENSFROMF, "In getAllEigsInterval, before F.gersgorin().");
  F.gersgorin(lambdaMinGers, lambdaMaxGers);
  output_current_memory_usage(LOG_AREA_DENSFROMF, "In getAllEigsInterval, after  F.gersgorin().");
  
  symmMatrix id;
  id.resetSizesAndBlocks(matrixSizesAndBlocks, matrixSizesAndBlocks);
  id = 1;

  // Use low accuracy here, otherwise this may take a long time.
  ergo_real acc = std::sqrt(std::sqrt(get_machine_epsilon()));

  // Write F to file before creating copy, to reduce memory usage.
  F.writeToFile();
  symmMatrix Fshifted(F); // Copy on file
  Fshifted.readFromFile();
  Fshifted += (ergo_real)(-1.0) * lambdaMinGers * id;

  output_current_memory_usage(LOG_AREA_DENSFROMF, "In getAllEigsInterval,  after creating Fshifted");

  ergo_real lambdaMax;
  int maxIter = 100;
  try {
    if(use_rand_perturbation)
      {
	add_random_diag_perturbation(n, Fshifted, 0.5 * acc);
	ergo_real newAcc = 0.5 * acc;
	lambdaMax = Fshifted.eucl(acc, maxIter) + lambdaMinGers + newAcc;
      }
    else
      lambdaMax = Fshifted.eucl(acc, maxIter) + lambdaMinGers + acc;
  }
  catch (mat::AcceptableMaxIter & e) {
    do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, 
	      "getAllEigsInterval, Lanczos failed to find extreme upper eigenvalue within maxiter, "
	      "using Gersgorin bound");
    lambdaMax = lambdaMaxGers;
  }

  /* Now we want to create Fshifted = ( F - lambdaMaxGers*id ) but we
     do this starting from the existing Fshifted, correcting it back
     to F and then subtracting lambdaMaxGers*id. */
  Fshifted += (ergo_real)( 1.0) * lambdaMinGers * id; // Now Fshifted = F.
  Fshifted += (ergo_real)(-1.0) * lambdaMaxGers * id;
  ergo_real lambdaMin;
  try {
    if(use_rand_perturbation)
      {
	add_random_diag_perturbation(n, Fshifted, 0.5 * acc);
	ergo_real newAcc = 0.5 * acc;
	lambdaMin = -Fshifted.eucl(newAcc, maxIter) + lambdaMaxGers - newAcc;
      }
    else
      lambdaMin = -Fshifted.eucl(acc, maxIter) + lambdaMaxGers - acc;
  }
  catch (mat::AcceptableMaxIter & e) {
    do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, 
	      "getAllEigsInterval, Lanczos failed to find extreme lower eigenvalue within maxiter, "
	      "using Gersgorin bound");
    lambdaMin = lambdaMinGers;
  }

  intervalType allEigsInterval(lambdaMin, lambdaMax);

  Fshifted.clear();
  F.readFromFile();

  timeMeter.print(LOG_AREA_DENSFROMF, "getAllEigsInterval");

  return allEigsInterval;
}


int get_dens_from_fock_sparse(int n, 
			      int noOfOccupiedOrbs, 
			      symmMatrix & resultDens, 
			      ergo_real factor,
			      symmMatrix const & Finput, // written to file
			      intervalType & homoInterval_Finput,
			      intervalType & lumoInterval_Finput,
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
			      std::string stats_prefix,
			      std::map<std::string, double> & puri_stats,
			      int do_sparsity_investigation,
			      int sparsity_plots_resolution_m,
			      int do_comparison_to_simple_purification,
			      int do_puri_mmul_tests,
			      generalVector * eigVecLUMO,
			      generalVector * eigVecHOMO)
{
  Util::TimeMeter timeMeterTot;

  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse() start!");

  std::string allocStatsStr1 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Before writeAndReadAll(): %s", allocStatsStr1.c_str());

  Util::TimeMeter timeMeterWriteAndReadAll;
  std::string sizesStr = mat::FileWritable::writeAndReadAll();
  timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());
  
  std::string allocStatsStr2 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "After writeAndReadAll(): %s", allocStatsStr2.c_str());

  if(noOfOccupiedOrbs == 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse: (noOfOccupiedOrbs == 0), skipping.");
      resultDens.clear();
      return 0;
    }

  // Select tolerated errors in the occupied subspace for the three truncations
  // and for purification.
  ergo_real subspaceThr_1     = 0.1 * subspaceErrorLimit;
  ergo_real subspaceThr_Puri  = 0.7 * subspaceErrorLimit;
  ergo_real subspaceThr_2     = 0.1 * subspaceErrorLimit;
  ergo_real subspaceThr_3     = 0.1 * subspaceErrorLimit;

  // Select tolerated errors in eigenvalues 
  ergo_real eigvalueThr_Puri  = 0.7  * eigvalueErrorLimit;
  ergo_real eigvalueThr_2     = 0.15 * eigvalueErrorLimit;
  ergo_real eigvalueThr_3     = 0.15 * eigvalueErrorLimit;
 
  symmMatrix F(Finput);
  F.readFromFile();
  std::string allocStatsStr3 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "After F.readFromFile(): %s", allocStatsStr3.c_str());
 
  {
    triangMatrix invCholFactor_tmp(invCholFactor);
    invCholFactor_tmp.readFromFile();
    output_current_memory_usage(LOG_AREA_DENSFROMF, "In get_dens_from_fock_sparse, before F = tr(Z) * F * Z");
    Util::TimeMeter timeMeterFortTransf;
    F = transpose(invCholFactor_tmp) * F * invCholFactor_tmp;
    timeMeterFortTransf.print(LOG_AREA_DENSFROMF, "F = transpose(invCholFactor) * F * invCholFactor");
    output_current_memory_usage(LOG_AREA_DENSFROMF, "In get_dens_from_fock_sparse,  after F = tr(Z) * F * Z");
  } // invCholFactor_tmp goes out of scope here

  // Now F contains F_ort. 

  if ( do_sparsity_investigation == 1 ) {
    output_magnitude_histogram( F, stats_prefix + "mag_histogram_F_ort", sparsity_plots_resolution_m);
  }

  F_ort_prev.readFromFile();
  output_current_memory_usage(LOG_AREA_DENSFROMF, 
			      "In get_dens_from_fock_sparse,  after F_ort_prev.readFromFile()");

  // Compare F_ort to F_ort_prev to check how far eigenvalues can have moved.
  {
    ergo_real maxEigValMovement_frob = symmMatrix::frob_diff(F, F_ort_prev);
    output_current_memory_usage(LOG_AREA_DENSFROMF, 
				"In get_dens_from_fock_sparse,  after getting maxEigValMovement_frob ");
    ergo_real acc = std::sqrt(get_machine_epsilon());
    Util::TimeMeter timeMeterMixedDiff;
    ergo_real maxEigValMovement_mixed = symmMatrix::mixed_diff(F, F_ort_prev, acc) + acc;
    timeMeterMixedDiff.print(LOG_AREA_DENSFROMF, "symmMatrix::mixed_diff for maxEigValMovement_mixed");
    output_current_memory_usage(LOG_AREA_DENSFROMF, 
				"In get_dens_from_fock_sparse,  after getting maxEigValMovement_mixed");
    Util::TimeMeter timeMeterEuclDiff;
    ergo_real maxEigValMovement_eucl = symmMatrix::eucl_diff(F, F_ort_prev, acc) + acc;
    timeMeterEuclDiff.print(LOG_AREA_DENSFROMF,  "symmMatrix::eucl_diff  for maxEigValMovement_eucl ");
    output_current_memory_usage(LOG_AREA_DENSFROMF, 
				"In get_dens_from_fock_sparse,  after getting maxEigValMovement_eucl ");
    
    // Increase HOMO/LUMO intervals so that they for sure contain the HOMO and LUMO eigenvalues of F_ort
    homoInterval_F_ort_prev.increase(maxEigValMovement_eucl);
    lumoInterval_F_ort_prev.increase(maxEigValMovement_eucl);
    
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxEigValMovement_frob  = %22.11f", (double)maxEigValMovement_frob);
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxEigValMovement_mixed = %22.11f", (double)maxEigValMovement_mixed);
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxEigValMovement_eucl  = %22.11f", (double)maxEigValMovement_eucl);
  }
  
  // Now, we will truncate F_ort (and update eigenvalue intervals):
  ergo_real truncError_1;
  {
    ergo_real gapMin = lumoInterval_F_ort_prev.low() - homoInterval_F_ort_prev.upp();
    ergo_real gapMax = lumoInterval_F_ort_prev.upp() - homoInterval_F_ort_prev.low();
    ergo_real threshold_1;
    // We consider the gap to be accurately known if the uncertainty is at most 10 %
    if ( gapMin > 0 && (gapMax - gapMin) / gapMin < 0.1 ) 
      // Gap is accurately known: we use gapMin
      threshold_1 = subspaceThr_1 * gapMin / (1+subspaceThr_1);
    else 
      // Gap is not accurately known. To avoid choosing a very tight
      // threshold value due to a small lower bound for the gap, we
      // use the largest of 'gap_expected_lower_bound' and calculated
      // 'gapMin':
      threshold_1 = gapMin > gap_expected_lower_bound ?
	subspaceThr_1 * gapMin / (1+subspaceThr_1) :
	subspaceThr_1 * gap_expected_lower_bound / (1+subspaceThr_1);
    
    double nnzF_before_trunc_pc  = (double)F.nnz()  * 100 / ((double)n*n); 
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Truncating F_ort ( %s ), selected threshold = %10.6g", 
	      mat::getNormTypeString(truncationNormPurification).c_str(), (double)threshold_1);
    Util::TimeMeter timeMeterFThresh;
    truncError_1 = F.thresh( threshold_1, truncationNormPurification );
    double nnzF_after_trunc_pc  = (double)F.nnz()  * 100 / ((double)n*n); 
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	      "Truncated F_ort ( %s ), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
	      mat::getNormTypeString(truncationNormPurification).c_str(),(double)threshold_1, (double)truncError_1, nnzF_before_trunc_pc, nnzF_after_trunc_pc);
    timeMeterFThresh.print(LOG_AREA_DENSFROMF, "Truncation of F_ort");
    puri_stats[stats_prefix + "nnz_percentage_F_ort"]   = nnzF_after_trunc_pc;
    
    
    // Increase HOMO and LUMO intervals so that they contain the eigenvalues of the truncated matrix:
    homoInterval_F_ort_prev.increase( truncError_1 );
    lumoInterval_F_ort_prev.increase( truncError_1 );
  }
  
  // Get interval containing all eigenvalues of F_ort.
  // We do this after truncation since getAllEigsInterval uses lots of memory (the memory usage can be reduced if needed).
  intervalType allEigsInterval = getAllEigsInterval(n, F, matrixSizesAndBlocks, use_rand_perturbation_for_alleigsint);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	    "allEigsInterval : [ %17.12f %17.12f ]",
	    (double)allEigsInterval.low(), (double)allEigsInterval.upp());


  
  // Now overwrite F_ort_prev with the new F_ort.
  F_ort_prev = F;
  F_ort_prev.writeToFile();
  // The HOMO and LUMO intervals now contain the HOMO and LUMO
  // eigenvalues of F_ort_prev but improved values will hopefully be
  // calculated in purification.
  
#if RUN_NEW_PURIFICATION
  /*******    CODE FOR NEW PURIFICATION     *******/
  intervalType homoIntervalSaved = homoInterval_F_ort_prev;
  intervalType lumoIntervalSaved = lumoInterval_F_ort_prev;
  /******* END OF CODE FOR NEW PURIFICATION *******/
#endif  

  
  if(noOfOccupiedOrbs == n)
    {
      // Special case: all occupied!
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse: (noOfOccupiedOrbs == n), setting resultDens = 1.");
      resultDens = 1;
    }
  else
    {
      // purification
      
      homoInterval_F_ort_prev.intersect(allEigsInterval);
      lumoInterval_F_ort_prev.intersect(allEigsInterval);
      
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"creating puriInfo object with "
		"n = %6i, nocc = %6i, maxMul = %3i, "
		"subspaceErrorLimit = %6g, "
		"eigvalueErrorLimit = %6g", 
		n, noOfOccupiedOrbs, maxMul, 
		(double)subspaceThr_Puri, (double)eigvalueThr_Puri);

      puriInfoType puriInfo(n,
			    noOfOccupiedOrbs,
			    allEigsInterval,
			    homoInterval_F_ort_prev,
			    lumoInterval_F_ort_prev,
			    eigvalueThr_Puri,
			    subspaceThr_Puri,
			    truncationNormPurification,
			    maxMul);
      
      mat::Gblas::timekeeping = true;
      mat::Gblas::time = 0;
    
      /* EMANUEL COMMENT: 
	 - truncationNorm can be set to any of 
	   mat::frobNorm, mat::mixedNorm, or mat::euclNorm
	   The best choice depends on a trade-off between spending
	   time in truncation and in matrix-matrix multiplication.
	 - XmX2Norm should probably be mat::euclNorm in order for 
	   purification to converge.
       */
      //      mat::normType const truncationNorm = mat::euclNorm;
      mat::normType const XmX2Norm       = mat::euclNorm;      
      purificationType purification(F, XmX2Norm, puriInfo);

      output_current_memory_usage(LOG_AREA_DENSFROMF, "Before purification.purify()");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"calling purification.purify(), number of threads = %i, trunc norm '%s'",
		mat::Params::getNProcs(), mat::getNormTypeString(truncationNormPurification).c_str());

      mat::FileWritable::resetStats();
      time_t puriStartWallTime;
      time(&puriStartWallTime);

      Util::TimeMeter timeMeterPurification;
      purification.purify();
      timeMeterPurification.print(LOG_AREA_DENSFROMF, "purification.purify()");

      {
	std::stringstream ss;
	ss << "Accumulated wall times for writeToFile in purify()          : " << mat::FileWritable::getStatsTimeWrite();
	do_output( LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, ss.str().c_str() );
      }
      {
	std::stringstream ss;
	ss << "Accumulated wall times for readFromFile in purify()         : " << mat::FileWritable::getStatsTimeRead();
	do_output( LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, ss.str().c_str() );
      }
      {
	std::stringstream ss;
	ss << "Accumulated wall times for copy and assign in purify()      : " << mat::FileWritable::getStatsTimeCopyAndAssign();
	do_output( LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, ss.str().c_str() );
      }


      {
	std::stringstream ss;
	ss << "Number of calls to writeToFile in purify()                  : " << mat::FileWritable::getStatsCountWrite();
	do_output( LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str() );
      }
      {
	std::stringstream ss;
	ss << "Number of calls to readFromFile in purify()                 : " << mat::FileWritable::getStatsCountRead();
	do_output( LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str() );
      }
      {
	std::stringstream ss;
	ss << "Number of calls to FileWritable copy and assign in purify() : " << mat::FileWritable::getStatsCountCopyAndAssign();
	do_output( LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str() );
      }


      do_output(LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, 
		"mat::Gblas::time after purification : %12.6f", 
		(double)mat::Gblas::time);

      int nMulTot = puriInfo.getNSteps();
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Purification finished, %3i multiplications.", nMulTot);

      output_current_memory_usage(LOG_AREA_DENSFROMF, "After  purification.purify()");

      if(!puriInfo.converged())
	{
	  if(!ignore_purification_failure)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, 
			"Error in purification: puriInfo.converged() "
			"returned false.");
	      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, 
			"Outputing info from %3i purification iterations.", 
			nMulTot);
	      for(int i = 0; i < nMulTot; i++)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "iter %3i: TraceX: %22.11f  TraceX2: %22.11f  XmX2EuclNorm_mid: %22.11f", 
			    i, (double)puriInfo(i).getTraceX(), (double)puriInfo(i).getTraceX2(), (double)puriInfo(i).getXmX2EuclNorm().midPoint());
		  do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "HOMO low %22.15f upp %22.15f  LUMO low %22.15f upp %22.15f", 
			    (double)puriInfo(i).getHomo().low(), (double)puriInfo(i).getHomo().upp(), 
			    (double)puriInfo(i).getLumo().low(), (double)puriInfo(i).getLumo().upp());
		}
	      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "End of purification error info.");
	      // Before returning, make sure matrices are read/written to file as expected by caller.
	      return -1;
	    }
	  else
	    {
	      do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, "WARNING: purification NOT converged, ignoring.");
	    }
	}
      else
	{      
	  if ( puriInfo. correct_occupation_was_forced() )
	    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		      "Warning: Correct occupation count could not be guaranteed in purification.");
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		    "Purification converged OK, subspaceError <= %22.11f", 
		    (double)puriInfo.subspaceError());
	}

      // Get some statistics from puriInfo object
      // use double instead of int to avoid integer overflow.
      double nnz_X_sum = 0;
      double nnz_X_min = (double)n*n+1;
      double nnz_X_max = 0;
      int    nnz_X_min_iter = -1;
      int    nnz_X_max_iter = -1;
      double nnz_X2_sum = 0;
      double nnz_X2_min = (double)n*n+1;
      double nnz_X2_max = 0;
      int    nnz_X2_min_iter = -1;
      int    nnz_X2_max_iter = -1;
      for(int i = 0; i < nMulTot; i++)
	{
	  /* Note: use size_t instead of int here, to avoid overflow. */
	  size_t nnz_X  = puriInfo(i).getNnzX();
	  size_t nnz_X2 = puriInfo(i).getNnzX2();
	  nnz_X_sum  += nnz_X;
	  nnz_X2_sum += nnz_X2;
          if(nnz_X < nnz_X_min)
            {
              nnz_X_min = nnz_X;
              nnz_X_min_iter = i;
            }
          if(nnz_X > nnz_X_max)
            {
              nnz_X_max = nnz_X;
              nnz_X_max_iter = i;
            }
          if(nnz_X2 < nnz_X2_min)
            {
              nnz_X2_min = nnz_X2;
              nnz_X2_min_iter = i;
            }
          if(nnz_X2 > nnz_X2_max)
            {
              nnz_X2_max = nnz_X2;
              nnz_X2_max_iter = i;
            }
	}
      double nnz_X_avg  = nnz_X_sum  / nMulTot;
      double nnz_X2_avg = nnz_X2_sum / nMulTot;
      double X_min_pc  = (double)nnz_X_min  * 100 / ((double)n*n); 
      double X_max_pc  = (double)nnz_X_max  * 100 / ((double)n*n);
      double X_avg_pc  = (double)nnz_X_avg  * 100 / ((double)n*n);
      double X2_min_pc = (double)nnz_X2_min * 100 / ((double)n*n);
      double X2_max_pc = (double)nnz_X2_max * 100 / ((double)n*n);
      double X2_avg_pc = (double)nnz_X2_avg * 100 / ((double)n*n);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                "X :  min %5.1f %%  max %5.1f %%  avg %5.1f %%  miniter %3i  maxiter %3i",
                X_min_pc, X_max_pc, X_avg_pc, nnz_X_min_iter, nnz_X_max_iter);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                "X2:  min %5.1f %%  max %5.1f %%  avg %5.1f %%  miniter %3i  maxiter %3i",
                X2_min_pc, X2_max_pc, X2_avg_pc, nnz_X2_min_iter, nnz_X2_max_iter);

      
      puri_stats[stats_prefix + "puri_acc_time_square"]   = puriInfo.getAccumulatedTimeSquare();
      puri_stats[stats_prefix + "puri_acc_time_thresh"]   = puriInfo.getAccumulatedTimeThresh();
      puri_stats[stats_prefix + "puri_acc_time_xmx2norm"] = puriInfo.getAccumulatedTimeXmX2Norm();
      puri_stats[stats_prefix + "puri_acc_time_total"]    = puriInfo.getAccumulatedTimeTotal();
      puri_stats[stats_prefix + "puri_n_steps"]           = (double)puriInfo.getNSteps();

      globalCounter++;
      if(create_m_files == 1)
	{
	  // Output matlab files with some statistics.
	  char ffname[888];
	  sprintf(ffname, "puriTime_%i.m", globalCounter);
	  std::ofstream ff(ffname);
	  puriInfo.mTimings(ff);
	  ff.close();
	  char ggname[888];
	  sprintf(ggname, "puriInfo_%i.m", globalCounter);
	  std::ofstream gg(ggname);
	  puriInfo.mInfo(gg);
	  gg.close();
	  char hhname[888];
	  sprintf(hhname, "puriMemUsage_%i.m", globalCounter);
	  std::ofstream hh(hhname);
	  puriInfo.mMemUsage(hh);
	  hh.close();
	}

      if (eigVecLUMO && eigVecHOMO) {
	generalVector lumo, homo;
	puriInfo.getHOMOandLUMOeigVecs(lumo, homo);
 	triangMatrix invCholFactor_tmp(invCholFactor);
 	invCholFactor_tmp.readFromFile();
	if ( !lumo.is_empty() ) {
	  *eigVecLUMO = lumo;
	  // perform congruence transformation 
	  (*eigVecLUMO) = invCholFactor_tmp * (*eigVecLUMO);
	}
	if ( !homo.is_empty() ) {
	  *eigVecHOMO = homo;
	  // perform congruence transformation 
	  (*eigVecHOMO) = invCholFactor_tmp * (*eigVecHOMO);
	}
      } // Note: invCholFactor_tmp goes out of scope

      intervalType homoIntervalNew = puriInfo.getHomoF();
      intervalType lumoIntervalNew = puriInfo.getLumoF();

      assert(!intervalType::intersect(homoInterval_F_ort_prev, homoIntervalNew).empty());
      assert(!intervalType::intersect(lumoInterval_F_ort_prev, lumoIntervalNew).empty());

      // Save the improved HOMO/LUMO intervals of F_ort:
      homoInterval_F_ort_prev = homoIntervalNew;
      lumoInterval_F_ort_prev = lumoIntervalNew;

      // Calculate HOMO_LUMO intervals of Finput. We need to expand
      // the F_ort intervals due to the truncation done earlier.
      homoInterval_Finput = homoInterval_F_ort_prev;
      lumoInterval_Finput = lumoInterval_F_ort_prev;
      homoInterval_Finput.increase( truncError_1 );
      lumoInterval_Finput.increase( truncError_1 );

      // Output info about gap.
      ergo_real gapMin = lumoInterval_F_ort_prev.low() - homoInterval_F_ort_prev.upp();
      ergo_real gapMax = lumoInterval_F_ort_prev.upp() - homoInterval_F_ort_prev.low();
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"E(LUMO) - E(HOMO) >= %22.11f = %22.11f eV", 
		(double)gapMin, (double)gapMin / UNIT_one_eV);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"E(LUMO) - E(HOMO) <= %22.11f = %22.11f eV", 
		(double)gapMax, (double)gapMax / UNIT_one_eV);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"HOMO interval : [ %17.12f %17.12f ]",
		(double)homoInterval_F_ort_prev.low(), (double)homoInterval_F_ort_prev.upp());
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"LUMO interval : [ %17.12f %17.12f ]",
		(double)lumoInterval_F_ort_prev.low(), (double)lumoInterval_F_ort_prev.upp());

      puri_stats[stats_prefix + "HOMO_LUMO_gap_lo_eV"]   = gapMin / UNIT_one_eV;
      puri_stats[stats_prefix + "HOMO_LUMO_gap_hi_eV"]   = gapMax / UNIT_one_eV;

      F.transfer(resultDens);

      if ( do_puri_mmul_tests == 1 ) { // Code to do some matrix-matrix multiplication tests
	// Copy F
	symmMatrix F_for_mmul_tests(F_ort_prev);
	F_for_mmul_tests.readFromFile();
	// Set up poly_choices vector
	std::vector<int> poly_choices;
	puriInfo.getPolys(poly_choices);
	// Set up thresh_values vector
	std::vector<ergo_real> thresh_values;
	puriInfo.getThreshValues(thresh_values);
	run_puri_mmul_tests(F_for_mmul_tests, 
			    allEigsInterval,
			    poly_choices,
			    thresh_values,
			    resultDens,
			    globalCounter);
      }
    }

  if ( do_comparison_to_simple_purification == 1 ) { // Code to test simple purification
    // Copy F
    symmMatrix F_simple(F_ort_prev);
    F_simple.readFromFile();
    run_comparison_to_simple_purification( F_simple, 
					   noOfOccupiedOrbs,
					   allEigsInterval,
					   resultDens,
					   globalCounter );    
  }

#if RUN_NEW_PURIFICATION
  /*******    CODE FOR NEW PURIFICATION     *******/
  if ( homoIntervalSaved.upp() < lumoIntervalSaved.low() ) {
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	      "**********************************");
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	      "Running new purification");
    for (int use_scaling = 0; use_scaling <=1; use_scaling++) {
      // Run new purification
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"************** RUN WITH use_scaling = %i.", use_scaling);    
      // Copy F
      symmMatrix F_and_D(F_ort_prev);
      F_and_D.readFromFile();
      Util::TimeMeter timeMeterNewPuri;
      pur::Purification_scaled<symmMatrix> myPuri( F_and_D,  
						   allEigsInterval,
						   homoIntervalSaved, 
						   lumoIntervalSaved,
						   eigvalueThr_Puri,
						   subspaceThr_Puri,
						   maxMul,
						   truncationNormPurification,
						   use_scaling );
      myPuri.purify();
      timeMeterNewPuri.print(LOG_AREA_DENSFROMF, "Puri: constr. + purify()");
      ergo_real densDiff = symmMatrix::eucl_diff(F_and_D, resultDens, 1e-10);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
		"Density difference: %22.11f.", (double)densDiff);    
      {
	std::stringstream ss;
	ss << "Original HOMO interval: " << homoIntervalSaved;
	do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str());
      }
      {
	std::stringstream ss;
	ss << "Original LUMO interval: " << lumoIntervalSaved;
	do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str());
      }
      {
	mat::Interval<ergo_real> hoF_new;
	mat::Interval<ergo_real> luF_new;
	myPuri.get_homo_lumo_intervals(hoF_new, luF_new);
	std::stringstream ss;
	ss << "Computed HOMO interval: " << hoF_new;
	do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str());
	std::stringstream ss2;
	ss2 << "Computed LUMO interval: " << luF_new;
	do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss2.str().c_str());
      }
      {
	std::stringstream ss;
	ss << "puriNewInfo_" << globalCounter << "_scal_" << use_scaling << ".m";
	std::ofstream ff( ss.str().c_str() );
	myPuri.mInfo(ff);
	ff.close();
      }
      {
	std::stringstream ss;
	ss << "puriNewTime_" << globalCounter << "_scal_" << use_scaling << ".m";
	std::ofstream ff( ss.str().c_str() );
	myPuri.mTime(ff);
	ff.close();
      }
    }
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	      "**********************************");
  } // end running new purification
  /******* END OF CODE FOR NEW PURIFICATION *******/
#endif
  

  if ( do_sparsity_investigation == 1 ) {
    output_magnitude_histogram( resultDens, stats_prefix + "mag_histogram_D_ort", sparsity_plots_resolution_m);
  }
  
  // Check trace of resulting density matrix, and scale it to force correct trace.
  ergo_real trace = resultDens.trace();
  ergo_real wantedTrace = noOfOccupiedOrbs;
  ergo_real traceError = trace - wantedTrace;
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	    "Trace of resulting density matrix is %22.11f, error is %18.14f.", 
	    (double)trace, (double)traceError);

  // Do truncation to speed up following multiplication operation.
  ergo_real threshold_2 = subspaceThr_2 * (1-2*eigvalueThr_Puri) / (1+subspaceThr_2);
  // Make sure that eigenvalue movement is not too large:
  threshold_2 = eigvalueThr_2 < threshold_2 ? eigvalueThr_2 : threshold_2;
  double nnzD_before_trunc_pc  = (double)resultDens.nnz()  * 100 / ((double)n*n); 
  ergo_real truncError_2 = resultDens.thresh( threshold_2, truncationNormPurification );
  double nnzD_after_trunc_pc  = (double)resultDens.nnz()  * 100 / ((double)n*n); 
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	    "Truncated D_ort ( %s ), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
	    mat::getNormTypeString(truncationNormPurification).c_str(), (double)threshold_2, (double)truncError_2, nnzD_before_trunc_pc, nnzD_after_trunc_pc);
  puri_stats[stats_prefix + "nnz_percentage_D_ort"]   = nnzD_after_trunc_pc;

  {
    triangMatrix invCholFactor_tmp(invCholFactor);
    invCholFactor_tmp.readFromFile();
    output_current_memory_usage(LOG_AREA_DENSFROMF, "Before D_S = Z * D_ort * ZT");
    Util::TimeMeter timeMeterWriteAndReadAll;
    std::string sizesStr = mat::FileWritable::writeAndReadAll();
    timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());
    Util::TimeMeter timeMeterDortTransf;
    resultDens = invCholFactor_tmp * resultDens * transpose(invCholFactor_tmp);
    timeMeterDortTransf.print(LOG_AREA_DENSFROMF, "resultDens = invCholFactor * resultDens * transpose(invCholFactor)");
    output_current_memory_usage(LOG_AREA_DENSFROMF, "After D_S = Z * D_ort * ZT");

    // Do truncation again, to reduce memory usage.
    ergo_real threshold_3 = subspaceThr_3 * (1-2*eigvalueThr_Puri-2*truncError_2) / (1+subspaceThr_3); 
    // Make sure that eigenvalue movement is not too large:
    threshold_3 = eigvalueThr_3 < threshold_3 ? eigvalueThr_3 : threshold_3;
    // Do truncation, taking into account that we are in 'non-orthogonal basis', passing invCholFactor to thresh 
    double nnzD_S_before_trunc_pc  = (double)resultDens.nnz()  * 100 / ((double)n*n); 
    ergo_real truncError_3 = resultDens.eucl_thresh( threshold_3, &invCholFactor_tmp );
    double nnzD_S_after_trunc_pc  = (double)resultDens.nnz()  * 100 / ((double)n*n); 
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	      "Truncated D_S (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
	      (double)threshold_3, (double)truncError_3, nnzD_S_before_trunc_pc, nnzD_S_after_trunc_pc);
    puri_stats[stats_prefix + "nnz_percentage_D_S"]   = nnzD_S_after_trunc_pc;    
  } // invCholFactor_tmp goes out of scope here
  
  resultDens *= factor;
  
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse ending OK");
  timeMeterTot.print(LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse");
  
  return 0;
}
