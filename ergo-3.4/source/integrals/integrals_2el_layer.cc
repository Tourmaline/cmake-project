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

#include <stdlib.h>
#include <math.h>

#include "integrals_2el_layer.h"
#include "integrals_2el.h"
#include "integrals_2el_boxed.h"
#include "integrals_2el_coulomb.h"
#include "integrals_2el_exchange.h"
#include "utilities.h"
#include "output.h"
#include "memorymanag.h"
#include "densityfitting.h"


static ergo_real* difdenSavedDensityMatrix    = NULL;
static ergo_real* difdenSavedResultFockMatrix = NULL;
static int difdensCount = 0;

int 
compute_2e_matrix_list_difden(const BasisInfoStruct & basisInfo,
			      const IntegralInfo & integralInfo,
			      const JK::ExchWeights & CAM_params,
			      ergo_real** resultList,
			      ergo_real** densList,
			      int noOfMatrices,
			      const JK::Params& J_K_params)
{
  if(noOfMatrices != 1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_2e_matrix_list_difden: (noOfMatrices != 1) not implemented!");
      return -1;
    }
  ergo_real* dens = densList[0];
  int n = basisInfo.noOfBasisFuncs;

  if(difdensCount == 10)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "restarting differential density procedure");
      ergo_free(difdenSavedDensityMatrix);
      difdenSavedDensityMatrix = NULL;
      ergo_free(difdenSavedResultFockMatrix);
      difdenSavedResultFockMatrix = NULL;
      difdensCount = 0;
    }

  if(difdenSavedDensityMatrix == NULL)
    {
      // first time
      if(compute_2e_matrix_list(basisInfo, integralInfo, CAM_params, resultList, densList,
				noOfMatrices, J_K_params) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_2e_matrix_list in compute_2e_matrix_list_difden");
	  return -1;
	}
      difdenSavedDensityMatrix    = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
      difdenSavedResultFockMatrix = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
      memcpy(difdenSavedDensityMatrix, densList[0], n*n*sizeof(ergo_real));
      memcpy(difdenSavedResultFockMatrix, resultList[0], n*n*sizeof(ergo_real));
      return 0;
    }
  
  // compute difference between current density matrix and saved density matrix
  ergo_real* densDiff = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  int i;
  for(i = 0; i < n*n; i++)
    densDiff[i] = dens[i] - difdenSavedDensityMatrix[i];

  // get max abs diff
  ergo_real maxabsdiff = 0;
  for(i = 0; i < n*n; i++)
    {
      ergo_real currAbs = std::fabs(densDiff[i]);
      if(currAbs > maxabsdiff)
	maxabsdiff = currAbs;
    }

  if(maxabsdiff > 0.1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "maxabsdiff = %22.11f, too large, not using differential density this time", (double)maxabsdiff);
      if(compute_2e_matrix_list(basisInfo, integralInfo, CAM_params, resultList, densList,
				noOfMatrices, J_K_params) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_2e_matrix_list in compute_2e_matrix_list_difden");
	  return -1;
	}
      memcpy(difdenSavedDensityMatrix, densList[0], n*n*sizeof(ergo_real));
      memcpy(difdenSavedResultFockMatrix, resultList[0], n*n*sizeof(ergo_real));
      ergo_free(densDiff);
      return 0;
    }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "using differential density, maxabsdiff = %22.11f", (double)maxabsdiff);

  // rescale densDiff
  ergo_real scaleFactor = maxabsdiff;
  for(i = 0; i < n*n; i++)
    densDiff[i] /= scaleFactor;

  // Create temporary J_K_params with modified thresholds.
  JK::Params J_K_params_mod = J_K_params;
  J_K_params_mod.threshold_J = J_K_params.threshold_J / maxabsdiff;
  J_K_params_mod.threshold_K = J_K_params.threshold_K / maxabsdiff;
  if(J_K_params_mod.threshold_J > 1e-5)
    J_K_params_mod.threshold_J = 1e-5;
  if(J_K_params_mod.threshold_K > 1e-5)
    J_K_params_mod.threshold_K = 1e-5;
  
  // Compute Fock matrix from densDiff
  if(compute_2e_matrix_list(basisInfo, integralInfo, CAM_params, resultList, &densDiff,
			    noOfMatrices, J_K_params_mod) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_2e_matrix_list in compute_2e_matrix_list_difden");
      return -1;
    }
  
  // add saved result
  for(i = 0; i < n*n; i++)
    resultList[0][i] = resultList[0][i] * scaleFactor + difdenSavedResultFockMatrix[i];

  memcpy(difdenSavedDensityMatrix, densList[0], n*n*sizeof(ergo_real));
  memcpy(difdenSavedResultFockMatrix, resultList[0], n*n*sizeof(ergo_real));

  difdensCount++;

  ergo_free(densDiff);

  return 0;
}


int 
compute_2e_matrix_list(const BasisInfoStruct & basisInfo,
		       const IntegralInfo & integralInfo,
		       const JK::ExchWeights & CAM_params,
		       ergo_real** resultList,
		       ergo_real** densList,
		       int noOfMatrices,
		       const JK::Params& J_K_params)
{
  const int maxNoOfMatrices = 100;
  int i, j;
  ergo_real* Jlist[maxNoOfMatrices];
  int nbast = basisInfo.noOfBasisFuncs;
  if(noOfMatrices > maxNoOfMatrices)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_2e_matrix_list: (noOfMatrices > maxNoOfMatrices)\n");
      return -1;
    }
  for(i = 0; i < noOfMatrices; i++)
    Jlist[i] = (ergo_real*)ergo_malloc(nbast*nbast*sizeof(ergo_real));

  // Use smallest of threshold_J and threshold_K in case of common J/K computation.
  ergo_real threshold_JK = J_K_params.threshold_J;
  if(J_K_params.threshold_K < threshold_JK)
    threshold_JK = J_K_params.threshold_K;

  if(noOfMatrices > 1)
    {
      if(CAM_params.computeRangeSeparatedExchange)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		    "Error : compute_JK_single_box cannot be used when CAM params are used."
		    "To avoid this problem, set FMM flag.");
	  return -1;
	}
      if(compute_JK_single_box(basisInfo,
			       integralInfo,
			       Jlist[0],
			       resultList[0],
			       densList[0],
			       threshold_JK) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_JK_single_box");
	  return -1;
	}
      if(compute_JK_single_box(basisInfo,
			       integralInfo,
			       Jlist[1],
			       resultList[1],
			       densList[1],
			       threshold_JK) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_JK_single_box");
	  return -1;
	}
    }
  else
    {
      if(J_K_params.use_fmm == 1)
	{
	  // FMM for J
	  if(compute_J_by_boxes(basisInfo,
				integralInfo,
				J_K_params,
				Jlist[0],
				densList[0]) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes");
	      return -1;
	    }
	  // exchange
	  // first set result buffer to zero
	  memset(resultList[0], 0, nbast*nbast*sizeof(ergo_real));
	  if(compute_K_by_boxes(basisInfo,
				integralInfo,
				CAM_params,
				J_K_params,
				resultList[0],
				NULL,
				densList[0],
				NULL,
				1) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_K_by_boxes");
	      return -1;
	    }
	}
      else
	{
	  // "single box"
	  if(CAM_params.computeRangeSeparatedExchange)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
			"error : compute_JK_single_box cannot be used when CAM params are used."
			"To avoid this problem, set FMM flag.");
	      return -1;
	    }
	  if(compute_JK_single_box(basisInfo,
				   integralInfo,
				   Jlist[0],
				   resultList[0],
				   densList[0],
				   threshold_JK) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_JK_single_box");
	      return -1;
	    }
	}
      // Now the matrix K is stored in resultList[0].
      // Multiply it by exch_weight.
      int n = basisInfo.noOfBasisFuncs;
      for(i = 0; i < n*n; i++)
	resultList[0][i] *= CAM_params.alpha;
    }
  

  if(noOfMatrices == 1)
    {
      // output Tr(DJ)
      ergo_real sum = 0;
      int i, k;
      for(i = 0; i < nbast; i++)
	for(k = 0; k < nbast; k++)
	  sum += densList[0][i*nbast+k] * Jlist[0][k*nbast+i];
      ergo_real TraceOfDJ = sum;
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Tr(DJ) = %22.11f", (double)TraceOfDJ);
    }

  // Add J and K to form final result.
  for(j = 0; j < noOfMatrices; j++)
    for(i = 0; i < nbast*nbast; i++)
      resultList[j][i] += Jlist[j][i];

  for(j = 0; j < noOfMatrices; j++)
    ergo_free(Jlist[j]);
  
  return 0;
}







int 
compute_2e_matrix_exchange(const BasisInfoStruct & basisInfo,
			   const IntegralInfo & integralInfo,
			   const JK::ExchWeights & CAM_params,
			   ergo_real* K,
			   ergo_real* dens,
			   ergo_real threshold)
{
  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error : compute_2e_matrix_exchange not implemented.\n");
  return -1;
}



int compute_2e_matrices_exchange(const BasisInfoStruct & basisInfo,
				 const IntegralInfo & integralInfo,
				 const JK::ExchWeights & CAM_params,
				 int noOfMatrices,
				 ergo_real** K_list,
				 ergo_real** D_list,
				 ergo_real threshold)
{
  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error : compute_2e_matrices_exchange not implemented.\n");
  return -1;
}



int 
compute_2e_matrix_coulomb(const BasisInfoStruct & basisInfo,
			  const BasisInfoStruct & basisInfoDensFit,
			  const IntegralInfo & integralInfo,
			  ergo_real* J,
			  ergo_real* dens,
			  const JK::Params& J_K_params,
			  DensfitData* df_data)
{
  if(J_K_params.use_densfit_for_J == 1)
    {
      ergo_real* gamma = ergo_new(basisInfoDensFit.noOfBasisFuncs,ergo_real);
      if(densfit_compute_gamma(&integralInfo,
			       basisInfo,
			       basisInfoDensFit,
			       dens,
			       gamma,
			       J_K_params.threshold_J) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in densfit_compute_gamma");
	  return -1;
	}
      ergo_real* c_vector = ergo_new(basisInfoDensFit.noOfBasisFuncs,
                                     ergo_real);
      if(densfit_compute_c_vector(&integralInfo,
				  basisInfoDensFit,
				  df_data,
				  gamma,
				  c_vector) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in densfit_compute_c_vector");
	  return -1;
	}
      if(densfit_compute_J(&integralInfo,
			   basisInfo,
			   basisInfoDensFit,
			   c_vector,
			   J,
			   J_K_params.threshold_J) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in densfit_compute_J");
	  return -1;
	}
      ergo_free(gamma);
      ergo_free(c_vector);
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "J matrix construction using density fitting complete.");
    }
  else
    {
      // FMM
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "calling compute_J_by_boxes");
      if(compute_J_by_boxes(basisInfo,
			    integralInfo,
			    J_K_params,
			    J,
			    dens) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes");
	  return -1;
	}
    }
  return 0;
}

