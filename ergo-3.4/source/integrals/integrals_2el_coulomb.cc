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


#include <cstring>
#include <cstdio>
#include <cassert>

#include <pthread.h>

#include "integrals_2el_coulomb.h"
#include "integrals_2el_utils.h"
#include "mm_limit_table.h"
#include "basis_func_pair_list.h"
#include "integrals_2el_repeating.h"
#include "integrals_hermite.h"
#include "integrals_general.h"
#include "utilities.h"
#include "pi.h"
#include "integrals_2el_util_funcs.h"


static const int HUGE_INTEGER_NUMBER = 2000000000;


typedef struct
{
  int boxIndex_1;
  int boxIndex_2;
  int branchIndex_1;
  int branchIndex_2;
} job_list_standard_entry_J_struct;

typedef struct
{
  int boxIndex;
  int multipoleBoxIndex;
  short int branchIndex;
  short int multipoleBranchIndex;
} job_list_multipole_entry_J_struct;


static int 
get_J_contribs_from_2_interacting_boxes(const BasisInfoStruct & basisInfo,
					const IntegralInfo & integralInfo,
					int maxNoOfMonomials,
					ergo_real* result_J_list,
					const distr_org_struct & distr_org_struct_1,
					const distr_org_struct & distr_org_struct_2,
					int interactionWithSelf,
					ergo_real threshold,
					JK_contribs_buffer_struct* bufferStructPtr)
{
  const JK::ExchWeights CAM_params_not_used;

  const ergo_real twoTimesPiToPow5half = 2 * pitopow52;
  ergo_real* summedIntegralList = bufferStructPtr->summedIntegralList;
  ergo_real* primitiveIntegralList = bufferStructPtr->primitiveIntegralList;

  const distr_group_struct* groupList_1 = &distr_org_struct_1.groupList[0];
  const distr_group_struct* groupList_2 = &distr_org_struct_2.groupList[0];
  const cluster_struct* clusterList_1 = &distr_org_struct_1.clusterList[0];
  const cluster_struct* clusterList_2 = &distr_org_struct_2.clusterList[0];
  const chunk_struct* chunkList_1 = &distr_org_struct_1.chunkList[0];
  int nChunks_1 = distr_org_struct_1.chunkCount;
  const chunk_struct* chunkList_2 = &distr_org_struct_2.chunkList[0];
  int nChunks_2 = distr_org_struct_2.chunkCount;
  const basis_func_pair_struct* basisFuncPairList_1 = &distr_org_struct_1.basisFuncPairList[0];
  const basis_func_pair_struct* basisFuncPairList_2 = &distr_org_struct_2.basisFuncPairList[0];
  const i_j_val_struct* spMatElementList_1 = &distr_org_struct_1.spMatElementList[0];
  const int* spMatCountList_1 = &distr_org_struct_1.spMatCountList[0];
  const int* spMatIdxList_1 = &distr_org_struct_1.spMatIdxList[0];
  const i_j_val_struct* spMatElementList_2 = &distr_org_struct_2.spMatElementList[0];
  const int* spMatCountList_2 = &distr_org_struct_2.spMatCountList[0];
  const int* spMatIdxList_2 = &distr_org_struct_2.spMatIdxList[0];

  for(int chunk_i = 0; chunk_i < nChunks_1; chunk_i++)
    {
      int chunk_j_start = 0;
      if(interactionWithSelf == 1)
	chunk_j_start = chunk_i;
      for(int chunk_j = chunk_j_start; chunk_j < nChunks_2; chunk_j++)
	{
	  int noOfBasisFuncPairs_1 = chunkList_1[chunk_i].noOfBasisFuncPairs;
	  int noOfBasisFuncPairs_2 = chunkList_2[chunk_j].noOfBasisFuncPairs;
	  // set integral list to zero
	  memset(summedIntegralList, 0, noOfBasisFuncPairs_1*noOfBasisFuncPairs_2*sizeof(ergo_real));

	  // get largest dmat element
	  ergo_real maxabsdmatelement = 0;
	  for(int i = 0; i < noOfBasisFuncPairs_1; i++)
	    for(int j = 0; j < noOfBasisFuncPairs_2; j++)
	      {
		ergo_real D_ab = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+i].dmatElement;
		ergo_real D_cd = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+j].dmatElement;
		ergo_real absval;
		absval = std::fabs(D_ab);
		if(absval > maxabsdmatelement)
		  maxabsdmatelement = absval;
		absval = std::fabs(D_cd);
		if(absval > maxabsdmatelement)
		  maxabsdmatelement = absval;
	      } // END FOR i j get largest dmat element

	  int cluster_i_start = chunkList_1[chunk_i].clusterStartIndex;
	  int clusterCount1 = chunkList_1[chunk_i].noOfClusters;
	  for(int cluster_i = cluster_i_start; cluster_i < cluster_i_start + clusterCount1; cluster_i++)
	    {
	      int cluster_j_start = chunkList_2[chunk_j].clusterStartIndex;
	      int clusterCount2 = chunkList_2[chunk_j].noOfClusters;
	      int cluterIndexEnd2 = cluster_j_start + clusterCount2;
	      if(interactionWithSelf == 1 && chunk_i == chunk_j)
		cluster_j_start = cluster_i;
	      for(int cluster_j = cluster_j_start; cluster_j < cluterIndexEnd2; cluster_j++)
		{
		  // check if we can skip this combination of clusters
		  if(clusterList_1[cluster_i].maxLimitingFactorForCluster * clusterList_2[cluster_j].maxLimitingFactorForCluster * maxabsdmatelement < threshold)
		    continue;

		  int group_i_start = clusterList_1[cluster_i].groupStartIndex;
		  int group_i_end = group_i_start + clusterList_1[cluster_i].noOfGroups;
		  int group_j_start = clusterList_2[cluster_j].groupStartIndex;
		  int group_j_end = group_j_start + clusterList_2[cluster_j].noOfGroups;

		  int n1max = clusterList_1[cluster_i].nmax;
		  int n2max = clusterList_2[cluster_j].nmax;

		  // Now we can precompute things that depend only on exponents
		  ergo_real alpha_1 = groupList_1[group_i_start].exponent;
		  ergo_real alpha_2 = groupList_2[group_j_start].exponent;
		  ergo_real alphasum = alpha_1 + alpha_2;
		  ergo_real alphaproduct = alpha_1 * alpha_2;
		  ergo_real alpha_0 = alphaproduct / alphasum;

		  ergo_real resultPreFactor = twoTimesPiToPow5half / (alphaproduct*std::sqrt(alphasum));

		  for(int group_i = group_i_start; group_i < group_i_end; group_i++)
		    {
		      if(interactionWithSelf == 1 && chunk_i == chunk_j && cluster_i == cluster_j)
			group_j_start = group_i;
		      for(int group_j = group_j_start; group_j < group_j_end; group_j++)
			{
			  // Only J is considered; we can use maxAbsDmatElementGroup
			  ergo_real maxabs_1 = groupList_1[group_i].maxAbsDmatElementGroup;
			  ergo_real maxabs_2 = groupList_2[group_j].maxAbsDmatElementGroup;
			  if((groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabs_1 < threshold) && 
			     (groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabs_2 < threshold))
			    continue;

			  // now we can do all integrals needed for this pair of groups
			  ergo_real dx = groupList_2[group_j].centerCoords[0] - groupList_1[group_i].centerCoords[0];
			  ergo_real dy = groupList_2[group_j].centerCoords[1] - groupList_1[group_i].centerCoords[1];
			  ergo_real dz = groupList_2[group_j].centerCoords[2] - groupList_1[group_i].centerCoords[2];

			  // now we have dx dy dz alpha0 alpha1 n1max n2max. Get all integrals for this case.
			  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1max];
			  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2max];

			  get_related_integrals_hermite(integralInfo,
							CAM_params_not_used,
							n1max, noOfMonomials_1,
							n2max, noOfMonomials_2,
							dx, dy, dz, alpha_0,
							resultPreFactor,
							primitiveIntegralList);

			  if(interactionWithSelf == 1 && group_j == group_i && chunk_i == chunk_j && cluster_i == cluster_j) {
			    do_summedIntegralList_contribs_self(&spMatElementList_1[spMatIdxList_1[group_i]], spMatCountList_1[group_i],
								&spMatElementList_2[spMatIdxList_2[group_j]], spMatCountList_2[group_j],
								noOfMonomials_1, noOfMonomials_2,
								primitiveIntegralList,
								noOfBasisFuncPairs_1, noOfBasisFuncPairs_2,
								summedIntegralList);
			  }
			  else {
			    do_summedIntegralList_contribs_std(&spMatElementList_1[spMatIdxList_1[group_i]], spMatCountList_1[group_i],
							       &spMatElementList_2[spMatIdxList_2[group_j]], spMatCountList_2[group_j],
							       noOfMonomials_1, noOfMonomials_2,
							       primitiveIntegralList,
							       noOfBasisFuncPairs_1, noOfBasisFuncPairs_2,
							       summedIntegralList);
			  }

			} // END FOR group_j
		    } // END FOR group_i
		} // END FOR cluster_j
	    } // END FOR cluster_i

	  for(int idx_1 = 0; idx_1 < noOfBasisFuncPairs_1; idx_1++)
	    for(int idx_2 = 0; idx_2 < noOfBasisFuncPairs_2; idx_2++)
	      {
		int a = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_1;
		int b = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_2;
		int c = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_1;
		int d = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_2;
		ergo_real integralValueCurr = summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_2];

		ergo_real D_ab = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].dmatElement;
		ergo_real D_cd = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].dmatElement;

		int J_list_index_ab = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].pairIndex;
		int J_list_index_cd = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].pairIndex;
		
		if(J_list_index_ab == J_list_index_cd)
		  integralValueCurr *= 2;

		if(std::fabs(integralValueCurr)*maxabsdmatelement < threshold)
		  continue;

		if(a != b && c != d && (a != c || b != d)) {
		  result_J_list[J_list_index_ab] += 2 * D_cd * integralValueCurr;
		  result_J_list[J_list_index_cd] += 2 * D_ab * integralValueCurr;
		}
		else if(a != b && c != d && a == c && b == d) {
		  result_J_list[J_list_index_ab] += 2 * D_cd * integralValueCurr;
		}
		else if(a == b && c != d) {
		  result_J_list[J_list_index_ab] += 2 * D_cd * integralValueCurr;
		  result_J_list[J_list_index_cd] += 1 * D_ab * integralValueCurr;
		}
		else if(a != b && c == d) {
		  result_J_list[J_list_index_ab] += 1 * D_cd * integralValueCurr;
		  result_J_list[J_list_index_cd] += 2 * D_ab * integralValueCurr;
		}
		else if(a == b && c == d && a != c) {
		  result_J_list[J_list_index_ab] += D_cd * integralValueCurr;
		  result_J_list[J_list_index_cd] += D_ab * integralValueCurr;
		}
		else if(a == b && c == d && a == c) {
		  result_J_list[J_list_index_ab] += D_cd * integralValueCurr;
		}
		else {
		  return -1; // This should never happen
		}

	      } // END FOR idx_1 idx_2
	} // END FOR chunk_j
    } // END FOR chunk_i

  return 0;
}




static int
do_multipole_interaction_between_2_boxes_branches(const distr_list_description_struct* distrDescription_1,
						  const multipole_struct_large* branchMultipole,
						  const multipole_struct_small* multipoleList_1,
						  ergo_real* result_J_list,
						  ergo_real threshold,
						  int* largest_L_used_so_far,
						  MMInteractor & interactor
						  )
{
  const chunk_struct* chunkList_1                 = &distrDescription_1->org.chunkList[0];
  const cluster_struct* clusterList_1             = &distrDescription_1->org.clusterList[0];
  const distr_group_struct* groupList_1           = &distrDescription_1->org.groupList[0];
  const minimal_distr_struct* minimalDistrList_1  = &distrDescription_1->org.minimalDistrList[0];
  int chunkCount_1                                = distrDescription_1->org.chunkCount;
  const basis_func_pair_struct* basisFuncPairList = &distrDescription_1->org.basisFuncPairList[0];

  int distrCountTot = 0;

  for(int chunkIndex_1 = 0; chunkIndex_1 < chunkCount_1; chunkIndex_1++)
    {
      int clusterCount_1 = chunkList_1[chunkIndex_1].noOfClusters;
      int cluster_start_1 = chunkList_1[chunkIndex_1].clusterStartIndex;

      for(int clusterIndex_1 = cluster_start_1; clusterIndex_1 < cluster_start_1 + clusterCount_1; clusterIndex_1++)
	{
	  int group_start_1 = clusterList_1[clusterIndex_1].groupStartIndex;
	  int group_end_1 = group_start_1 + clusterList_1[clusterIndex_1].noOfGroups;
	  for(int groupIndex_1 = group_start_1; groupIndex_1 < group_end_1; groupIndex_1++)
	    {
	      const distr_group_struct* currGroup_1 = &groupList_1[groupIndex_1];

	      ergo_real dx = branchMultipole->centerCoords[0] - currGroup_1->centerCoords[0];
	      ergo_real dy = branchMultipole->centerCoords[1] - currGroup_1->centerCoords[1];
	      ergo_real dz = branchMultipole->centerCoords[2] - currGroup_1->centerCoords[2];
	      ergo_real r = std::sqrt(dx*dx + dy*dy + dz*dz);

	      // loop over distrs of 1 (and at the same time over multipoles for those distrs)
	      // in order to find largest norm for each subvector.
	      int distr_start = currGroup_1->startIndex;
	      int distr_end = distr_start + currGroup_1->distrCount;
	      ergo_real maxMomentVectorNormForDistrsList[MAX_MULTIPOLE_DEGREE_BASIC+1];
	      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		maxMomentVectorNormForDistrsList[l] = 0;
	      int maxDegreeForDistrs = 0;
	      for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++)
		{
		  const multipole_struct_small* distrMultipole = &multipoleList_1[distrCountTot + distrIndex - distr_start];
		  if(distrMultipole->degree > maxDegreeForDistrs)
		    maxDegreeForDistrs = distrMultipole->degree;
		  for(int l = 0; l <= distrMultipole->degree; l++)
		    {
		      int startIndex = l*l;
		      int endIndex = (l+1)*(l+1);
		      ergo_real sum = 0;
		      for(int A = startIndex; A < endIndex; A++)
			sum += distrMultipole->momentList[A]*distrMultipole->momentList[A];
		      ergo_real subNorm = std::sqrt(sum);
		      if(subNorm > maxMomentVectorNormForDistrsList[l])
			maxMomentVectorNormForDistrsList[l] = subNorm;
		    }
		}

	      // check which degree is needed
	      int degreeNeeded = mm_limits_get_minimum_multipole_degree_needed(r, 
									       branchMultipole, 
									       maxDegreeForDistrs, 
									       maxMomentVectorNormForDistrsList, 
									       threshold);
	      if(degreeNeeded < 0)
		return -1;

	      if(largest_L_used_so_far != NULL)
		{
		  if(degreeNeeded > *largest_L_used_so_far)
		    *largest_L_used_so_far = degreeNeeded;
		}

	      int branchNoOfMoments = (degreeNeeded+1)*(degreeNeeded+1);

	      // create interaction matrix
	      ergo_real T[currGroup_1->multipolePtr->noOfMoments * branchNoOfMoments];
	      interactor.getInteractionMatrix(dx, dy, dz, currGroup_1->multipolePtr->degree, degreeNeeded, T);

	      ergo_real tempVector[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	      for(int A = 0; A < currGroup_1->multipolePtr->noOfMoments; A++)
		{
		  ergo_real sum = 0;
		  for(int B = 0; B < branchNoOfMoments; B++)
		    sum += branchMultipole->momentList[B] * T[A*branchNoOfMoments+B];
		  tempVector[A] = sum;
		}
	      
	      // loop over distrs of 1 (and at the same time over multipoles for those distrs)
	      for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++)
		{
		  const multipole_struct_small* distrMultipole = &multipoleList_1[distrCountTot];
		  distrCountTot++;
		  ergo_real sum = 0;
		  for(int A = 0; A < distrMultipole->noOfMoments; A++)
		    sum += tempVector[A] * distrMultipole->momentList[A];
		  int basisFuncPairIndex = minimalDistrList_1[distrIndex].basisFuncPairIndex;
		  int pairIndex = basisFuncPairList[chunkList_1[chunkIndex_1].basisFuncPairListIndex+basisFuncPairIndex].pairIndex;
		  result_J_list[pairIndex] += sum;
		} // END FOR distrIndex
	    }
	}
    }

    return 0;
}





static int
add_multipole_jobs_for_2_boxes_branches_recursive(int multipoleBoxIndex,
						  int multipoleBranchIndex,
						  int n, 
						  const box_struct* boxList,
						  int boxIndex,
						  int branchIndex,
						  int numberOfLevels,
						  int currLevel,
						  job_list_multipole_entry_J_struct* jobList_multipole,
						  int maxNoOfJobs_multipole
						  )
{
  if(currLevel == numberOfLevels - 1)
    {
      if(maxNoOfJobs_multipole <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in add_multipole_jobs_for_2_boxes_branches_recursive: (maxNoOfJobs_multipole <= 0)");
	  return -1;
	}
      if(jobList_multipole != NULL)
	{
	  jobList_multipole[0].boxIndex = boxIndex;
	  jobList_multipole[0].branchIndex = branchIndex;
	  jobList_multipole[0].multipoleBoxIndex = multipoleBoxIndex;
	  jobList_multipole[0].multipoleBranchIndex = multipoleBranchIndex;
	}
      return 1;
    }
  // go through children
  int noOfChildren = boxList[boxIndex].basicBox.noOfChildBoxes;
  int noOfNewJobs = 0;
  for(int i = 0; i < noOfChildren; i++)
    {
      int childIndex = boxList[boxIndex].basicBox.firstChildBoxIndex + i;
      job_list_multipole_entry_J_struct* jobListPtr = NULL;      
      if(jobList_multipole != NULL)
	jobListPtr = &jobList_multipole[noOfNewJobs];
      int nJobs = add_multipole_jobs_for_2_boxes_branches_recursive(multipoleBoxIndex,
								    multipoleBranchIndex,
								    n,
								    boxList,
								    childIndex,
								    branchIndex,
								    numberOfLevels,
								    currLevel + 1,
								    jobListPtr,
								    maxNoOfJobs_multipole - noOfNewJobs
								    );
	if(nJobs < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in add_multipole_jobs_for_2_boxes_branches_recursive");
	  return -1;
	}
      noOfNewJobs += nJobs;
    }  
  return noOfNewJobs;
}






static ergo_real 
get_min_distance_from_point_to_box(const ergo_real* boxCenterCoords,
				   ergo_real halfwidth,
				   const ergo_real* point)
{
  ergo_real dxList[3];
  for(int k = 0; k < 3; k++)
    {
      ergo_real dx = std::fabs(boxCenterCoords[k] - point[k]);
      if(dx > halfwidth)
	dxList[k] = dx - halfwidth;
      else
	dxList[k] = 0;
    }
  ergo_real sum = 0;
  for(int k = 0; k < 3; k++)
    sum += dxList[k] * dxList[k];
  return std::sqrt(sum);
}






static int
get_joblists_J_for_two_boxes_recursive(const BasisInfoStruct & basisInfo,
				       const IntegralInfo & integralInfo,
				       int maxNoOfMonomials,
				       ergo_real threshold,
				       const box_struct* boxList,
				       int numberOfLevels,
				       int currLevel,
				       int boxIndex_1,
				       int boxIndex_2,
				       int branchIndex_1,
				       int branchIndex_2,

				       job_list_standard_entry_J_struct* jobList_standard,
				       int maxNoOfJobs_standard,
				       int* noOfNewJobs_standard,

				       job_list_multipole_entry_J_struct* jobList_multipole,
				       int maxNoOfJobs_multipole,
				       int* noOfNewJobs_multipole
				       )
{
  // check if multipoles can be used.
  // start by computing the minimum distance between the boxes.
  // We assume that both boxes have the same width.
  ergo_real dxList[3];
  for(int coordIndex = 0; coordIndex< 3; coordIndex++)
    {
      ergo_real x1 = boxList[boxIndex_1].basicBox.centerCoords[coordIndex];
      ergo_real x2 = boxList[boxIndex_2].basicBox.centerCoords[coordIndex];
      ergo_real dx = std::fabs(x1 - x2);
      ergo_real width = boxList[boxIndex_1].basicBox.width;
      if(dx > width)
	dxList[coordIndex] = dx - width;
      else
	dxList[coordIndex] = 0;
    }
  ergo_real sumOfSquares = 0;
  for(int coordIndex = 0; coordIndex< 3; coordIndex++)
    sumOfSquares += dxList[coordIndex] * dxList[coordIndex];
  ergo_real distance = std::sqrt(sumOfSquares);
  
  ergo_real maxDistanceOutsideBox_1 = boxList[boxIndex_1].branchList[branchIndex_1].org.maxDistanceOutsideBox;
  ergo_real maxDistanceOutsideBox_2 = boxList[boxIndex_2].branchList[branchIndex_2].org.maxDistanceOutsideBox;
  
  int n = basisInfo.noOfBasisFuncs;

  int useMultipoleDescription = 0;

  if(boxIndex_1 != boxIndex_2 && distance >= maxDistanceOutsideBox_1 + maxDistanceOutsideBox_2)
    {
      // The distance is OK.
      // We also want to check that the multipole degree needed is not too high.
      // For that we need max norms of subvectors for distrs of both branches.
      
      // First the case with distrs of 1 interacting with multipole of 2
      ergo_real r_1 = get_min_distance_from_point_to_box(boxList[boxIndex_1].basicBox.centerCoords,
							 boxList[boxIndex_1].basicBox.width / 2,
							 boxList[boxIndex_2].branchList[branchIndex_2].multipole.centerCoords);
      int degreeNeeded_1 = 
	mm_limits_get_minimum_multipole_degree_needed(r_1,
						      &boxList[boxIndex_2].branchList[branchIndex_2].multipole,
						      MAX_MULTIPOLE_DEGREE_BASIC,
						      boxList[boxIndex_1].branchList[branchIndex_1].maxMomentVectorNormForDistrsList,
						      threshold);
      if(degreeNeeded_1 < 0)
	return -1;
      
      // Now the case with distrs of 2 interacting with multipole of 1
      ergo_real r_2 = get_min_distance_from_point_to_box(boxList[boxIndex_2].basicBox.centerCoords,
							 boxList[boxIndex_2].basicBox.width / 2,
							 boxList[boxIndex_1].branchList[branchIndex_1].multipole.centerCoords);
      int degreeNeeded_2 = 
	mm_limits_get_minimum_multipole_degree_needed(r_2,
						      &boxList[boxIndex_1].branchList[branchIndex_1].multipole,
						      MAX_MULTIPOLE_DEGREE_BASIC,
						      boxList[boxIndex_2].branchList[branchIndex_2].maxMomentVectorNormForDistrsList,
						      threshold);
      if(degreeNeeded_2 < 0)
	return -1;
      
      // TODO: check if this is really safe. In the V-matrix case it
      // turned out we needed two degrees margin compared to
      // MAX_MULTIPOLE_DEGREE, because in some cases the box multipole
      // is alternating between odd/even large/small
      // elements. Probably that can happen here too, and then we
      // could get large errors from that?
      if(degreeNeeded_1 < MAX_MULTIPOLE_DEGREE && degreeNeeded_2 < MAX_MULTIPOLE_DEGREE)
	useMultipoleDescription = 1;
    }

  if(useMultipoleDescription == 1)
    {
      // Use multipole description
      int noOfNewJobs_1 = add_multipole_jobs_for_2_boxes_branches_recursive(boxIndex_2,
									    branchIndex_2,
									    n,
									    boxList,
									    boxIndex_1,
									    branchIndex_1,
									    numberOfLevels,
									    currLevel,
									    jobList_multipole,
									    maxNoOfJobs_multipole
									    );
      job_list_multipole_entry_J_struct* secondPtr = NULL;
      if(jobList_multipole != NULL)
	secondPtr = &jobList_multipole[noOfNewJobs_1];
      int noOfNewJobs_2 = add_multipole_jobs_for_2_boxes_branches_recursive(boxIndex_1,
									    branchIndex_1,
									    n,
									    boxList,
									    boxIndex_2,
									    branchIndex_2,
									    numberOfLevels,
									    currLevel,
									    secondPtr,
									    maxNoOfJobs_multipole - noOfNewJobs_1
									    );
      if(noOfNewJobs_1 < 0 || noOfNewJobs_2 < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in add_multipole_jobs_for_2_boxes_branches_recursive");
	  return -1;
	}
      *noOfNewJobs_standard = 0;
      *noOfNewJobs_multipole = noOfNewJobs_1 + noOfNewJobs_2;

      return 0;
    }

  // Multipoles could not be used. We must either go to the next level or compute integrals explicitly.
  if(currLevel == numberOfLevels-1)
    {
      // We are at the level of smallest boxes. Add standard job to job list.

      if(boxIndex_1 == boxIndex_2 && branchIndex_1 > branchIndex_2)
	{
	  *noOfNewJobs_standard = 0;
	  *noOfNewJobs_multipole = 0;
	  return 0;
	}
      if(maxNoOfJobs_standard <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive: (maxNoOfJobs_standard <= 0)");
	  return -1;
	}
      if(jobList_standard != NULL)
	{
	  jobList_standard[0].boxIndex_1 = boxIndex_1;
	  jobList_standard[0].branchIndex_1 = branchIndex_1;
	  jobList_standard[0].boxIndex_2 = boxIndex_2;
	  jobList_standard[0].branchIndex_2 = branchIndex_2;
	}
      *noOfNewJobs_standard = 1;
      *noOfNewJobs_multipole = 0;
      return 0;      
    }

  // Go to next level. Do interaction between all pairs of children of the two boxes.
  int noOfChildren_1 = boxList[boxIndex_1].basicBox.noOfChildBoxes;
  int noOfChildren_2 = boxList[boxIndex_2].basicBox.noOfChildBoxes;

  if(noOfChildren_1 <= 0 || noOfChildren_2 <= 0)
    exit(0);

  int noOfNewJobs_standard_count = 0;
  int noOfNewJobs_multipole_count = 0;
  for(int i = 0; i < noOfChildren_1; i++)
    {
      int start_j = 0;
      if(boxIndex_1 == boxIndex_2)
	start_j = i;
      for(int j = start_j; j < noOfChildren_2; j++)
	{
	  int childIndex_1 = boxList[boxIndex_1].basicBox.firstChildBoxIndex + i;
	  int childIndex_2 = boxList[boxIndex_2].basicBox.firstChildBoxIndex + j;
	  int nJobs_standard = 0;
	  int nJobs_multipole = 0;
	  job_list_multipole_entry_J_struct* jobList_multipole_mod = NULL;
	  if(jobList_multipole != NULL)
	    jobList_multipole_mod = &jobList_multipole[noOfNewJobs_multipole_count];
	  job_list_standard_entry_J_struct* jobList_standard_mod = NULL;
	  if(jobList_standard != NULL)
	    jobList_standard_mod = &jobList_standard[noOfNewJobs_standard_count];
	  if(get_joblists_J_for_two_boxes_recursive(basisInfo,
						    integralInfo,
						    maxNoOfMonomials,
						    threshold,
						    boxList,
						    numberOfLevels,
						    currLevel + 1,
						    childIndex_1,
						    childIndex_2,
						    branchIndex_1,
						    branchIndex_2,

						    jobList_standard_mod,
						    maxNoOfJobs_standard - noOfNewJobs_standard_count,
						    &nJobs_standard,

						    jobList_multipole_mod,
						    maxNoOfJobs_multipole - noOfNewJobs_multipole_count,
						    &nJobs_multipole
						    ) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive for child boxes");
	      return -1;
	    }
	  noOfNewJobs_standard_count += nJobs_standard;
	  noOfNewJobs_multipole_count += nJobs_multipole;
	} // END FOR j
    } // END FOR i
  
  *noOfNewJobs_standard = noOfNewJobs_standard_count;
  *noOfNewJobs_multipole = noOfNewJobs_multipole_count;
  return 0;
}




static int
get_list_of_labeled_distrs_maxLimitingFactor_linear(const BasisInfoStruct & basisInfo,
						    const IntegralInfo & integralInfo,
						    ergo_real threshold,
						    const basis_func_index_pair_struct* basisFuncIndexPairList,
						    int basisFuncIndexPairCount,
						    ergo_real* resultMaxLimitingFactor)
{ 
  IntegratorWithMemory integrator(&integralInfo);

  ergo_real maxLimitingFactor = 0;
  for(int kk = 0; kk < basisFuncIndexPairCount; kk++)
    {
      int i = basisFuncIndexPairList[kk].index_1;
      int j = basisFuncIndexPairList[kk].index_2;
	  
      const int maxCountProduct = POLY_PRODUCT_MAX_DISTRS;
      DistributionSpecStruct psi_list[maxCountProduct];
      /* form product of basisfuncs i and j, store product in psi_list */
      int n_psi = get_product_simple_primitives(basisInfo, i,
						basisInfo, j,
						psi_list,
						maxCountProduct,
						0);
      if(n_psi < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives");
	  return -1;
	}
      for(int k = 0; k < n_psi; k++)
	{
	  ergo_real limitingFactor = std::sqrt(integrator.do_2e_integral(&psi_list[k]));
	  if(limitingFactor > maxLimitingFactor)
	    maxLimitingFactor = limitingFactor;
	} // END FOR k
    } // END FOR kk
  *resultMaxLimitingFactor = maxLimitingFactor;

  return 0;
}




static int
get_list_of_labeled_distrs_linear(const BasisInfoStruct & basisInfo,
				  const IntegralInfo & integralInfo,
				  ergo_real threshold,
				  DistributionSpecStructLabeled* resultList,
				  int maxCountDistrs,
				  ergo_real maxLimitingFactor,
				  const basis_func_index_pair_struct* basisFuncIndexPairList,
				  int basisFuncIndexPairCount,
				  const ergo_real* D_list)
{
  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(basisFuncIndexPairCount, D_list);

  IntegratorWithMemory integrator(&integralInfo);

  // create list of product primitives, with labels
  int distrCount = 0;
  for(int kk = 0; kk < basisFuncIndexPairCount; kk++)
    {
      int i = basisFuncIndexPairList[kk].index_1;
      int j = basisFuncIndexPairList[kk].index_2;

      ergo_real dmatElement = D_list[kk];

      const int maxCountProduct = POLY_PRODUCT_MAX_DISTRS;
      DistributionSpecStruct psi_list[maxCountProduct];
      /* form product of basisfuncs i and j, store product in psi_list */
      int n_psi = get_product_simple_primitives(basisInfo, i,
						basisInfo, j,
						psi_list,
						maxCountProduct,
						0);
      if(n_psi < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives");
	  return -1;
	}
      for(int k = 0; k < n_psi; k++)
	{
	  ergo_real limitingFactor = std::sqrt(integrator.do_2e_integral(&psi_list[k]));
	  if(limitingFactor*maxLimitingFactor*maxDensityMatrixElement > threshold)
	    {
	      if(maxCountDistrs > 0 && distrCount >= maxCountDistrs)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs: (maxCountDistrs > 0 && distrCount >= maxCountDistrs)");
		  return -1;
		}
	      if(resultList != NULL)
		{
		  resultList[distrCount].distr = psi_list[k];
		  resultList[distrCount].basisFuncIndex_1 = i;
		  resultList[distrCount].basisFuncIndex_2 = j;
		  resultList[distrCount].pairIndex = kk;
		  resultList[distrCount].limitingFactor = limitingFactor;
		  resultList[distrCount].dmatElement = dmatElement;
		}
	      distrCount++;
	    } // END IF above threshold
	} // END FOR k
    } // END FOR kk

  return distrCount;
}





static int
compare_multipole_jobs(const void* p1, const void* p2)
{
  job_list_multipole_entry_J_struct* job_1 = (job_list_multipole_entry_J_struct*)p1;
  job_list_multipole_entry_J_struct* job_2 = (job_list_multipole_entry_J_struct*)p2;
  if(job_1->boxIndex > job_2->boxIndex)
    return 1;
  if(job_1->boxIndex < job_2->boxIndex)
    return -1;
  // now we know that boxIndex is the same for both
  if(job_1->branchIndex > job_2->branchIndex)
    return 1;
  if(job_1->branchIndex < job_2->branchIndex)
    return -1;
  // now we know that boxIndex and branchIndex are the same for both
  if(job_1->multipoleBoxIndex > job_2->multipoleBoxIndex)
    return 1;
  if(job_1->multipoleBoxIndex < job_2->multipoleBoxIndex)
    return -1;
  // now we know that boxIndex and branchIndex and multipoleBoxIndex are the same for both
  // we do not care about the order of different multipoleBranchIndex
  return 0;
}





static void
get_largest_and_smallest_extent_for_list_of_distributions(int n, 
							  const DistributionSpecStructLabeled* distrList, 
							  ergo_real* result_extent_min, 
							  ergo_real* result_extent_max)
{
  ergo_real extent_min = distrList[0].distr.extent;
  ergo_real extent_max = distrList[0].distr.extent;
  for(int i = 0; i < n; i++)
    {
      ergo_real extent = distrList[i].distr.extent;
      if(extent > extent_max)
	extent_max = extent;
      if(extent < extent_min)
	extent_min = extent;
    }
  *result_extent_min = extent_min;
  *result_extent_max = extent_max;
}





static int
get_max_no_of_monomials_for_list_of_distributions(int n, 
						  const DistributionSpecStructLabeled* distrList, 
						  const IntegralInfo & integralInfo)
{
  int maxNoOfMonomials = 0;
  for(int i = 0; i < n; i++)
    {
      int degree = 0;
      for(int j = 0; j < 3; j++)
	degree += distrList[i].distr.monomialInts[j];
      int noOfMonomials = integralInfo.monomial_info.no_of_monomials_list[degree];
      if(noOfMonomials > maxNoOfMonomials)
	maxNoOfMonomials = noOfMonomials;
    } // END FOR ABcount
  return maxNoOfMonomials;
}





static int
get_branch_splitter_info(ergo_real* branchSplitterList,
			 int maxNoOfBranches,
			 const JK::Params& J_K_params,
			 ergo_real toplevelBoxSize,
			 ergo_real extent_max)
{
  int noOfBranches = 0;
  if(J_K_params.fmm_no_of_branches > 0)
    {
      // Use branches as specified in input parameters
      noOfBranches = J_K_params.fmm_no_of_branches;
      if(noOfBranches >= maxNoOfBranches)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_branch_splitter_info: (noOfBranches >= maxNoOfBranches)");
	  return -1;
	}
      for(int i = 0; i < noOfBranches-1; i++)
	{
	  ergo_real splitterValue = 0;
	  switch(i)
	    {
	    case 0: splitterValue = J_K_params.fmm_branch_splitter_extent_1; break;
	    case 1: splitterValue = J_K_params.fmm_branch_splitter_extent_2; break;
	    case 2: splitterValue = J_K_params.fmm_branch_splitter_extent_3; break;
	    case 3: splitterValue = J_K_params.fmm_branch_splitter_extent_4; break;
	    case 4: splitterValue = J_K_params.fmm_branch_splitter_extent_5; break;
	    default:
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_branch_splitter_info: default reached.");
	      return -1;
	    }
	  branchSplitterList[i] = splitterValue;
	}  
    }
  else
    {
      // Use default branch settings based on box size and extent_max
      ergo_real splitterValue = toplevelBoxSize / 2;
      noOfBranches = 2;
      while(splitterValue < extent_max)
	{
	  noOfBranches++;
	  splitterValue *= 2;
	}
      // Now we know how many branches we need. Create splitter list.
      int count = 0;
      splitterValue = 0;
      branchSplitterList[noOfBranches-2-count] = splitterValue;
      count++;
      splitterValue = toplevelBoxSize / 2;
      for(count = 1; count < noOfBranches-1; count++)
	{
	  branchSplitterList[noOfBranches-2-count] = splitterValue;
	  splitterValue *= 2;
	}
    }
  return noOfBranches;
}









static int
create_branches(int noOfBranches,
		const ergo_real* branchSplitterList,
		int distrCount,
		DistributionSpecStructLabeled* distrListOrdered,
		int noOfBoxesTopLevel,
		box_struct* boxListTopLevel
		)
{
  // Start by finding out largest number of distrs per box.
  int maxNoOfDistrsPerBox = 0;
  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    int distrCountCurrBox = boxListTopLevel[i].basicBox.noOfItems;
    if(distrCountCurrBox > maxNoOfDistrsPerBox)
      maxNoOfDistrsPerBox = distrCountCurrBox;
  }
  
  std::vector<int> branchBucketIndexList[MAX_NO_OF_BRANCHES];
  int branchBucketCountList[MAX_NO_OF_BRANCHES];
  for(int i = 0; i < noOfBranches; i++)
    branchBucketIndexList[i].resize(maxNoOfDistrsPerBox);
  
  std::vector<DistributionSpecStructLabeled> distrListTemp(maxNoOfDistrsPerBox);
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating distrListTemp");

  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    DistributionSpecStructLabeled* distrListCurrBox  = &distrListOrdered[boxListTopLevel[i].basicBox.firstItemIndex];
    int distrCountCurrBox = boxListTopLevel[i].basicBox.noOfItems;
    memcpy(&distrListTemp[0], distrListCurrBox, distrCountCurrBox*sizeof(DistributionSpecStructLabeled));
    DistributionSpecStructLabeled* distrListCurrBox2 = &distrListTemp[0];
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
      branchBucketCountList[branchIndex] = 0;
    for(int j = 0; j < distrCountCurrBox; j++) {
      int branchIndex; // declare here because value is used after loop after break.
      for(branchIndex = noOfBranches-1; branchIndex > 0; branchIndex--) {
	ergo_real extent = distrListCurrBox[j].distr.extent;
	ergo_real width = boxListTopLevel[i].basicBox.width;
	// get minWallDist : minimum wall distance
	ergo_real minWallDist = width;
	for(int coordIndex = 0; coordIndex< 3; coordIndex++) {
	  // get wall distance for this coordinate
	  ergo_real dx = distrListCurrBox[j].distr.centerCoords[coordIndex] - boxListTopLevel[i].basicBox.centerCoords[coordIndex];
	  ergo_real wallDist = width - std::fabs(dx);
	  if(wallDist < minWallDist)
	    minWallDist = wallDist;
	} // END FOR coordIndex
	if(minWallDist < 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (minWallDist < 0)");
	  return -1;
	}
	if((extent - minWallDist) < branchSplitterList[branchIndex-1])
	  break;
      }
      branchBucketIndexList[branchIndex][branchBucketCountList[branchIndex]] = j;
      branchBucketCountList[branchIndex]++;
    } // END FOR j
    int newCount = 0;
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
      boxListTopLevel[i].branchIndexList[branchIndex] = boxListTopLevel[i].basicBox.firstItemIndex + newCount;	  
      boxListTopLevel[i].branchCountList[branchIndex] = branchBucketCountList[branchIndex];
      for(int k = 0; k < branchBucketCountList[branchIndex]; k++) {
	distrListCurrBox[newCount] = distrListCurrBox2[branchBucketIndexList[branchIndex][k]];
	newCount++;
      }
    } // END FOR branchIndex
  } // END FOR i divide distrs into branches according to extent.

  return 0;
}









static int
execute_joblist_J_std_serial(int noOfJobs_J_standard,
			     const job_list_standard_entry_J_struct* jobList_J_standard,
			     const BasisInfoStruct & basisInfo,
			     const IntegralInfo & integralInfo,
			     int maxNoOfMonomials,
			     ergo_real* result_J_list,
			     const box_struct* boxList,
			     ergo_real threshold)
{
  Util::TimeMeter timeMeter;

  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(integralInfo, maxNoOfMonomials, 0, &bufferStruct);

  for(int jobIndex = 0; jobIndex < noOfJobs_J_standard; jobIndex++)
    {
      int boxIndex_1 = jobList_J_standard[jobIndex].boxIndex_1;
      int boxIndex_2 = jobList_J_standard[jobIndex].boxIndex_2;
      int branchIndex_1 = jobList_J_standard[jobIndex].branchIndex_1;
      int branchIndex_2 = jobList_J_standard[jobIndex].branchIndex_2;
      int self = 0;
      if(boxIndex_1 == boxIndex_2 && branchIndex_1 == branchIndex_2)
	self = 1;
      if(get_J_contribs_from_2_interacting_boxes(basisInfo,
						 integralInfo,
						 maxNoOfMonomials,
						 result_J_list,
						 boxList[boxIndex_1].branchList[branchIndex_1].org,
						 boxList[boxIndex_2].branchList[branchIndex_2].org,
						 self,
						 threshold,
						 &bufferStruct) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_J_contribs_from_2_interacting_boxes");
	  return -1;
	}
    } // END FOR jobIndex

  free_buffers_needed_by_integral_code(&bufferStruct);

  timeMeter.print(LOG_AREA_INTEGRALS, "execute_joblist_J_std_serial");

  return 0;
}


struct J_std_joblist_thread_struct
{
  pthread_t thread;
  const BasisInfoStruct & basisInfo;
  const IntegralInfo* integralInfo;
  ergo_real* result_J_list;
  int maxNoOfMonomials;
  ergo_real threshold;
  const box_struct* boxList;
  const job_list_standard_entry_J_struct* jobList_J_standard;
  int noOfJobs_J_standard;
  int thread_ID;
  int noOfThreads;
  int resultCode;
  explicit J_std_joblist_thread_struct(const BasisInfoStruct & basisInfoIn) :
    basisInfo(basisInfoIn) { }
};


static void*
execute_joblist_J_std_thread_func(void* arg)
{
  J_std_joblist_thread_struct* params = (J_std_joblist_thread_struct*)arg;

  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(*params->integralInfo, params->maxNoOfMonomials, 0, &bufferStruct);

  const box_struct* boxList = params->boxList;

  for(int jobIndex = 0; jobIndex < params->noOfJobs_J_standard; jobIndex++)
    {
      if(jobIndex % params->noOfThreads != params->thread_ID)
	continue;

      int boxIndex_1 = params->jobList_J_standard[jobIndex].boxIndex_1;
      int boxIndex_2 = params->jobList_J_standard[jobIndex].boxIndex_2;
      int branchIndex_1 = params->jobList_J_standard[jobIndex].branchIndex_1;
      int branchIndex_2 = params->jobList_J_standard[jobIndex].branchIndex_2;
      int self = 0;
      if(boxIndex_1 == boxIndex_2 && branchIndex_1 == branchIndex_2)
	self = 1;
      if(get_J_contribs_from_2_interacting_boxes(params->basisInfo,
						 *params->integralInfo,
						 params->maxNoOfMonomials,
						 params->result_J_list,
						 boxList[boxIndex_1].branchList[branchIndex_1].org,
						 boxList[boxIndex_2].branchList[branchIndex_2].org,
						 self,
						 params->threshold,
						 &bufferStruct) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_J_contribs_from_2_interacting_boxes");
	  params->resultCode = -1;
	  return NULL;
	}
    } // END FOR jobIndex

  free_buffers_needed_by_integral_code(&bufferStruct);

  params->resultCode = 0;
  return NULL;
}


static int
execute_joblist_J_std_threaded(int noOfThreads,
			       int noOfJobs_J_standard,
			       const job_list_standard_entry_J_struct* jobList_J_standard,
			       const BasisInfoStruct & basisInfo,
			       const IntegralInfo & integralInfo,
			       int maxNoOfMonomials,
			       ergo_real* result_J_list,
			       int noOfBasisFuncIndexPairs,
			       const box_struct* boxList,
			       ergo_real threshold)
{
  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "execute_joblist_J_std_threaded, noOfThreads = %2i", noOfThreads);
  
  J_std_joblist_thread_struct* threadParamsList[noOfThreads];
  
  // Set common parameters for all threads
  for(int i = 0; i < noOfThreads; i++)
    {
      threadParamsList[i] = new J_std_joblist_thread_struct(basisInfo);
      //threadParamsList[i].basisInfo = basisInfo;
      threadParamsList[i]->integralInfo = &integralInfo;
      threadParamsList[i]->maxNoOfMonomials = maxNoOfMonomials;
      threadParamsList[i]->boxList = boxList;
      threadParamsList[i]->jobList_J_standard = jobList_J_standard;
      threadParamsList[i]->noOfJobs_J_standard = noOfJobs_J_standard;
      threadParamsList[i]->noOfThreads = noOfThreads;
      threadParamsList[i]->resultCode = -1; // initialize to error code
      threadParamsList[i]->threshold = threshold;
    } // END FOR i
  
  // Set result pointer for thread 0
  // Thread 0 uses the original result_J_list pointer.
  threadParamsList[0]->result_J_list = result_J_list;

  // Set result_J_list pointer for other threads
  for(int i = 1; i < noOfThreads; i++)
    {
      threadParamsList[i]->result_J_list = new ergo_real[noOfBasisFuncIndexPairs];
      memset(threadParamsList[i]->result_J_list, 0, noOfBasisFuncIndexPairs * sizeof(ergo_real));
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating memory for threads.");


  // Set ID number for all threads
  for(int i = 0; i < noOfThreads; i++)
    threadParamsList[i]->thread_ID = i;

  /* start threads */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(pthread_create(&threadParamsList[i]->thread, 
			NULL, 
			execute_joblist_J_std_thread_func, 
			threadParamsList[i]) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_create for thread %i", i);
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
	  for(int j = 0; j < i; j++)
	    {
	      if(pthread_join(threadParamsList[j]->thread, NULL) != 0)
		do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", j);
	    } /* END FOR j */
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "all threads finished, returning error code");
	  return -1;
	}
    } /* END FOR i */
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "%i threads started OK.", noOfThreads);

  /* wait for threads to finish */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(pthread_join(threadParamsList[i]->thread, NULL) != 0)
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
    } /* END FOR i */
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "all %i threads have finished.", noOfThreads);
  
  /* now all threads have finished, check for errors */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(threadParamsList[i]->resultCode != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_std_thread_func"
		    " for thread %i", i);
	  return -1;
	}
    } /* END FOR i */
  
  
  // add contributions from other threads
  for(int i = 1; i < noOfThreads; i++)
    {
      for(int j = 0; j < noOfBasisFuncIndexPairs; j++)
	result_J_list[j] += threadParamsList[i]->result_J_list[j];
    }
  
  // Free extra result_J_list buffers used by threads.
  // Note that this loop must start with 1, not 0.
  for(int i = 1; i < noOfThreads; i++)
    delete [] threadParamsList[i]->result_J_list;

  for(int i = 0; i < noOfThreads; i++)
    delete threadParamsList[i];

  timeMeter.print(LOG_AREA_INTEGRALS, "execute_joblist_J_std_threaded");
  
  return 0;
}


static int 
sort_list_of_multipole_jobs_fixed_boxIndex(job_list_multipole_entry_J_struct* jobList, int n)
{
  // Start by bucket-sort by branchIndex.
  const int maxNoOfBranches = 10;  
  job_list_multipole_entry_J_struct* bucketList[maxNoOfBranches];
  // Get number of branches
  int branchIndex_min = maxNoOfBranches;
  int branchIndex_max = 0;
  for(int i = 0; i < n; i++)
    {
      int currBranchIndex = jobList[i].branchIndex;
      if(currBranchIndex > branchIndex_max)
	branchIndex_max = currBranchIndex;
      if(currBranchIndex < branchIndex_min)
	branchIndex_min = currBranchIndex;
    }
  assert(branchIndex_min >= 0);
  assert(branchIndex_max < maxNoOfBranches);
  int noOfBranches = branchIndex_max + 1;
  for(int i = 0; i < noOfBranches; i++)
    bucketList[i] = new job_list_multipole_entry_J_struct[n];

  int counterList[maxNoOfBranches];
  for(int i = 0; i < maxNoOfBranches; i++)
    counterList[i] = 0;
  
  for(int i = 0; i < n; i++)
    {
      int currBranchIndex = jobList[i].branchIndex;
      assert(currBranchIndex < noOfBranches);
      int count = counterList[currBranchIndex];
      bucketList[currBranchIndex][count] = jobList[i];
      counterList[currBranchIndex]++;
    }

  // OK, bucket-sort done. Now sort contents of each bucket.
  int count = 0;
  for(int i = 0; i < maxNoOfBranches; i++)
    {
      int currCount = counterList[i];
      if(currCount == 0) continue;

      qsort(bucketList[i], currCount, sizeof(job_list_multipole_entry_J_struct),
            compare_multipole_jobs);

      // check qsort result
      for(int j  = 0; j < currCount-1; j++)
      {
        job_list_multipole_entry_J_struct* curr = &bucketList[i][j];
        job_list_multipole_entry_J_struct* next = &bucketList[i][j+1];
        if(compare_multipole_jobs(curr, next) > 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: qsort result not sorted.");
          return -1;
        }
      }

      // Copy result
      memcpy(&jobList[count], bucketList[i], currCount * sizeof(job_list_multipole_entry_J_struct));
      count += currCount;
    }

  for(int i = 0; i < noOfBranches; i++)
    delete [] bucketList[i];
  
  return 0;
}


static int 
sort_list_of_multipole_jobs(job_list_multipole_entry_J_struct* jobList, int n)
{
  if(n == 0)
    return 0;

  Util::TimeMeter timeMeterInit;

  // Start by bucket-sort by boxIndex.
  // Go through list once to find max boxIndex.
  int boxIndex_max = jobList[0].boxIndex;
  for(int i = 0; i < n; i++)
    {
      int currBoxIndex = jobList[i].boxIndex;
      if(currBoxIndex > boxIndex_max)
	boxIndex_max = currBoxIndex;
    }

  // Go through list once more to find maxNoOfJobsWithSameBoxIndex
  int noOfBoxIndexes = boxIndex_max + 1;
  std::vector<int> counterList1(noOfBoxIndexes);
  for(int i = 0; i < noOfBoxIndexes; i++)
    counterList1[i] = 0;
  for(int i = 0; i < n; i++)
    {
      int currBoxIndex = jobList[i].boxIndex;
      counterList1[currBoxIndex]++;
    }
  int maxNoOfJobsWithSameBoxIndex = 0;
  for(int i = 0; i < noOfBoxIndexes; i++)
    {
      if(counterList1[i] > maxNoOfJobsWithSameBoxIndex)
	maxNoOfJobsWithSameBoxIndex = counterList1[i];
    }

  timeMeterInit.print(LOG_AREA_INTEGRALS, "sort_list_of_multipole_jobs init part");

  //std::vector<job_list_multipole_entry_J_struct> bucketList(noOfBoxIndexes*maxNoOfJobsWithSameBoxIndex);
  std::vector< std::vector<job_list_multipole_entry_J_struct> > bucketList(noOfBoxIndexes);
  for(int i = 0; i < noOfBoxIndexes; i++)
    bucketList[i].resize(counterList1[i]);

  std::vector<int> counterList2(noOfBoxIndexes);
  for(int i = 0; i < noOfBoxIndexes; i++)
    counterList2[i] = 0;
  for(int i = 0; i < n; i++)
    {
      int currBoxIndex = jobList[i].boxIndex;
      if(counterList2[currBoxIndex] >= counterList1[currBoxIndex])
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in sort_list_of_multipole_jobs: (counterList2[currBoxIndex] >= counterList1[currBoxIndex])");
	  return -1;
	}
      bucketList[currBoxIndex][counterList2[currBoxIndex]] = jobList[i];
      counterList2[currBoxIndex]++;
    }

  // OK, bucket-sort done. Now sort contents of each bucket.
  int count = 0;
  for(int i = 0; i < noOfBoxIndexes; i++)
    {
      int currCount = counterList2[i];
      if(sort_list_of_multipole_jobs_fixed_boxIndex(&bucketList[i][0], currCount) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in sort_list_of_multipole_jobs_fixed_boxIndex");
	  return -1;
	}
      // Copy result
      memcpy(&jobList[count], &bucketList[i][0], currCount * sizeof(job_list_multipole_entry_J_struct));
      count += currCount;
    }

  // check that list is sorted
  for(int i = 0; i < n-1; i++)
    {
      if(jobList[i].boxIndex > jobList[i+1].boxIndex)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: list not sorted!");
	  return -1;
	}
      if(jobList[i].boxIndex == jobList[i+1].boxIndex)
	{
	  if(jobList[i].branchIndex > jobList[i+1].branchIndex)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: list not sorted!");
	      return -1;
	    }
	  if(jobList[i].branchIndex == jobList[i+1].branchIndex)
	    {
	      if(jobList[i].multipoleBoxIndex > jobList[i+1].multipoleBoxIndex)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: list not sorted!");
		  return -1;
		}
	    }
	}
    }

  timeMeterInit.print(LOG_AREA_INTEGRALS, "sort_list_of_multipole_jobs complete");

  return 0;
}

/** executes given jobList using FMM.
    @param jobIndexLo the first jobindex for which this thread is responsible.
    @param jobIndexHi the last jobindex for which this thread is responsible is jobIndexHi-1.
    @param integralInfo info needed for evaluation of integrals of Gaussian functions.
    @param basisInfo info about the used basis set.
    @param J_K_params includes various parameters for J and K matrix construction.
    @param jobList_J_multipole list of multipole-jobs.
    @param boxList list of boxes.
    @param maxnoOfMinimalDistrsPerBoxBranch needed to determine size of work buffer.
    @param result_J_list the list of matrix elements to be updated.
    @param largest_L_used largest L-value used (output).
*/
static int
execute_joblist_J_fmm_shared(int jobIndexLo, int jobIndexHi,
                             const IntegralInfo& integralInfo,
                             const BasisInfoStruct & basisInfo,
                             const JK::Params& J_K_params,
                             const job_list_multipole_entry_J_struct* jobList_J_multipole,
                             const box_struct *boxList,
                             int maxnoOfMinimalDistrsPerBoxBranch,
                             ergo_real* result_J_list,
                             int* largest_L_used)
{
  // Execute multipole job list for J
  int boxIndexSaved = -1;
  int branchIndexSaved = -1;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "execute_joblist_J_fmm_shared: Allocating multipoleList_4, maxnoOfMinimalDistrsPerBoxBranch = %9i", 
	    maxnoOfMinimalDistrsPerBoxBranch);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "execute_joblist_J_fmm_shared: jobIndexLo = %12d, jobIndexHi = %12d", 
	    jobIndexLo, jobIndexHi);

  std::vector<multipole_struct_small> multipoleList_4(maxnoOfMinimalDistrsPerBoxBranch);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating multipoleList_4");

  MMInteractor interactor;
  
  int jobIndex = jobIndexLo;
  while(jobIndex<jobIndexHi)
    {
      // check how many of the following jobs that differ only in multipoleBranchIndex
      int jobIndex2 = jobIndex;

      int boxIndex = jobList_J_multipole[jobIndex].boxIndex;
      int branchIndex = jobList_J_multipole[jobIndex].branchIndex;
      int multipoleBoxIndex = jobList_J_multipole[jobIndex].multipoleBoxIndex;

      while(++jobIndex2 < jobIndexHi)
        {
          if(jobList_J_multipole[jobIndex2].boxIndex != boxIndex)
            break;
          if(jobList_J_multipole[jobIndex2].branchIndex != branchIndex)
            break;
          if(jobList_J_multipole[jobIndex2].multipoleBoxIndex != multipoleBoxIndex)
            break;
        }
      int nJobs = jobIndex2 - jobIndex;

      // check if we need to create new list of multipoles
      if(boxIndex != boxIndexSaved || branchIndex != branchIndexSaved)
        {
          // create list of multipoles
          const chunk_struct* chunkList = &boxList[boxIndex].branchList[branchIndex].org.chunkList[0];
          const cluster_struct* clusterList = &boxList[boxIndex].branchList[branchIndex].org.clusterList[0];
          const distr_group_struct* groupList = &boxList[boxIndex].branchList[branchIndex].org.groupList[0];
          const minimal_distr_struct* minimalDistrList = &boxList[boxIndex].branchList[branchIndex].org.minimalDistrList[0];
          int chunkCount = boxList[boxIndex].branchList[branchIndex].org.chunkCount;
          int count_temp = 0;
          for(int chunkIndex = 0; chunkIndex < chunkCount; chunkIndex++)
            {
              int clusterCount = chunkList[chunkIndex].noOfClusters;
              int cluster_start = chunkList[chunkIndex].clusterStartIndex;
              for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++)
                {
                  int group_start = clusterList[clusterIndex].groupStartIndex;
                  int group_end = group_start + clusterList[clusterIndex].noOfGroups;
                  for(int groupIndex = group_start; groupIndex < group_end; groupIndex++)
                    {
                      const distr_group_struct* currGroup = &groupList[groupIndex];
		      
                      int distr_start = currGroup->startIndex;
                      int distr_end = distr_start + currGroup->distrCount;
                      for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++)
                        {
                          int monomialIndex = minimalDistrList[distrIndex].monomialIndex;
                          ergo_real coeff = minimalDistrList[distrIndex].coeff;
                          // get monomialInts from monomialIndex
                          DistributionSpecStruct distr;
                          distr.monomialInts[0] = integralInfo.monomial_info.monomial_list[monomialIndex].ix;
                          distr.monomialInts[1] = integralInfo.monomial_info.monomial_list[monomialIndex].iy;
                          distr.monomialInts[2] = integralInfo.monomial_info.monomial_list[monomialIndex].iz;
                          distr.coeff = coeff;
                          distr.exponent = currGroup->exponent;
                          distr.centerCoords[0] = currGroup->centerCoords[0];
                          distr.centerCoords[1] = currGroup->centerCoords[1];
                          distr.centerCoords[2] = currGroup->centerCoords[2];

                          multipole_struct_small multipole;
                          if(compute_multipole_moments(integralInfo, &distr, &multipole) != 0)
                            {
                              do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_multipole_moments");
                              return -1;
                            }

                          multipoleList_4[count_temp] = multipole;
                          count_temp++;
                        } // END FOR distrIndex

                    } // END FOR groupIndex
                } // END FOR clusterIndex
            } // END FOR chunkIndex
	  // save these boxIndex and branchIndex values, so that we do not need to recompute multipoles until we reach the next box/branch
          boxIndexSaved = boxIndex;
          branchIndexSaved = branchIndex;
        } // END IF need to create new list of multipoles


      // OK, now we have nJobs, which is at least 1
      int first_multipoleBranchIndex = jobList_J_multipole[jobIndex].multipoleBranchIndex;
      multipole_struct_large multipoleSum = boxList[multipoleBoxIndex].branchList[first_multipoleBranchIndex].multipole;
      memset(multipoleSum.momentList, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE*sizeof(ergo_real));
      for(int jobNo = 0; jobNo < nJobs; jobNo++)
	{
	  int multipoleBranchIndex = jobList_J_multipole[jobIndex+jobNo].multipoleBranchIndex;
	  const multipole_struct_large* branchMultipole = &boxList[multipoleBoxIndex].branchList[multipoleBranchIndex].multipole;
	  for(int mm = 0; mm < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; mm++)
	    multipoleSum.momentList[mm] += branchMultipole->momentList[mm];
	} // END FOR jobNo
      setup_multipole_maxAbsMomentList(&multipoleSum);

      do_multipole_interaction_between_2_boxes_branches(&boxList[boxIndex].branchList[branchIndex],
							&multipoleSum,
							&multipoleList_4[0],
							result_J_list,
							J_K_params.threshold_J * J_K_params.multipole_threshold_factor,
							largest_L_used,
							interactor
							);

      jobIndex = jobIndex2;
      
    } // END WHILE (jobIndex < jobIndexHi)

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "execute_joblist_J_fmm_shared: done!");
  
  return 0;
}

static int
execute_joblist_J_fmm_serial(const IntegralInfo& integralInfo,
                             const BasisInfoStruct & basisInfo,
                             const JK::Params& J_K_params,
                             int noOfJobs_J_multipole,
                             const job_list_multipole_entry_J_struct* jobList_J_multipole,
                             const box_struct *boxList,
                             int maxnoOfMinimalDistrsPerBoxBranch,
                             ergo_real* result_J_list)
{
  Util::TimeMeter timeMeterJmul;
  int largest_L_used = 0;
  int rc = execute_joblist_J_fmm_shared(0, noOfJobs_J_multipole,
                                        integralInfo,
                                        basisInfo,
                                        J_K_params,
                                        jobList_J_multipole,
                                        boxList,
                                        maxnoOfMinimalDistrsPerBoxBranch,
                                        result_J_list,
                                        &largest_L_used);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "multipole job list for J executed, largest L used: %2i", 
	    largest_L_used);
  timeMeterJmul.print(LOG_AREA_INTEGRALS, "Executing multipole job list for J");
  return rc;
}

struct JFMMWorkerData {
  const IntegralInfo* integralInfo;
  const BasisInfoStruct* basisInfo;
  const JK::Params* J_K_params;
  const job_list_multipole_entry_J_struct* jobList_J_multipole;
  const box_struct *boxList;
  ergo_real* result_J_list;
  pthread_t threadID;
  int jobIndexLo, jobIndexHi;
  int maxnoOfMinimalDistrsPerBoxBranch;
  int result;
  int largest_L_used;
};

static void*
execute_J_fmm_worker(void *arg)
{
  JFMMWorkerData *data = static_cast<JFMMWorkerData*>(arg);
  
  data->result = 
    execute_joblist_J_fmm_shared(data->jobIndexLo, data->jobIndexHi,
                                 *data->integralInfo, *data->basisInfo, *data->J_K_params,
                                 data->jobList_J_multipole,
                                 data->boxList,
                                 data->maxnoOfMinimalDistrsPerBoxBranch,
                                 data->result_J_list,
                                 &data->largest_L_used);
  return NULL;
}



/** Compute the FMM part of the Coulomb matrix using threads. 0th
    thread reuses result_J_list, all the other threads need to have temporary
    memory allocated.
*/
static int
execute_joblist_J_fmm_thread(int noOfThreads, int noOfBasisFuncIndexPairs,
                             const IntegralInfo& integralInfo,
                             const BasisInfoStruct & basisInfo,
                             const JK::Params& J_K_params,
                             int noOfJobs_J_multipole,
                             const job_list_multipole_entry_J_struct* jobList_J_multipole,
                             const box_struct *boxList,
                             int maxnoOfMinimalDistrsPerBoxBranch,
                             ergo_real* result_J_list)
{
  std::vector<JFMMWorkerData> threadData(noOfThreads);
  int lastJob = 0;
  Util::TimeMeter timeMeterJmul;

  int th; // declare here because value used after loop after break.
  for(th = 0; th < noOfThreads; th++) {
    threadData[th].integralInfo = &integralInfo;
    threadData[th].basisInfo    = &basisInfo;
    threadData[th].J_K_params   = &J_K_params;
    threadData[th].jobList_J_multipole  = jobList_J_multipole;
    threadData[th].boxList = boxList;
    threadData[th].maxnoOfMinimalDistrsPerBoxBranch = maxnoOfMinimalDistrsPerBoxBranch;
    threadData[th].jobIndexLo = lastJob;
    /* Now we want to compute jobIndexHi as
       (((th+1)*noOfJobs_J_multipole)/noOfThreads) but we need to be
       careful to avoid integer overflow if noOfJobs_J_multipole is
       large. */
    size_t noOfJobs_J_multipole_as_size_t = noOfJobs_J_multipole;
    threadData[th].jobIndexHi = lastJob = ((th+1)*noOfJobs_J_multipole_as_size_t)/noOfThreads;
    threadData[th].largest_L_used = 0;
    if(th) {
      threadData[th].result_J_list = new ergo_real[noOfBasisFuncIndexPairs];
      if( threadData[th].result_J_list == NULL) {
        do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error allocating data for thread %i", th);
        do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
        th++; /* Correct the wait loop upper index. */
        break;
      }
      memset(threadData[th].result_J_list, 0, noOfBasisFuncIndexPairs*sizeof(ergo_real));
    } else {
      threadData[th].result_J_list = result_J_list;
    }
    if(pthread_create(&threadData[th].threadID, NULL,
                      execute_J_fmm_worker, &threadData[th]) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_create for thread %i", th);
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
      th++; /* Correct the wait loop upper index. */
      break;
    }
  }
  int myResult = 0, largest_L_used = 0;
  for(int i = 0; i < th; i++) {
    if(pthread_join(threadData[i].threadID, NULL) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
      myResult++;
    }
    myResult += threadData[i].result;
    if(i == 0) {
      threadData[i].result_J_list = NULL;
    } else {
      for(int idx=0; idx<noOfBasisFuncIndexPairs; idx++)
        result_J_list[idx] += threadData[i].result_J_list[idx];
      delete [] threadData[i].result_J_list;
    }
    if(threadData[i].largest_L_used > largest_L_used)
      largest_L_used = threadData[i].largest_L_used;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "multipole job list for J executed, largest L used: %2i",
	    largest_L_used);
  timeMeterJmul.print(LOG_AREA_INTEGRALS, "Executing multipole job list for J");
  return myResult;
}





/** Computes the Coulomb interaction.
    @param basisInfo
    @param integralInfo
    @param J_K_params the evaluation parameters, thresholds and all.
    @param basisFuncIndexPairList
    @param basisFuncIndexPairCount the length of basisFuncIndexPairList.
    @param D_list basisFuncIndexPairCount elements, with indices
    matching basisFuncIndexPairList.
    @param result_J_list preallocated list that will contain the results.
    @param noOfBasisFuncIndexPairs the length of result_J_list.
    happens to be always equal to basisFuncIndexPairCount
 */
int
compute_J_by_boxes_linear(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  const JK::Params& J_K_params,
			  const basis_func_index_pair_struct* basisFuncIndexPairList,
			  int basisFuncIndexPairCount,
			  const ergo_real* D_list,
			  ergo_real* result_J_list,
			  int noOfBasisFuncIndexPairs)
{
  Util::TimeMeter timeMeterTot;
  
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "entering compute_J_by_boxes_linear, no of basis funcs = %5i, threshold_J = %7.3g", 
	    n, (double)J_K_params.threshold_J);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "use_fmm = %i, fmm_box_size = %6.2f", 
	    J_K_params.use_fmm, (double)J_K_params.fmm_box_size);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "beginning of compute_J_by_boxes_linear");
  
  init_multipole_code();

  Util::TimeMeter timeMeterDistrList;

  // get largest limiting factor
  ergo_real maxLimitingFactor = 0;
  if(get_list_of_labeled_distrs_maxLimitingFactor_linear(basisInfo,
							 integralInfo,
							 J_K_params.threshold_J,
							 basisFuncIndexPairList,
							 basisFuncIndexPairCount,
							 &maxLimitingFactor) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs_maxLimitingFactor_linear");
      return -1;
    }

  // Get number of distributions
  int distrCount = get_list_of_labeled_distrs_linear(basisInfo,
						     integralInfo,
						     J_K_params.threshold_J,
						     NULL,
						     0,
						     maxLimitingFactor,
						     basisFuncIndexPairList,
						     basisFuncIndexPairCount,
						     D_list);
  if(distrCount <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear: (distrCount <= 0)");
      return -1;
    }

  std::vector<DistributionSpecStructLabeled> distrList(distrCount);

  // create list of product primitives, with labels
  int distrCountTemp = get_list_of_labeled_distrs_linear(basisInfo,
							 integralInfo,
							 J_K_params.threshold_J,
							 &distrList[0],
							 distrCount,
							 maxLimitingFactor,
							 basisFuncIndexPairList,
							 basisFuncIndexPairCount,
							 D_list);
  if(distrCountTemp != distrCount)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear: (distrCountTemp != distrCount)");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating list of primitive distributions");

  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(basisFuncIndexPairCount, D_list);

  // compute extent for all distrs
  Util::TimeMeter timeMeterComputeExtentForAllDistrs;
  compute_extent_for_list_of_distributions(distrCount, 
					   &distrList[0], 
					   J_K_params.threshold_J,
					   maxLimitingFactor,
					   maxDensityMatrixElement);
  timeMeterComputeExtentForAllDistrs.print(LOG_AREA_INTEGRALS, "Compute extent for all distrs");

  // get largest and smallest extent
  ergo_real extent_min, extent_max;
  get_largest_and_smallest_extent_for_list_of_distributions(distrCount, &distrList[0], &extent_min, &extent_max);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "extent_min = %8.3f, extent_max = %8.3f", (double)extent_min, (double)extent_max);


  // get maximum number of monomials
  int maxNoOfMonomials = get_max_no_of_monomials_for_list_of_distributions(distrCount, &distrList[0], integralInfo);

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Creating list of distributions etc done, distrCount = %9i", distrCount);
  timeMeterDistrList.print(LOG_AREA_INTEGRALS, "Creating list of distributions etc");


  //
  // This is where we start to worry about the box system
  //

  Util::TimeMeter timeMeterBoxes;
  BoxSystem boxSystem;
  const ergo_real toplevelBoxSize = J_K_params.fmm_box_size;
  if(create_box_system_and_reorder_distrs(distrCount, 
					  &distrList[0],
					  toplevelBoxSize,
					  boxSystem) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system_and_reorder_distrs");
      return -1;
    }



  // Now we have the box system.
  // Create new list of boxes (more advanced boxes this time)

  std::vector<box_struct> boxList(boxSystem.totNoOfBoxes);

  for(int i = 0; i < boxSystem.totNoOfBoxes; i++)
    boxList[i].basicBox = boxSystem.boxList[i];


  int numberOfLevels = boxSystem.noOfLevels;
  int levelCounterList[numberOfLevels];
  int levelStartIndexList[numberOfLevels];
  for(int i = 0; i < numberOfLevels; i++)
    {
      levelCounterList[i] = boxSystem.levelList[i].noOfBoxes;
      levelStartIndexList[i] = boxSystem.levelList[i].startIndexInBoxList;
    }

  
  // OK, boxes created.

  timeMeterBoxes.print(LOG_AREA_INTEGRALS, "Creating boxes");


  int noOfBoxesTopLevel = levelCounterList[numberOfLevels-1];
  box_struct* boxListTopLevel = &boxList[levelStartIndexList[numberOfLevels-1]];
  

  // within each box, divide distrs into branches according to how far they penetrate outside the box.
  ergo_real branchSplitterList[MAX_NO_OF_BRANCHES];

  int noOfBranches = get_branch_splitter_info(branchSplitterList,
					      MAX_NO_OF_BRANCHES,
					      J_K_params,
					      toplevelBoxSize,
					      extent_max);
  if(noOfBranches <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_branch_splitter_info");
      return -1;
    }

  char s[888];
  s[0] = '\0';
  for(int i = 0; i < noOfBranches-1; i++)
    {
      char ss[888];
      sprintf(ss, " %5.2f", (double)branchSplitterList[i]);
      strcat(s, ss);
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "noOfBranches = %i, splitters: %s", noOfBranches, s);


  if(create_branches(noOfBranches,
		     branchSplitterList,
		     distrCount,
		     &distrList[0],
		     noOfBoxesTopLevel,
		     boxListTopLevel) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_branches");
      return -1;
    }


  
  Util::TimeMeter timeMeterJorg;
  int groupCount = 0;
  int maxnoOfMinimalDistrsPerBoxBranch = 0;
  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
	{
	  DistributionSpecStructLabeled* distrListCurrBox = &distrList[boxListTopLevel[i].branchIndexList[branchIndex]];
	  int distrCountCurrBox = boxListTopLevel[i].branchCountList[branchIndex];

	  if(organize_distributions(integralInfo,
				    distrListCurrBox, 
				    distrCountCurrBox,
				    &boxListTopLevel[i].branchList[branchIndex].org,
				    boxListTopLevel[i].basicBox.centerCoords,
				    boxListTopLevel[i].basicBox.width) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in organize_distributions for box %i branch %i", i, branchIndex);
	      return -1;
	    }
	  groupCount += boxListTopLevel[i].branchList[branchIndex].org.groupCount;
	  if(boxListTopLevel[i].branchList[branchIndex].org.minimalDistrCount > maxnoOfMinimalDistrsPerBoxBranch)
	    maxnoOfMinimalDistrsPerBoxBranch = boxListTopLevel[i].branchList[branchIndex].org.minimalDistrCount;
	} // END FOR branchIndex
    } // END FOR i
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "J org done, groupCount = %8i", groupCount);
  timeMeterJorg.print(LOG_AREA_INTEGRALS, "J org");
  

  // Generate multipole for each group, and find center-of-charge for each branch
  Util::TimeMeter timeMeterGenerateGr;
  std::vector<multipole_struct_small> multipoleListForGroups(groupCount);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating multipoleListForGroups");

  int count = 0;
  ergo_real totChargeWholeSystem = 0;
  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      ergo_real centerOfChargeList[3];
      ergo_real averagePosList[3];
      for(int kk = 0; kk < 3; kk++)
	{
	  centerOfChargeList[kk] = 0;
	  averagePosList[kk] = 0;
	}
      int avgPosCounter = 0;
      for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
	{
	  ergo_real chargeSum = 0;

	  const chunk_struct* chunkList = &boxListTopLevel[i].branchList[branchIndex].org.chunkList[0];
	  const cluster_struct* clusterList = &boxListTopLevel[i].branchList[branchIndex].org.clusterList[0];
	  distr_group_struct* groupList = &boxListTopLevel[i].branchList[branchIndex].org.groupList[0];
	  const minimal_distr_struct* minimalDistrList = &boxListTopLevel[i].branchList[branchIndex].org.minimalDistrList[0];
	  int chunkCount = boxListTopLevel[i].branchList[branchIndex].org.chunkCount;
	  const basis_func_pair_struct* basisFuncPairList = &boxListTopLevel[i].branchList[branchIndex].org.basisFuncPairList[0];

	  ergo_real* maxMomentVectorNormForDistrsList = boxListTopLevel[i].branchList[branchIndex].maxMomentVectorNormForDistrsList;
	  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	    maxMomentVectorNormForDistrsList[l] = 0;

	  for(int chunkIndex = 0; chunkIndex < chunkCount; chunkIndex++)
	    {
	      int clusterCount = chunkList[chunkIndex].noOfClusters;
	      int cluster_start = chunkList[chunkIndex].clusterStartIndex;
	      for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++)
		{
		  int group_start = clusterList[clusterIndex].groupStartIndex;
		  int group_end = group_start + clusterList[clusterIndex].noOfGroups;
		  for(int groupIndex = group_start; groupIndex < group_end; groupIndex++)
		    {
		      distr_group_struct* currGroup = &groupList[groupIndex];

		      // Now create a single multipole description of the density of this group.
		      multipole_struct_small* multipoleCurrGroup = &multipoleListForGroups[count];
		      currGroup->multipolePtr = multipoleCurrGroup;
		      multipoleCurrGroup->degree = -1;
		      multipoleCurrGroup->noOfMoments = 0;
		      multipoleCurrGroup->centerCoords[0] = currGroup->centerCoords[0];
		      multipoleCurrGroup->centerCoords[1] = currGroup->centerCoords[1];
		      multipoleCurrGroup->centerCoords[2] = currGroup->centerCoords[2];
		      memset(multipoleCurrGroup->momentList, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC*sizeof(ergo_real));
		      count++;

		      int distr_start = currGroup->startIndex;
		      int distr_end = distr_start + currGroup->distrCount;
		      for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++)
			{
			  int basisFuncPairIndex = minimalDistrList[distrIndex].basisFuncPairIndex;
			  int monomialIndex = minimalDistrList[distrIndex].monomialIndex;
			  ergo_real coeff = minimalDistrList[distrIndex].coeff;
			  // get monomialInts from monomialIndex
			  DistributionSpecStruct distr;
			  distr.monomialInts[0] = integralInfo.monomial_info.monomial_list[monomialIndex].ix;
			  distr.monomialInts[1] = integralInfo.monomial_info.monomial_list[monomialIndex].iy;
			  distr.monomialInts[2] = integralInfo.monomial_info.monomial_list[monomialIndex].iz;
			  distr.coeff = coeff;
			  distr.exponent = currGroup->exponent;
			  distr.centerCoords[0] = currGroup->centerCoords[0];
			  distr.centerCoords[1] = currGroup->centerCoords[1];
			  distr.centerCoords[2] = currGroup->centerCoords[2];

			  multipole_struct_small multipole;
			  if(compute_multipole_moments(integralInfo, &distr, &multipole) != 0)
			    {
			      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_multipole_moments");
			      return -1;
			    }

			  // add this multipole to multipole for group.
			  int a = basisFuncPairList[chunkList[chunkIndex].basisFuncPairListIndex+basisFuncPairIndex].index_1;
			  int b = basisFuncPairList[chunkList[chunkIndex].basisFuncPairListIndex+basisFuncPairIndex].index_2;
			  ergo_real factor = basisFuncPairList[chunkList[chunkIndex].basisFuncPairListIndex+basisFuncPairIndex].dmatElement;
			  if(a != b)
			    factor *= 2;

			  for(int l = 0; l <= multipole.degree; l++)
			    {
			      int startIndex = l*l;
			      int endIndex = (l+1)*(l+1);
			      ergo_real sum = 0;
			      for(int A = startIndex; A < endIndex; A++)
				sum += multipole.momentList[A]*multipole.momentList[A];
			      ergo_real subNorm = std::sqrt(sum);
			      if(subNorm > maxMomentVectorNormForDistrsList[l])
				maxMomentVectorNormForDistrsList[l] = subNorm;
			    }
			  
			  for(int kk = 0; kk < multipole.noOfMoments; kk++)
			    multipoleCurrGroup->momentList[kk] += factor * multipole.momentList[kk];
			  if(multipole.degree > multipoleCurrGroup->degree)
			    multipoleCurrGroup->degree = multipole.degree;
			  if(multipole.noOfMoments > multipoleCurrGroup->noOfMoments)
			    multipoleCurrGroup->noOfMoments = multipole.noOfMoments;
			} // END FOR distrIndex

		      // OK, multipoleCurrGroup is complete.
		      chargeSum += multipoleCurrGroup->momentList[0];
		      for(int kk = 0; kk < 3; kk++)
			{
			  centerOfChargeList[kk] += multipoleCurrGroup->centerCoords[kk] * multipoleCurrGroup->momentList[0];
			  averagePosList[kk] += multipoleCurrGroup->centerCoords[kk];
			}
		      avgPosCounter++;

		    } // END FOR groupIndex
		} // END FOR clusterIndex
	    } // END FOR chunkIndex

	  totChargeWholeSystem += chargeSum;
	  boxListTopLevel[i].branchList[branchIndex].totCharge = chargeSum;
	} // END FOR branchIndex

      // use average position instead of center-of-charge, 
      // because center-of-charge is ill-defined when some charges are negative.
      if(avgPosCounter == 0)
	{
	  for(int kk = 0; kk < 3; kk++)
	    boxListTopLevel[i].multipolePoint[kk] = boxListTopLevel[i].basicBox.centerCoords[kk];
	}
      else
	{
	  for(int kk = 0; kk < 3; kk++)
	    boxListTopLevel[i].multipolePoint[kk] = averagePosList[kk] / avgPosCounter;
	}

      // check that "multipolePoint" is not too far from box center.
      ergo_real sumofsquares = 0;
      for(int kk = 0; kk < 3; kk++)
	{
	  ergo_real dx = boxListTopLevel[i].multipolePoint[kk] - boxListTopLevel[i].basicBox.centerCoords[kk];
	  sumofsquares += dx*dx;
	}
      ergo_real distFromCenter = std::sqrt(sumofsquares);
      if(distFromCenter > boxListTopLevel[i].basicBox.width)
	{
	  printf("error: (distFromCenter > boxListTopLevel[i].width)\n");
	  printf("distFromCenter = %22.11f\n", (double)distFromCenter);
	  printf("boxListTopLevel[i].basicBox.width = %22.11f\n", (double)boxListTopLevel[i].basicBox.width);
	  printf("avgPosCounter = %i\n", avgPosCounter);
	  exit(0);
	  return -1;
	}

      // Now we have set boxListTopLevel[i].multipolePoint
      // Copy it to each branch
      for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
	{
	  for(int kk = 0; kk < 3; kk++)
	    boxListTopLevel[i].branchList[branchIndex].multipolePoint[kk] = boxListTopLevel[i].multipolePoint[kk];
	}

    } // END FOR i

  timeMeterGenerateGr.print(LOG_AREA_INTEGRALS, "Generate group multipoles");
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "totChargeWholeSystem = %22.11f", (double)totChargeWholeSystem);

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Generating multipole for each branch at top level, MAX_MULTIPOLE_DEGREE = %2i", (int)MAX_MULTIPOLE_DEGREE);

  // Generate multipole for each branch at top level (smallest boxes)
  Util::TimeMeter timeMeterTranslate1;
  MMTranslator translator;
  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
	{
	  ergo_real* multipolePointCoords = boxListTopLevel[i].branchList[branchIndex].multipolePoint;
	  multipole_struct_large branchMultipole;
	  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
	    branchMultipole.momentList[A] = 0;
	  for(int kk = 0; kk < 3; kk++)
	    branchMultipole.centerCoords[kk] = multipolePointCoords[kk];
	  branchMultipole.degree = MAX_MULTIPOLE_DEGREE;
	  branchMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

	  const chunk_struct* chunkList = &boxListTopLevel[i].branchList[branchIndex].org.chunkList[0];
	  const cluster_struct* clusterList = &boxListTopLevel[i].branchList[branchIndex].org.clusterList[0];
	  const distr_group_struct* groupList = &boxListTopLevel[i].branchList[branchIndex].org.groupList[0];
	  int chunkCount = boxListTopLevel[i].branchList[branchIndex].org.chunkCount;
	  for(int chunkIndex = 0; chunkIndex < chunkCount; chunkIndex++)
	    {
	      int clusterCount = chunkList[chunkIndex].noOfClusters;
	      int cluster_start = chunkList[chunkIndex].clusterStartIndex;
	      for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++)
		{
		  int group_start = clusterList[clusterIndex].groupStartIndex;
		  int group_end = group_start + clusterList[clusterIndex].noOfGroups;
		  for(int groupIndex = group_start; groupIndex < group_end; groupIndex++)
		    {
		      const distr_group_struct* currGroup = &groupList[groupIndex];

		      // take multipole for this group, and translate it to center-of-charge point
		      ergo_real dx = currGroup->multipolePtr->centerCoords[0] - multipolePointCoords[0];
		      ergo_real dy = currGroup->multipolePtr->centerCoords[1] - multipolePointCoords[1];
		      ergo_real dz = currGroup->multipolePtr->centerCoords[2] - multipolePointCoords[2];

		      ergo_real W[MAX_NO_OF_MOMENTS_PER_MULTIPOLE*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
		      translator.getTranslationMatrix
                        (dx, dy, dz, MAX_MULTIPOLE_DEGREE,
                         currGroup->multipolePtr->degree, W);

		      multipole_struct_large translatedMultipole;
		      for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
			{
			  ergo_real sum = 0;
			  for(int B = 0; B < currGroup->multipolePtr->noOfMoments; B++)
			    sum += W[A*currGroup->multipolePtr->noOfMoments+B] * currGroup->multipolePtr->momentList[B];
			  translatedMultipole.momentList[A] = sum;
			} // END FOR A
		      for(int kk = 0; kk < 3; kk++)
			translatedMultipole.centerCoords[kk] = multipolePointCoords[kk];
		      translatedMultipole.degree = MAX_MULTIPOLE_DEGREE;
		      translatedMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

		      // add translated multipole to branch multipole
		      for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
			branchMultipole.momentList[A] += translatedMultipole.momentList[A];
		    } // END FOR groupIndex
		} // END FOR clusterIndex
	    } // END FOR chunkIndex
	  setup_multipole_maxAbsMomentList(&branchMultipole);
	  boxListTopLevel[i].branchList[branchIndex].multipole = branchMultipole;
	} // END FOR branchIndex
    } // END FOR i
  timeMeterTranslate1.print(LOG_AREA_INTEGRALS, "Translate multipoles (step 1)");


  // OK, multipoles created for top level.
  // Now go through the other levels, joining multipoles from child boxes to a single multipole (one per branch) in parent box

  Util::TimeMeter timeMeterTranslate2;
  for(int levelNumber = numberOfLevels-2; levelNumber >= 0; levelNumber--)
    {
      int noOfBoxesCurrLevel = levelCounterList[levelNumber];

      box_struct* boxListCurrLevel = &boxList[levelStartIndexList[levelNumber]];
      for(int boxIndex = 0; boxIndex < noOfBoxesCurrLevel; boxIndex++)
	{
	  box_struct* currBox = &boxListCurrLevel[boxIndex];
	  int noOfChildren = currBox->basicBox.noOfChildBoxes;

	  if(noOfChildren == 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ERROR: (noOfChildren == 0)");
	      return -1;
	    }

	  for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
	    {
	      multipole_struct_large* newMultipole = &currBox->branchList[branchIndex].multipole;
	      for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
		newMultipole->momentList[A] = 0;

	      // get average position of child multipoles
	      ergo_real avgPosList[3];
	      for(int kk = 0; kk < 3; kk++)
		avgPosList[kk] = 0;

	      ergo_real* maxMomentVectorNormForDistrsList = currBox->branchList[branchIndex].maxMomentVectorNormForDistrsList;
	      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		maxMomentVectorNormForDistrsList[l] = 0;

	      for(int childIndex = 0; childIndex < noOfChildren; childIndex++)
		{
		  int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
		  box_struct* childBox = &boxList[childIndexInBoxList];
		  for(int kk = 0; kk < 3; kk++)
		    avgPosList[kk] += childBox->branchList[branchIndex].multipole.centerCoords[kk];
		} // END FOR childIndex

	      for(int kk = 0; kk < 3; kk++)
		newMultipole->centerCoords[kk] = avgPosList[kk] / noOfChildren;
	      newMultipole->degree = MAX_MULTIPOLE_DEGREE;
	      newMultipole->noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

	      // We also want to get maxExtent and maxDistanceOutsideBox for parent box (use largest values found among the children).
	      ergo_real maxExtent = 0;
	      ergo_real maxDistanceOutsideBox = 0;

	      // Now translate child multipoles and add to parent multipole
	      for(int childIndex = 0; childIndex < noOfChildren; childIndex++)
		{
		  int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
		  box_struct* childBox = &boxList[childIndexInBoxList];
		  multipole_struct_large* childMultipole = &childBox->branchList[branchIndex].multipole;

		  if(childBox->branchList[branchIndex].org.maxExtent > maxExtent)
		    maxExtent = childBox->branchList[branchIndex].org.maxExtent;

		  if(childBox->branchList[branchIndex].org.maxDistanceOutsideBox > maxDistanceOutsideBox)
		    maxDistanceOutsideBox = childBox->branchList[branchIndex].org.maxDistanceOutsideBox;

		  ergo_real dx = childMultipole->centerCoords[0] - newMultipole->centerCoords[0];
		  ergo_real dy = childMultipole->centerCoords[1] - newMultipole->centerCoords[1];
		  ergo_real dz = childMultipole->centerCoords[2] - newMultipole->centerCoords[2];

		  ergo_real W[MAX_NO_OF_MOMENTS_PER_MULTIPOLE*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
		  translator.getTranslationMatrix(dx, dy, dz,
                                                  MAX_MULTIPOLE_DEGREE,
                                                  MAX_MULTIPOLE_DEGREE, W);

		  multipole_struct_large translatedMultipole;
		  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
		    {
		      ergo_real sum = 0;
		      for(int B = 0; B < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; B++)
			sum += W[A*MAX_NO_OF_MOMENTS_PER_MULTIPOLE+B] * childMultipole->momentList[B];
		      translatedMultipole.momentList[A] = sum;
		    } // END FOR A
		  for(int kk = 0; kk < 3; kk++)
		    translatedMultipole.centerCoords[kk] = newMultipole->centerCoords[kk];
		  translatedMultipole.degree = MAX_MULTIPOLE_DEGREE;
		  translatedMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

		  // add translated multipole to parent multipole
		  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
		    newMultipole->momentList[A] += translatedMultipole.momentList[A];
		  
		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    {
		      ergo_real childValue = childBox->branchList[branchIndex].maxMomentVectorNormForDistrsList[l];
		      if(childValue > maxMomentVectorNormForDistrsList[l])
			maxMomentVectorNormForDistrsList[l] = childValue;
		    }
		  
		} // END FOR childIndex

	      setup_multipole_maxAbsMomentList(newMultipole);
 
	      currBox->branchList[branchIndex].org.maxExtent = maxExtent;
	      currBox->branchList[branchIndex].org.maxDistanceOutsideBox = maxDistanceOutsideBox;

	    } // END FOR branchIndex
	} // END FOR boxIndex
    } // END FOR levelNumber
  timeMeterTranslate2.print(LOG_AREA_INTEGRALS, "Translate multipoles (step 2)");

  // Set J to zero
  memset(result_J_list, 0, basisFuncIndexPairCount*sizeof(ergo_real));


  // set up list of upper limits for interaction matrix elements
  Util::TimeMeter timeMeterIntMatLimits;
  ergo_real maxDistance = getSafeMaxDistance(basisInfo);
  mm_limits_init(maxDistance);
  timeMeterIntMatLimits.print(LOG_AREA_INTEGRALS, "mm_limits_init");


  // Create job lists for J

  Util::TimeMeter timeMeterJjoblist;
  int noOfJobs_J_standard_firstCount = 0;
  int noOfJobs_J_multipole_firstCount = 0;
  for(int branch_i = 0; branch_i < noOfBranches; branch_i++)
    {
      for(int branch_j = 0; branch_j < noOfBranches; branch_j++)
	{
	  int noOfNewJobs_standard = 0;
	  int noOfNewJobs_multipole = 0;
	  if(get_joblists_J_for_two_boxes_recursive(basisInfo,
						    integralInfo,
						    maxNoOfMonomials,
						    J_K_params.threshold_J,
						    &boxList[0],
						    numberOfLevels,
						    0,
						    0,
						    0,
						    branch_i,
						    branch_j,
						    NULL,
						    HUGE_INTEGER_NUMBER,
						    &noOfNewJobs_standard,
						    NULL,
						    HUGE_INTEGER_NUMBER,
						    &noOfNewJobs_multipole
						    ) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive");
	      return -1;
	    }
	  noOfJobs_J_standard_firstCount += noOfNewJobs_standard;
	  noOfJobs_J_multipole_firstCount += noOfNewJobs_multipole;
	}
    }
  std::vector<job_list_standard_entry_J_struct> jobList_J_standard(noOfJobs_J_standard_firstCount);
  std::vector<job_list_multipole_entry_J_struct> jobList_J_multipole(noOfJobs_J_multipole_firstCount);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating jobLists for J");

  int noOfJobs_J_standard = 0;
  int noOfJobs_J_multipole = 0;
  for(int branch_i = 0; branch_i < noOfBranches; branch_i++)
    {
      for(int branch_j = 0; branch_j < noOfBranches; branch_j++)
	{
	  int noOfNewJobs_standard = 0;
	  int noOfNewJobs_multipole = 0;
	  if(get_joblists_J_for_two_boxes_recursive(basisInfo,
						    integralInfo,
						    maxNoOfMonomials,
						    J_K_params.threshold_J,
						    &boxList[0],
						    numberOfLevels,
						    0,
						    0,
						    0,
						    branch_i,
						    branch_j,
						    &jobList_J_standard[noOfJobs_J_standard],
						    noOfJobs_J_standard_firstCount - noOfJobs_J_standard,
						    &noOfNewJobs_standard,
						    &jobList_J_multipole[noOfJobs_J_multipole],
						    noOfJobs_J_multipole_firstCount - noOfJobs_J_multipole,
						    &noOfNewJobs_multipole
						    ) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive");
	      return -1;
	    }
	  noOfJobs_J_standard += noOfNewJobs_standard;
	  noOfJobs_J_multipole += noOfNewJobs_multipole;
	}
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "job lists for J created OK, noOfJobs_J_standard = %8i, noOfJobs_J_multipole = %8i",
	    noOfJobs_J_standard, noOfJobs_J_multipole);
  timeMeterJjoblist.print(LOG_AREA_INTEGRALS, "Creating job lists for J");

  // Execute standard job list for J

  int noOfThreads = J_K_params.noOfThreads_J;
  
  if(noOfThreads <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear: (noOfThreads <= 0)");
      return -1;
    }
  if(noOfThreads == 1)
    {
      if(execute_joblist_J_std_serial(noOfJobs_J_standard,
				      &jobList_J_standard[0],
				      basisInfo,
				      integralInfo,
				      maxNoOfMonomials,
				      result_J_list,
				      &boxList[0],
				      J_K_params.threshold_J) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_std_serial");
	  return -1;
	}
    }
  else
    {
      if(execute_joblist_J_std_threaded(noOfThreads,
					noOfJobs_J_standard,
					&jobList_J_standard[0],
					basisInfo,
					integralInfo,
					maxNoOfMonomials,
					result_J_list,
					noOfBasisFuncIndexPairs,
					&boxList[0],
					J_K_params.threshold_J) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_std_threaded");
	  return -1;
	}
    }


  // sort multipole job list by boxindex and branchindex.
  Util::TimeMeter timeMeterJmulSort;
  if(sort_list_of_multipole_jobs(&jobList_J_multipole[0], noOfJobs_J_multipole) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in sort_list_of_multipole_jobs");
      return -1;
    }
  timeMeterJmulSort.print(LOG_AREA_INTEGRALS, "sort_list_of_multipole_jobs");


  /* Execute multipole list */
  if(noOfThreads == 1)
    {
      if( execute_joblist_J_fmm_serial(integralInfo, basisInfo, J_K_params,
                                       noOfJobs_J_multipole, &jobList_J_multipole[0],
                                       &boxList[0], maxnoOfMinimalDistrsPerBoxBranch,
                                       result_J_list) != 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_fmm_serial");
          return -1;
        }
    }
  else
    {
      if( execute_joblist_J_fmm_thread(noOfThreads, noOfBasisFuncIndexPairs,
                                       integralInfo, basisInfo, J_K_params,
                                       noOfJobs_J_multipole, &jobList_J_multipole[0],
                                       &boxList[0], maxnoOfMinimalDistrsPerBoxBranch,
                                       result_J_list) != 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_fmm_thread");
          return -1;
        }
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_J_by_boxes_linear ending OK.");
  timeMeterTot.print(LOG_AREA_INTEGRALS, "compute_J_by_boxes_linear");
  
  return 0;
}




int
compute_J_by_boxes(const BasisInfoStruct & basisInfo,
		   const IntegralInfo & integralInfo,
		   const JK::Params& J_K_params,
		   ergo_real* J,
		   const ergo_real* dens)
{
  int n = basisInfo.noOfBasisFuncs;

  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(n*n, dens);

  std::vector<basis_func_index_pair_struct> basisFuncIndexPairList;

  int noOfBasisFuncIndexPairs = get_basis_func_pair_list_2el(basisInfo,
							     integralInfo,
							     J_K_params.threshold_J,
							     maxDensityMatrixElement,
							     basisFuncIndexPairList);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basis_func_pair_list");
    return -1;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "noOfBasisFuncIndexPairs = %i ==> storing %6.2f %% of a full matrix", 
	    noOfBasisFuncIndexPairs, (double)100*noOfBasisFuncIndexPairs/((double)n*n));

  std::vector<ergo_real> D_list(noOfBasisFuncIndexPairs);
  std::vector<ergo_real> J_list(noOfBasisFuncIndexPairs);

  // Setup D_list
  for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
    {
      int a = basisFuncIndexPairList[i].index_1;
      int b = basisFuncIndexPairList[i].index_2;
      D_list[i] = dens[a*n+b];
    }
  
  if(compute_J_by_boxes_linear(basisInfo,
			       integralInfo,
			       J_K_params,
			       &basisFuncIndexPairList[0],
			       noOfBasisFuncIndexPairs,
			       &D_list[0],
			       &J_list[0],
			       noOfBasisFuncIndexPairs) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear");
      return -1;
    }

  // Now transfer result from J_list to J
  // First set all of J to zero (otherwise there may be something in the non-relevant part of the matrix)
  for(int i = 0; i < n*n; i++)
    J[i] = 0;
  for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
    {
      int a = basisFuncIndexPairList[i].index_1;
      int b = basisFuncIndexPairList[i].index_2;
      J[a*n+b] = J_list[i];
      J[b*n+a] = J_list[i];
    }

  return 0;
}



/*
compute_J_by_boxes_nosymm does the same as compute_J_by_boxes,
but without assuming the density matrix to be symmetric.
*/
int
compute_J_by_boxes_nosymm(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  const JK::Params& J_K_params,
			  ergo_real* J,
			  const ergo_real* dens)
{
  int n = basisInfo.noOfBasisFuncs;

  // Create symmetrized density matrix P_ab = 0.5*(D_ab + D_ba)
  std::vector<ergo_real> P(n*n);
  for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
      P[a*n+b] = 0.5 * ( dens[a*n+b] + dens[b*n+a] );

  int rc = compute_J_by_boxes(basisInfo,
			      integralInfo,
			      J_K_params,
			      J,
			      &P[0]);
  return rc;
}






