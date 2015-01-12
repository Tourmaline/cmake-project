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
#include <memory.h>
#include <algorithm>

#include "organize_distrs.h"
#include "pi.h"

#include <cstdio>


static void
do_sort_int_list(int* list, int n)
{
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n-i-1; j++)
      {
	if(list[j+1] < list[j])
	  {
	    int temp = list[j];
	    list[j] = list[j+1];
	    list[j+1] = temp;
	  }
      } // END FOR i j
}



static void get_conversion_matrix_for_group(
				    const IntegralInfo & integralInfo,
				    const distr_group_struct & group,
				    int n1max,
				    const minimal_distr_struct* minimalDistrList_1,
				    int noOfBasisFuncPairs_1, 
				    const i_j_val_struct* convMat1_sp,
				    int convMat1_nnz,
				    i_j_val_struct* BB1_x_Ai1_x_convMat1_sp_result, // result
				    int & BB1_x_Ai1_x_convMat1_nnz_result)  // result
{
  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1max];
  int distrCount_i = group.distrCount;
  i_j_val_struct Ai1_sp[distrCount_i];
  for(int i = 0; i < distrCount_i; i++) {
    int monomialIndex = minimalDistrList_1[group.startIndex+i].monomialIndex;
    ergo_real value = minimalDistrList_1[group.startIndex+i].coeff;
    Ai1_sp[i].i = i;
    Ai1_sp[i].j = monomialIndex;
    Ai1_sp[i].value = value;
  }
  i_j_val_struct BB1_sp[distrCount_i];
  for(int kk = 0; kk < distrCount_i; kk++) {
    int idx = minimalDistrList_1[kk+group.startIndex].basisFuncPairIndex;
    BB1_sp[kk].i = idx;
    BB1_sp[kk].j = kk;
    BB1_sp[kk].value = 1;
  }
  // Multiply Ai1 by convMat1. Dimensions: (distrCount_i*noOfMonomials_1) x (noOfMonomials_1*noOfMonomials_1)
  i_j_val_struct* Ai1_x_convMat1_sp = new i_j_val_struct[distrCount_i*noOfMonomials_1];
  int Ai1_x_convMat1_nnz = spmat_multiply_matrices(Ai1_sp, distrCount_i, convMat1_sp, convMat1_nnz, Ai1_x_convMat1_sp, distrCount_i, noOfMonomials_1);
  // Multiply BB1 by Ai1_x_convMat1. Dimensions: (noOfBasisFuncPairs_1*distrCount_i) x (distrCount_i*noOfMonomials_1)
  BB1_x_Ai1_x_convMat1_nnz_result = spmat_multiply_matrices(BB1_sp, distrCount_i, Ai1_x_convMat1_sp, Ai1_x_convMat1_nnz, BB1_x_Ai1_x_convMat1_sp_result, noOfBasisFuncPairs_1, noOfMonomials_1);
  delete [] Ai1_x_convMat1_sp;
}



int
organize_distributions(const IntegralInfo & integralInfo,
		       DistributionSpecStructLabeled* distrList_in, 
		       int distrCount, 
		       distr_org_struct* result,
		       const ergo_real* boxCenterCoords,
		       ergo_real boxWidth)
{
  std::vector<DistributionSpecStructLabeled> distrList(distrCount);

  // sort list of distributions by center, type and exponent
  // first group the ones that have same center and same exponent.
  std::vector<int> groupCountList(distrCount);
  std::vector<int> groupIndexList(distrCount);


  // start by bucket sort based on "best" coordinate.
  const ergo_real HUGE_NUMBER = 888888888;
  ergo_real xminList[3];
  ergo_real xmaxList[3];
  ergo_real xdiffList[3];
  for(int kk = 0; kk < 3; kk++)
    {
      xminList[kk] =  HUGE_NUMBER;
      xmaxList[kk] = -HUGE_NUMBER;
    }
  for(int i = 0; i < distrCount; i++)
    {
      for(int kk = 0; kk < 3; kk++)
	{
	  ergo_real x = distrList_in[i].distr.centerCoords[kk];
	  if(x < xminList[kk])
	    xminList[kk] = x;
	  if(x > xmaxList[kk])
	    xmaxList[kk] = x;
	}
    } // END FOR i
  int bestCoordIndex = 0;
  for(int kk = 0; kk < 3; kk++)
    {
      xdiffList[kk] = xmaxList[kk] - xminList[kk];
      if(xdiffList[kk] > xdiffList[bestCoordIndex])
	bestCoordIndex = kk;
    }
#define NO_OF_SORT_BUCKETS 30
  ergo_real splitterList[NO_OF_SORT_BUCKETS-1];
  for(int i = 0; i < NO_OF_SORT_BUCKETS-1; i++)
    splitterList[i] = xminList[bestCoordIndex] + ((ergo_real)i + 1) * xdiffList[bestCoordIndex] / NO_OF_SORT_BUCKETS;
  int* bucketList[NO_OF_SORT_BUCKETS];
  int bucketCounterList[NO_OF_SORT_BUCKETS];
  for(int i = 0; i < NO_OF_SORT_BUCKETS; i++)
    {
      bucketList[i] = new int[distrCount];
      bucketCounterList[i] = 0;
    }
  for(int i = 0; i < distrCount; i++)
    {
      int bucketIndex = -1;
      for(int j = 0; j < NO_OF_SORT_BUCKETS-1; j++)
	{
	  if(distrList_in[i].distr.centerCoords[bestCoordIndex] < splitterList[j])
	    {
	      bucketIndex = j;
	      break;
	    }
	}
      if(bucketIndex == -1)
	bucketIndex = NO_OF_SORT_BUCKETS-1;
      bucketList[bucketIndex][bucketCounterList[bucketIndex]] = i;
      bucketCounterList[bucketIndex]++;
    } // END FOR i

  int destCount = 0;
  int groupCount = 0;

  // create groups for one bucket at a time
  for(int bucketIndex = 0; bucketIndex < NO_OF_SORT_BUCKETS; bucketIndex++)
    {
      int nLeft = bucketCounterList[bucketIndex];
      while(nLeft > 0)
	{
	  int i = 0;
	  int remainingIndex = 0;
	  int destCountSaved = destCount;
	  distrList[destCount] = distrList_in[bucketList[bucketIndex][i]];
	  destCount++;
	  // now find all that belong to same group
	  for(int k = i+1; k < nLeft; k++)
	    {
	      ergo_real dx, dy, dz;
	      dx = distrList_in[bucketList[bucketIndex][k]].distr.centerCoords[0] - distrList[destCountSaved].distr.centerCoords[0];
	      dy = distrList_in[bucketList[bucketIndex][k]].distr.centerCoords[1] - distrList[destCountSaved].distr.centerCoords[1];
	      dz = distrList_in[bucketList[bucketIndex][k]].distr.centerCoords[2] - distrList[destCountSaved].distr.centerCoords[2];
	      ergo_real r2 = dx*dx + dy*dy + dz*dz;
	      ergo_real absExponentDiff = distrList_in[bucketList[bucketIndex][k]].distr.exponent - distrList[destCountSaved].distr.exponent;
	      if(absExponentDiff < 0)
		absExponentDiff *= -1;
	      if(absExponentDiff < 1e-11 && r2 < 1e-10)
		{
		  // OK, close enough, we regard this as being same center and same exponent.
		  // add to distrList, and remove from distrList_in.
		  distrList[destCount] = distrList_in[bucketList[bucketIndex][k]];
		  destCount++;
		}
	      else
		{
		  // no, different center or exponent
		  if(remainingIndex != k)
		    bucketList[bucketIndex][remainingIndex] = bucketList[bucketIndex][k];
		  remainingIndex++;
		}
	    } // END FOR k find all that belong to same group      
	  int noOfDistrsInGroup = destCount - destCountSaved;
	  nLeft -= noOfDistrsInGroup;
	  groupCountList[groupCount] = noOfDistrsInGroup;
	  groupCount++;
	  if(remainingIndex == 0)
	    break;
	} // END WHILE group the ones that have same center and same exponent.      
    } // END FOR bucketIndex

  for(int i = 0; i < NO_OF_SORT_BUCKETS; i++)
    {
      delete [] bucketList[i];
      bucketList[i] = NULL;
    }

  // set groupIndexList
  int currGroupIndex = 0;
  for(int i = 0; i < groupCount; i++)
    {
      groupIndexList[i] = currGroupIndex;
      currGroupIndex += groupCountList[i];
    }

  // Set groupID
  for(int i = 0; i < groupCount; i++)
    {
      DistributionSpecStructLabeled* groupPtr = &distrList[groupIndexList[i]];
      int currCount = groupCountList[i];
      for(int j = 0; j < currCount; j++)
	groupPtr[j].groupID = i + 1;
    } // END FOR i 

  // Within each group, sort by monomialInts and basisFuncIndeces
  for(int i = 0; i < groupCount; i++)
    {
      DistributionSpecStructLabeled* groupPtr = &distrList[groupIndexList[i]];
      int currCount = groupCountList[i];
      for(int k = 0; k < currCount; k++)
	for(int m = 0; m < currCount - 1 - k; m++)
	  {
	    int doSwitch = 0;
	    if(doSwitch == 0 && groupPtr[m].distr.monomialInts[0] > groupPtr[m+1].distr.monomialInts[0])
	      doSwitch = 1;
	    else
	      doSwitch = -1;
	    if(doSwitch == 0 && groupPtr[m].distr.monomialInts[1] > groupPtr[m+1].distr.monomialInts[1])
	      doSwitch = 1;
	    else
	      doSwitch = -1;
	    if(doSwitch == 0 && groupPtr[m].distr.monomialInts[2] > groupPtr[m+1].distr.monomialInts[2])
	      doSwitch = 1;
	    else
	      doSwitch = -1;
	    if(doSwitch == 0 && groupPtr[m].basisFuncIndex_1 > groupPtr[m+1].basisFuncIndex_1)
	      doSwitch = 1;
	    else
	      doSwitch = -1;
	    if(doSwitch == 0 && groupPtr[m].basisFuncIndex_2 > groupPtr[m+1].basisFuncIndex_2)
	      doSwitch = 1;
	    else
	      doSwitch = -1;
	    if(doSwitch == 1)
	      {
		// switch
		DistributionSpecStructLabeled temp;
		temp = groupPtr[m];
		groupPtr[m] = groupPtr[m+1];
		groupPtr[m+1] = temp;
	      }
	  } // END FOR k m
    } // END FOR i

  
  result->groupList.resize(groupCount);
  distr_group_struct* groupList = &result->groupList[0];

  for(int i = 0; i < groupCount; i++)
    {
      groupList[i].distrCount = groupCountList[i];
      groupList[i].startIndex = groupIndexList[i];
      // get nmax
      int nmax = 0;
      for(int ii = groupIndexList[i]; ii < groupIndexList[i] + groupCountList[i]; ii++)
	{
	  int sum = 0;
	  for(int kk = 0; kk < 3; kk++)
	    sum += distrList[ii].distr.monomialInts[kk];
	  if(sum > nmax)
	    nmax = sum;
	}
      groupList[i].nmax = nmax;
      // get centerCoords and exponent
      for(int ii = 0; ii < 3; ii++)
	groupList[i].centerCoords[ii] = distrList[groupIndexList[i]].distr.centerCoords[ii];
      groupList[i].exponent = distrList[groupIndexList[i]].distr.exponent;
      // get maxSize, maxLimitingFactor, maxExtent for this group.
      ergo_real maxSize = 0;
      ergo_real maxLimitingFactor = 0;
      ergo_real maxExtent = 0;
      for(int ii = groupIndexList[i]; ii < groupIndexList[i] + groupCountList[i]; ii++)
	{
	  ergo_real size = std::fabs(std::pow((ergo_real)pi/distrList[ii].distr.exponent, (ergo_real)1.5) * distrList[ii].distr.coeff);
	  if(size > maxSize)
	    maxSize = size;
	  ergo_real limitingFactor = distrList[ii].limitingFactor;
	  if(limitingFactor > maxLimitingFactor)
	    maxLimitingFactor = limitingFactor;
	  ergo_real extent = distrList[ii].distr.extent;
	  if(extent > maxExtent)
	    maxExtent = extent;
	}
      groupList[i].maxSizeGroup = maxSize;
      groupList[i].maxLimitingFactorGroup = maxLimitingFactor;
      groupList[i].maxExtentGroup = maxExtent;

      // Get maxAbsDmatElementGroup
      ergo_real maxabs = 0;
      for(int ii = groupIndexList[i]; ii < groupIndexList[i] + groupCountList[i]; ii++)
	{
	  ergo_real absval = std::fabs(distrList[ii].dmatElement);
	  if(absval > maxabs)
	    maxabs = absval;
	}
      groupList[i].maxAbsDmatElementGroup = maxabs;
    } // END FOR i

#define MAX_NO_OF_GROUPS_PER_CLUSTER 10


  // create clusters and chunks.
  // move groups into new list, one cluster at a time.
  int chunkCount = 0;
  int clusterCount = 0;
  int basisFuncPairCount = 0;

  std::vector<distr_group_struct> groupList2(groupCount);

  std::vector<cluster_struct> clusterList(groupCount);
  std::vector<chunk_struct> chunkList(groupCount);

  std::vector<basis_func_pair_struct> basisFuncPairList(distrCount);
  
  int noOfGroupsInNewList = 0;
  int noOfGroupsLeftInOldList = groupCount;
  while(noOfGroupsInNewList < groupCount)
    {
      // the group that is first now will define the beginning of a new cluster, and a new chunk.
      chunk_struct newChunk;
      memset(&newChunk, 0, sizeof(chunk_struct));
      clusterList[clusterCount].groupStartIndex = noOfGroupsInNewList;
      newChunk.clusterStartIndex = clusterCount;
      newChunk.basisFuncPairListIndex = basisFuncPairCount;

      // add basisFuncPairs for first group to newChunk
      for(int i = groupList[0].startIndex; i < groupList[0].startIndex + groupList[0].distrCount; i++)
	{
	  int alreadyInList = 0;
	  for(int kk = 0; kk < newChunk.noOfBasisFuncPairs; kk++)
	    {
	      if(distrList[i].basisFuncIndex_1 == basisFuncPairList[newChunk.basisFuncPairListIndex+kk].index_1 &&
		 distrList[i].basisFuncIndex_2 == basisFuncPairList[newChunk.basisFuncPairListIndex+kk].index_2)
		{
		  alreadyInList = 1;
		  break;
		}
	    } // END FOR kk
	  if(alreadyInList == 0)
	    {
	      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].index_1 = distrList[i].basisFuncIndex_1;
	      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].index_2 = distrList[i].basisFuncIndex_2;
	      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].pairIndex = distrList[i].pairIndex;
	      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].dmatElement = distrList[i].dmatElement;
	      newChunk.noOfBasisFuncPairs++;
	      basisFuncPairCount++;
	      if(newChunk.noOfBasisFuncPairs >= MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (newChunk.noOfBasisFuncPairs >= MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK)");
		  return -1;
		}
	    }
	} // END FOR i add basisFuncPairs for first group to newChunk

      int noOfClustersInCurrChunk = 1;
      int oldListIndex = 0;
      memcpy(&groupList2[noOfGroupsInNewList], &groupList[0], sizeof(distr_group_struct));
      noOfGroupsInNewList++;
      int noOfGroupsInCurrCluster = 1;
      // now find other groups with same exponent and same nmax
      ergo_real exponent = groupList[0].exponent;
      int nmax = groupList[0].nmax;
      for(int i = 1; i < noOfGroupsLeftInOldList; i++)
	{
	  ergo_real absexponentDiff = std::fabs(exponent - groupList[i].exponent);
	  if(absexponentDiff < 1e-11 && groupList[i].nmax == nmax && noOfGroupsInCurrCluster < MAX_NO_OF_GROUPS_PER_CLUSTER)
	    {
	      // same exponent and nmax found, add this group to cluster
	      memcpy(&groupList2[noOfGroupsInNewList], &groupList[i], sizeof(distr_group_struct));
	      noOfGroupsInNewList++;
	      noOfGroupsInCurrCluster++;
	      // add basisFuncPairs for group to newChunk
	      for(int ii = groupList[i].startIndex; ii < groupList[i].startIndex + groupList[i].distrCount; ii++)
		{
		  int alreadyInList = 0;
		  for(int kk = 0; kk < newChunk.noOfBasisFuncPairs; kk++)
		    {
		      if(distrList[ii].basisFuncIndex_1 == basisFuncPairList[newChunk.basisFuncPairListIndex+kk].index_1 &&
			 distrList[ii].basisFuncIndex_2 == basisFuncPairList[newChunk.basisFuncPairListIndex+kk].index_2)
			{
			  alreadyInList = 1;
			  break;
			}
		    }
		  if(alreadyInList == 0)
		    {
		      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].index_1 = distrList[ii].basisFuncIndex_1;
		      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].index_2 = distrList[ii].basisFuncIndex_2;
		      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].pairIndex = distrList[ii].pairIndex;
		      basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].dmatElement = distrList[ii].dmatElement;
		      newChunk.noOfBasisFuncPairs++;
		      basisFuncPairCount++;
		      if(newChunk.noOfBasisFuncPairs >= MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK)
			{
			  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (newChunk.noOfBasisFuncPairs >= MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK)");
			  return -1;
			}
		    }
		}
	    }
	  else
	    {
	      memcpy(&groupList[oldListIndex], &groupList[i], sizeof(distr_group_struct));
	      oldListIndex++;
	    }
	} // END FOR i
      noOfGroupsLeftInOldList -= noOfGroupsInCurrCluster;
      clusterList[clusterCount].noOfGroups = noOfGroupsInCurrCluster;
      clusterCount++;
      // the cluster just created is the first one in a new chunk.
      // if possible, we want to add more clusters for that chunk.
      int definingClusterStartGrIndex = clusterList[clusterCount-1].groupStartIndex;
      int definingClusterGrCount = clusterList[clusterCount-1].noOfGroups;
      // look for other clusters to put in the same chunk.
      noOfGroupsInCurrCluster = 0;
      while(noOfGroupsInNewList < groupCount)// && noOfGroupsInCurrCluster < MAX_NO_OF_GROUPS_PER_CLUSTER)
	{
	  // look for a group that has the right basis funcs.
	  int foundIndex = -1;
	  for(int i = 0; i < noOfGroupsLeftInOldList; i++)
	    {
	      // we demand that all basisfuncpairs must be present in the chunk (defined by first cluster)
	      int allPresentSoFar = 1;
	      for(int ii = 0; ii < groupList[i].distrCount; ii++)
		{
		  // check if this distr is present in the chunk
		  int bfidx1 = distrList[groupList[i].startIndex+ii].basisFuncIndex_1;
		  int bfidx2 = distrList[groupList[i].startIndex+ii].basisFuncIndex_2;
		  int found = 0;
		  for(int gr = definingClusterStartGrIndex; gr < definingClusterStartGrIndex + definingClusterGrCount; gr++)
		    {
		      int idistr;
		      for(idistr = 0; idistr < groupList2[gr].distrCount; idistr++)
			{
			  if(distrList[groupList2[gr].startIndex+idistr].basisFuncIndex_1 == bfidx1 && distrList[groupList2[gr].startIndex+idistr].basisFuncIndex_2 == bfidx2)
			    {
			      found = 1;
			      break;
			    }
			}
		      if(found == 1)
			break;
		    }
		  if(found == 0)
		    {
		      allPresentSoFar = 0;
		      break;
		    }
		} // END FOR ii
	      if(allPresentSoFar == 1)
		{
		  // OK, use this group
		  foundIndex = i;
		  break;
		}
	    } // END FOR i look for a group that has the right basis funcs.
	  if(foundIndex == -1)
	    break;
	  // OK, we have a group with accepted basis funcs.
	  // This group will be the first in a new cluster.
	  
	  clusterList[clusterCount].groupStartIndex = noOfGroupsInNewList;
	  int oldListIndex = 0;
	  memcpy(&groupList2[noOfGroupsInNewList], &groupList[foundIndex], sizeof(distr_group_struct));
	  noOfGroupsInNewList++;
	  noOfGroupsInCurrCluster = 1;

	  // add basisFuncPairs for group to newChunk
	  for(int ii = groupList[foundIndex].startIndex; ii < groupList[foundIndex].startIndex + groupList[foundIndex].distrCount; ii++)
	    {
	      int alreadyInList = 0;
	      for(int kk = 0; kk < newChunk.noOfBasisFuncPairs; kk++)
		{
		  if(distrList[ii].basisFuncIndex_1 == basisFuncPairList[newChunk.basisFuncPairListIndex+kk].index_1 &&
		     distrList[ii].basisFuncIndex_2 == basisFuncPairList[newChunk.basisFuncPairListIndex+kk].index_2)
		    {
		      alreadyInList = 1;
		      break;
		    }
		}
	      if(alreadyInList == 0)
		{
		  basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].index_1 = distrList[ii].basisFuncIndex_1;
		  basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].index_2 = distrList[ii].basisFuncIndex_2;
		  basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].pairIndex = distrList[ii].pairIndex;
		  basisFuncPairList[newChunk.basisFuncPairListIndex+newChunk.noOfBasisFuncPairs].dmatElement = distrList[ii].dmatElement;
		  newChunk.noOfBasisFuncPairs++;
		  basisFuncPairCount++;
		  if(newChunk.noOfBasisFuncPairs >= MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK)
		    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (newChunk.noOfBasisFuncPairs >= MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK)");
		      return -1;
		    }
		}
	    }

	  ergo_real exponent = groupList[foundIndex].exponent;
	  int nmax = groupList[foundIndex].nmax;
	  
	  // we have copied the entry at foundIndex to new list, all after that must be moved one step.
	  for(int i = foundIndex+1; i < noOfGroupsLeftInOldList; i++)
	    memcpy(&groupList[i-1], &groupList[i], sizeof(distr_group_struct));
	  noOfGroupsLeftInOldList--;

	  int noOfGroupsInCurrCluster = 1;
	  // now find other groups with same exponent and same nmax and accepted basis funcs
	  oldListIndex = 0;
	  for(int i = 0; i < noOfGroupsLeftInOldList; i++)
	    {
	      int addToCluster = 0;
	      ergo_real absexponentDiff = std::fabs(exponent - groupList[i].exponent);
	      if(absexponentDiff < 1e-11 && groupList[i].nmax == nmax && noOfGroupsInCurrCluster < MAX_NO_OF_GROUPS_PER_CLUSTER)
		{
		  // same exponent and nmax found, now check basis funcs
		  int allPresentSoFar = 1;
		  for(int ii = 0; ii < groupList[i].distrCount; ii++)
		    {
		      // check if this distr is present in the chunk
		      int bfidx1 = distrList[groupList[i].startIndex+ii].basisFuncIndex_1;
		      int bfidx2 = distrList[groupList[i].startIndex+ii].basisFuncIndex_2;
		      int found = 0;
		      for(int gr = definingClusterStartGrIndex; gr < definingClusterStartGrIndex + definingClusterGrCount; gr++)
			{
			  for(int idistr = 0; idistr < groupList2[gr].distrCount; idistr++)
			    {
			      if(distrList[groupList2[gr].startIndex+idistr].basisFuncIndex_1 == bfidx1 && 
				 distrList[groupList2[gr].startIndex+idistr].basisFuncIndex_2 == bfidx2)
				{
				  found = 1;
				  break;
				}
			    }
			  if(found == 1)
			    break;
			}
		      if(found == 0)
			{
			  allPresentSoFar = 0;
			  break;
			}
		    } // END FOR ii
		  if(allPresentSoFar == 1)
		    addToCluster = 1;
		}
	      if(addToCluster == 1)
		{
		  // same exponent and nmax found and accepted funcs, add this group to cluster
		  memcpy(&groupList2[noOfGroupsInNewList], &groupList[i], sizeof(distr_group_struct));
		  noOfGroupsInNewList++;
		  noOfGroupsInCurrCluster++;
		}
	      else
		{
		  if(i != oldListIndex)
		    memcpy(&groupList[oldListIndex], &groupList[i], sizeof(distr_group_struct));
		  oldListIndex++;
		}
	    } // END FOR i
	  noOfGroupsLeftInOldList -= noOfGroupsInCurrCluster-1;
	  clusterList[clusterCount].noOfGroups = noOfGroupsInCurrCluster;
	  clusterCount++;
	  noOfClustersInCurrChunk++;
	} // END WHILE look for other clusters to put in the same chunk
      
      newChunk.noOfClusters = noOfClustersInCurrChunk;
      chunkList[chunkCount] = newChunk;
      chunkCount++;
      
    } // END WHILE create clusters

  // check all chunks
  for(int i = 0; i < chunkCount; i++)
    {
      for(int j = 0; j < chunkList[i].noOfBasisFuncPairs; j++)
	{
	  for(int k = 0; k < chunkList[i].noOfBasisFuncPairs; k++)
	    {
	      if(j != k &&
		 basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_1 == basisFuncPairList[chunkList[i].basisFuncPairListIndex+k].index_1 &&
		 basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_2 == basisFuncPairList[chunkList[i].basisFuncPairListIndex+k].index_2)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: basisFuncPairs not unique in chunk");
		  return -1;
		}
	    }
	}
    }



  memcpy(&groupList[0], &groupList2[0], groupCount*sizeof(distr_group_struct));


  // OK, clusters and chunks done.
  


  // set nmax and exponent for all clusters
  for(int i = 0; i < clusterCount; i++)
    {
      int groupStartIndex = clusterList[i].groupStartIndex;
      int nGroups = clusterList[i].noOfGroups;
      int nmax = 0;
      ergo_real exponent = groupList[groupStartIndex].exponent;
      for(int j = groupStartIndex; j < groupStartIndex + nGroups; j++)
	{
	  if(groupList[j].nmax > nmax)
	    nmax = groupList[j].nmax;
	  ergo_real exponentdiff = std::fabs(groupList[j].exponent - exponent);
	  if(exponentdiff > 1e-11)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: different exponents found in same cluster");
	      return -1;
	    }
	} // END FOR j
      clusterList[i].nmax = nmax;
      clusterList[i].exponent = exponent;
    } // END FOR i set nmax for all clusters




  // Sort clusters according to chunks
  std::vector<cluster_struct> tempClusterList(clusterCount);
  int count = 0;
  for(int i = 0; i < chunkCount; i++)
    {
      int savedCount = count;
      for(int j = chunkList[i].clusterStartIndex; j < chunkList[i].clusterStartIndex + chunkList[i].noOfClusters; j++)
	{
	  tempClusterList[count] = clusterList[j];
	  count++;
	} // END FOR j
      chunkList[i].clusterStartIndex = savedCount;
    } // END FOR i
  memcpy(&clusterList[0], &tempClusterList[0], clusterCount*sizeof(cluster_struct));
  tempClusterList.clear();

  // Sort groups according to clusters, and set maxLimitingFactorForCluster
  std::vector<distr_group_struct> tempGroupList(groupCount);
  count = 0;
  for(int i = 0; i < clusterCount; i++)
    {
      ergo_real maxLimitingFactorForCluster = 0;
      int savedCount = count;
      for(int j = clusterList[i].groupStartIndex; j < clusterList[i].groupStartIndex + clusterList[i].noOfGroups; j++)
	{
	  ergo_real maxLimitingFactor = groupList[j].maxLimitingFactorGroup;
	  if(maxLimitingFactor > maxLimitingFactorForCluster)
	    maxLimitingFactorForCluster = maxLimitingFactor;
	  tempGroupList[count] = groupList[j];
	  count++;
	} // END FOR j
      clusterList[i].groupStartIndex = savedCount;
      clusterList[i].maxLimitingFactorForCluster = maxLimitingFactorForCluster;
    } // END FOR i
  memcpy(&groupList[0], &tempGroupList[0], groupCount*sizeof(distr_group_struct));
  tempGroupList.clear();

  // Sort distrs according to groups
  std::vector<DistributionSpecStructLabeled> tempDistrList(distrCount);
  //output_current_memory_usage("organize_distributions after allocating tempDistrList");
  count = 0;
  for(int i = 0; i < groupCount; i++)
    {
      int savedCount = count;
      for(int j = groupList[i].startIndex; j < groupList[i].startIndex + groupList[i].distrCount; j++)
	{
	  tempDistrList[count] = distrList[j];
	  count++;
	} // END FOR j
      groupList[i].startIndex = savedCount;
    } // END FOR i
  memcpy(&distrList[0], &tempDistrList[0], distrCount*sizeof(DistributionSpecStructLabeled));
  tempDistrList.clear();
  
  
  result->minimalDistrList.resize(distrCount);
  minimal_distr_struct* minimalDistrList = &result->minimalDistrList[0];
  for(int i = 0; i < distrCount; i++)
    {
      minimalDistrList[i].coeff = distrList[i].distr.coeff;
      minimalDistrList[i].monomialIndex = integralInfo.monomial_info.monomial_index_list
	[(int)distrList[i].distr.monomialInts[0]]
	[(int)distrList[i].distr.monomialInts[1]]
	[(int)distrList[i].distr.monomialInts[2]];
    }

  
  // get maxExtent
  ergo_real maxExtent = 0;
  for(int i = 0; i < distrCount; i++)
    {
      if(distrList[i].distr.extent > maxExtent)
	maxExtent = distrList[i].distr.extent;
    }
  result->maxExtent = maxExtent;

  // get maxDistanceOutsideBox
  ergo_real maxDistanceOutsideBox = 0;
  for(int i = 0; i < distrCount; i++)
    {
      // get minWallDist : minimum wall distance
      ergo_real minWallDist = boxWidth;
      int coordIndex;
      for(coordIndex = 0; coordIndex< 3; coordIndex++)
	{
	  // get wall distance for this coordinate
	  ergo_real dx = distrList[i].distr.centerCoords[coordIndex] - boxCenterCoords[coordIndex];
	  ergo_real wallDist = boxWidth / 2 - std::fabs(dx);
	  if(wallDist < minWallDist)
	    minWallDist = wallDist;
	} // END FOR coordIndex
      if(minWallDist < -0.00001)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (minWallDist < -0.00001)");
	  return -1;
	}
      ergo_real distanceOutsideBox = distrList[i].distr.extent - minWallDist;
      if(distanceOutsideBox > maxDistanceOutsideBox)
	maxDistanceOutsideBox = distanceOutsideBox;
    }
  result->maxDistanceOutsideBox = maxDistanceOutsideBox;


  for(int i = 0; i < chunkCount; i++)
    for(int j = chunkList[i].clusterStartIndex; j < chunkList[i].clusterStartIndex + chunkList[i].noOfClusters; j++)
      {
	int k_start = clusterList[j].groupStartIndex;
	int k_end = k_start + clusterList[j].noOfGroups;
	for(int k = k_start; k < k_end; k++)
	  {
	    int m_start = groupList[k].startIndex;
	    int m_end = m_start + groupList[k].distrCount;
	    for(int m = m_start; m < m_end; m++)
	      {
		int foundIndex = -1;
		for(int kk = 0; kk < chunkList[i].noOfBasisFuncPairs; kk++)
		  {
		    if(basisFuncPairList[chunkList[i].basisFuncPairListIndex+kk].index_1 == distrList[m].basisFuncIndex_1 &&
		       basisFuncPairList[chunkList[i].basisFuncPairListIndex+kk].index_2 == distrList[m].basisFuncIndex_2)
		      {
			foundIndex = kk;
			break;
		      }
		  } // END FOR kk
		if(foundIndex < 0)
		  {
		    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error setting basisFuncPairIndex");
		    return -1;
		  }
		minimalDistrList[m].basisFuncPairIndex = foundIndex;
	      }
	  }
      }


  // within each group, sort minimalDistrList by monomialIndex and basisFuncPairIndex, 
  // join distrs that differ only in coefficient.
  for(int i = 0; i < groupCount; i++)
    {
      minimal_distr_struct* p = & minimalDistrList[groupList[i].startIndex];
      int count = groupList[i].distrCount;
      for(int j = 0; j < count; j++)
	for(int k = 0; k < count - 1 - j; k++)
	  {
	    int doSwitch = 0;
	    if(p[k].monomialIndex > p[k+1].monomialIndex)
	      doSwitch = 1;
	    if(p[k].monomialIndex == p[k+1].monomialIndex)
	      {
		if(p[k].basisFuncPairIndex > p[k+1].basisFuncPairIndex)
		  doSwitch = 1;
	      }
	    if(doSwitch == 1)
	      {
		minimal_distr_struct temp;
		temp = p[k];
		p[k] = p[k+1];
		p[k+1] = temp;
	      }
	  } // END FOR j k
      // OK, list sorted.
      // We want to join together any entries that differ only in coefficient.
      int j = 0;
      int ii = 0;
      while(ii < count)
	{
	  ergo_real coeffSum = p[ii].coeff;
	  int k = ii + 1;
	  while(k < count)
	    {
	      if(p[k].monomialIndex != p[ii].monomialIndex || p[k].basisFuncPairIndex != p[ii].basisFuncPairIndex)
		break;
	      coeffSum += p[k].coeff;
	      k++;
	    }
	  p[j] = p[ii];
	  p[j].coeff = coeffSum;
	  j++;
	  int nResult = k - ii;
	  ii += nResult;
	}
      groupList[i].distrCount = j;
    } // END FOR i
  // Now go through groups again to move the distrs together now that the groups are smaller.
  count = 0;
  for(int i = 0; i < groupCount; i++)
    {
      int oldStartIndex = groupList[i].startIndex;
      groupList[i].startIndex = count;
      int distrCount = groupList[i].distrCount;
      for(int j = 0; j < distrCount; j++)
	{
	  minimalDistrList[count] = minimalDistrList[oldStartIndex+j];
	  count++;
	}
    } // END FOR i
  // check that no group contains repeating distrs
  for(int i = 0; i < groupCount; i++)
    {
      minimal_distr_struct* p = & minimalDistrList[groupList[i].startIndex];
      int distrCount = groupList[i].distrCount;
      for(int j = 0; j < distrCount; j++)
	for(int k = j+1; k < distrCount; k++)
	  {
	    if(p[k].monomialIndex == p[j].monomialIndex && p[k].basisFuncPairIndex == p[j].basisFuncPairIndex)
	      {
		do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: identical distrs found in same group.");
		return -1;
	      }
	  }
    }


  int basisFuncForChunksCount = 0;
  // Now get list of basis func indeces occurring in each chunk, store in basisFuncListForChunks.
  std::vector<int> basisFuncListForChunks(2*distrCount);

  for(int i = 0; i < chunkCount; i++)
    {

      int count = 0;
      for(int j = 0; j < chunkList[i].noOfBasisFuncPairs; j++)
	{
	  int i1 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_1;
	  int i2 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_2;
	  // Check if i1 and i2 are already present.
	  int i1_found = 0;
	  int i2_found = 0;
	  for(int k = 0; k < count; k++)
	    {
	      if(basisFuncListForChunks[basisFuncForChunksCount+k] == i1)
		i1_found = 1;
	      if(basisFuncListForChunks[basisFuncForChunksCount+k] == i2)
		i2_found = 1;
	    } // END FOR k
	  if(i1_found == 0)
	    {
	      basisFuncListForChunks[basisFuncForChunksCount+count] = i1;
	      count++;
	    }
	  if(i2_found == 0 && i1 != i2)
	    {
	      basisFuncListForChunks[basisFuncForChunksCount+count] = i2;
	      count++;
	    }
	} // END FOR j
	  // sort list for this chunk
      do_sort_int_list(&basisFuncListForChunks[basisFuncForChunksCount], count);
      // now "rename" index_1 and index_2 using basisFuncListForChunks.
      for(int j = 0; j < chunkList[i].noOfBasisFuncPairs; j++)
	{
	  int i1 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_1;
	  int i2 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_2;
	  // find positions of i1 and i2..
	  int i1_index = -1;
	  int i2_index = -1;
	  for(int k = 0; k < count; k++)
	    {
	      if(basisFuncListForChunks[basisFuncForChunksCount+k] == i1)
		i1_index = k;
	      if(basisFuncListForChunks[basisFuncForChunksCount+k] == i2)
		i2_index = k;
	    } // END FOR k
	  if(i1_index < 0 || i2_index < 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: error 1!!!");
	      return -1;
	    }
	  basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_1_mod = i1_index;
	  basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_2_mod = i2_index;
	} // END FOR j

      chunkList[i].basisFuncForChunksIndex = basisFuncForChunksCount;
      chunkList[i].basisFuncForChunkCount = count;
      basisFuncForChunksCount += count;
    } // END FOR i

  // Check result
  for(int i = 0; i < chunkCount; i++)
    {
      for(int j = 0; j < chunkList[i].noOfBasisFuncPairs; j++)
	{
	  int i1 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_1;
	  int i2 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_2;
	  // Check if i1 and i2 are present.
	  int i1_found = 0;
	  int i2_found = 0;
	  for(int k = 0; k < chunkList[i].basisFuncForChunkCount; k++)
	    {
	      if(basisFuncListForChunks[chunkList[i].basisFuncForChunksIndex+k] == i1)
		i1_found = 1;
	      if(basisFuncListForChunks[chunkList[i].basisFuncForChunksIndex+k] == i2)
		i2_found = 1;
	    } // END FOR k
	  if(i1_found == 0 || i2_found == 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: error !!");
	      return -1;
	    }
	} // END FOR j
    } // END FOR i



  int basisFuncListCount = 0;
  // Now get list of basis func indices occurring, store in basisFuncList.
  // Use basisFuncListForChunks to do this.
  std::vector<int> basisFuncList(basisFuncForChunksCount);
  memcpy(&basisFuncList[0], &basisFuncListForChunks[0], basisFuncForChunksCount*sizeof(int));
  std::sort(&basisFuncList[0], &basisFuncList[basisFuncForChunksCount]);
      
  int prevIndex = -1;
  int i = 0;
  while(i < basisFuncForChunksCount)
    {
      // now i points to a new basis func index.
      // check that sort order is OK.
      if(basisFuncList[i] < prevIndex)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error 11! i = %i, basisFuncList[i] = %i, prevIndex = %i", i, basisFuncList[i], prevIndex);
	  return -1;
	}
      basisFuncList[basisFuncListCount] = basisFuncList[i];
      basisFuncListCount++;
      prevIndex = basisFuncList[i];
      do i++; while(i < basisFuncForChunksCount &&
		    basisFuncList[i] == prevIndex);
    }

  // Now go through chunks again to "rename" indices according to basisFuncList.
  for(int i = 0; i < chunkCount; i++)
    {
      for(int j = 0; j < chunkList[i].noOfBasisFuncPairs; j++)
	{
	  int i1 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_1;
	  int i2 = basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_2;
	  // find positions of i1 and i2..
	  int i1_index = -1;
	  int i2_index = -1;
	  for(int k = 0; k < basisFuncListCount; k++)
	    {
	      if(basisFuncList[k] == i1)
		i1_index = k;
	      if(basisFuncList[k] == i2)
		i2_index = k;
	    } // END FOR k
	  if(i1_index < 0 || i2_index < 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error 3!!!");
	      return -1;
	    }
	  basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_inbox_1 = i1_index;
	  basisFuncPairList[chunkList[i].basisFuncPairListIndex+j].index_inbox_2 = i2_index;
	} // END FOR j

    } // END FOR i
  
  // take care of basisFuncListForChunks_map
  result->basisFuncListForChunks_map.resize(basisFuncForChunksCount);
  for(int i = 0; i < basisFuncForChunksCount; i++)
    {
      for(int k = 0; k < basisFuncListCount; k++)
	{
	  if(basisFuncListForChunks[i] == basisFuncList[k])
	    result->basisFuncListForChunks_map[i] = k;
	}
    }


  // take care of integral conversion matrices
  int noOfBasisFuncPairs_max = 0;
  for(int i = 0; i < chunkCount; i++) {
    int noOfBasisFuncPairs = chunkList[i].noOfBasisFuncPairs;
    if(noOfBasisFuncPairs > noOfBasisFuncPairs_max)
      noOfBasisFuncPairs_max = noOfBasisFuncPairs;
  }
  int noOfMonomials_max = 0;
  for(int i = 0; i < chunkCount; i++)
    for(int j = chunkList[i].clusterStartIndex; j < chunkList[i].clusterStartIndex + chunkList[i].noOfClusters; j++) {
      int noOfMonomials = integralInfo.monomial_info.no_of_monomials_list[clusterList[j].nmax];
      if(noOfMonomials > noOfMonomials_max)
	noOfMonomials_max = noOfMonomials;
    }
  // We do not know the final size of the spMatElementList yet. We give it some size to start with, and then increase the size when needed.
  std::vector<i_j_val_struct> spMatElementList(5*noOfBasisFuncPairs_max*noOfMonomials_max);
  int spMatElementListCount = 0;
  result->spMatCountList.resize(groupCount);
  result->spMatIdxList.resize(groupCount);
  for(int i = 0; i < chunkCount; i++)
    for(int j = chunkList[i].clusterStartIndex; j < chunkList[i].clusterStartIndex + chunkList[i].noOfClusters; j++) {
      // Now we are dealing with cluster j
      int noOfBasisFuncPairs = chunkList[i].noOfBasisFuncPairs;
      int nmax = clusterList[j].nmax;
      int noOfMonomials = integralInfo.monomial_info.no_of_monomials_list[nmax];
      int group_k_start = clusterList[j].groupStartIndex;
      int group_k_end = group_k_start + clusterList[j].noOfGroups;
      // Get conversion matrices
      ergo_real alpha = groupList[group_k_start].exponent;
      i_j_val_struct convMat_sp[noOfMonomials*noOfMonomials];
      int convMat_nnz = integralInfo.get_hermite_conversion_matrix_right_sparse(nmax, 1.0/alpha, convMat_sp);
      for(int group_k = group_k_start; group_k < group_k_end; group_k++) {
	i_j_val_struct BB1_x_Ai1_x_convMat1_sp[noOfBasisFuncPairs*noOfMonomials];
	int BB1_x_Ai1_x_convMat1_nnz = 0;
	get_conversion_matrix_for_group(integralInfo,
					groupList[group_k],
					nmax,
					minimalDistrList,
					noOfBasisFuncPairs, 
					convMat_sp,
					convMat_nnz,
					BB1_x_Ai1_x_convMat1_sp, // result
					BB1_x_Ai1_x_convMat1_nnz);  // result
	spmat_sort_elements(BB1_x_Ai1_x_convMat1_sp, BB1_x_Ai1_x_convMat1_nnz);
	// Check if the size of spMatElementList needs to be extended.
	while((int)(spMatElementList.size()) < spMatElementListCount + BB1_x_Ai1_x_convMat1_nnz)
	  spMatElementList.resize(2*spMatElementList.size());
	memcpy(&spMatElementList[spMatElementListCount], BB1_x_Ai1_x_convMat1_sp, BB1_x_Ai1_x_convMat1_nnz*sizeof(i_j_val_struct));
	result->spMatCountList[group_k] = BB1_x_Ai1_x_convMat1_nnz;
	result->spMatIdxList[group_k] = spMatElementListCount;
	spMatElementListCount += BB1_x_Ai1_x_convMat1_nnz;
	if(spMatElementListCount > groupCount*noOfBasisFuncPairs_max*noOfMonomials_max)
	  return -1;
      }
    }



  result->spMatElementList.resize(spMatElementListCount);
  memcpy(&result->spMatElementList[0], &spMatElementList[0], spMatElementListCount*sizeof(i_j_val_struct));

  result->clusterList.resize(clusterCount);
  memcpy(&result->clusterList[0], &clusterList[0], clusterCount*sizeof(cluster_struct));

  result->chunkList.resize(chunkCount);
  memcpy(&result->chunkList[0], &chunkList[0], chunkCount*sizeof(chunk_struct));

  result->basisFuncPairList.resize(basisFuncPairCount);
  memcpy(&result->basisFuncPairList[0], &basisFuncPairList[0], basisFuncPairCount*sizeof(basis_func_pair_struct));

  result->basisFuncListForChunks.resize(basisFuncForChunksCount);
  memcpy(&result->basisFuncListForChunks[0], &basisFuncListForChunks[0], basisFuncForChunksCount*sizeof(int));

  result->basisFuncList.resize(basisFuncListCount);
  memcpy(&result->basisFuncList[0], &basisFuncList[0], basisFuncListCount*sizeof(int));

  result->chunkCount = chunkCount;
  result->clusterCount = clusterCount;
  result->minimalDistrCount = distrCount;
  result->groupCount = groupCount;
  result->basisFuncPairCount = basisFuncPairCount;
  result->basisFuncForChunksCount = basisFuncForChunksCount;
  result->basisFuncListCount = basisFuncListCount;

  memcpy(&distrList_in[0], &distrList[0], distrCount*sizeof(DistributionSpecStructLabeled));

  return 0;
} // END organize_distributions

