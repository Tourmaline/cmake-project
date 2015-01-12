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


#include <string.h>
#include <stdio.h>

#include "integrals_2el_exchange.h"
#include "integrals_2el_utils.h"
#include "integrals_hermite.h"
#include "mm_limit_table.h"
#include "pi.h"
#include "pthread.h"
#include "utilities.h"
#include "matrix_algebra.h"
#include "integrals_2el_util_funcs.h"


static const int HUGE_INTEGER_NUMBER = 2000000000;


typedef struct
{
  int boxIndex_1;
  int boxIndex_2;
  int useMultipole;
  ergo_real distance;
} job_list_entry_K_struct;



static inline int
ergo_csr_find_index_inline(const csr_matrix_struct* csr, int row, int col)
{
  int n = csr->rowList[row].noOfElementsInRow;
  int baseIndex = csr->rowList[row].firstElementIndex;
  int* colList = &csr->columnIndexList[baseIndex];
  int lo = 0;
  int hi = n-1;
  while(lo < hi - 1)
    {
      int mid = (lo + hi) / 2;
      if(colList[mid] < col)
	lo = mid;
      else
	hi = mid;
    }
  if(colList[lo] == col)
    return baseIndex + lo;
  if(colList[hi] == col)
    return baseIndex + hi;
  // Not found
  return -1;
}



static inline ergo_real 
ergo_CSR_get_element_inline(const csr_matrix_struct* csr, 
			    int row,
			    int col)
{
  int row2 = row;
  int col2 = col;
  if(csr->symmetryFlag)
    {
      if(row > col)
	{
	  row2 = col;
	  col2 = row;
	}
    }
  int i = ergo_csr_find_index_inline(csr, row2, col2);
  if(i < 0)
    return 0;
  return csr->elementList[i];
}


typedef struct
{
  int a, b, c, d;
  int poly_ab_index;
  int poly_cd_index;
  int idx1;
  int idx2;
  ergo_real densValue;
} abcd_struct;

#define set_abcd_list_item_macro(i,A,B,C,D,v,i1,i2)			\
  list[i].a = A; list[i].b = B; list[i].c = C; list[i].d = D; list[i].densValue = v; list[i].idx1 = i1; list[i].idx2 = i2; 


pthread_mutex_t K_CSR_shared_access_mutex = PTHREAD_MUTEX_INITIALIZER;

static int 
get_K_contribs_from_2_interacting_boxes(const BasisInfoStruct & basisInfo,
					const IntegralInfo & integralInfo,
					const JK::ExchWeights & CAM_params,
					int maxNoOfMonomials,
					ergo_real* K,
					csr_matrix_struct* K_CSR_shared,
					const ergo_real* dens,
					const csr_matrix_struct* dens_CSR,
					int symmetryFlag,
					const distr_org_struct & distr_org_struct_1,
					const distr_org_struct & distr_org_struct_2,
					int interactionWithSelf,
					ergo_real threshold,
					JK_contribs_buffer_struct* bufferStructPtr,
					int use_multipole_screening_for_clusters,
					ergo_real boxDistance)
{
  int n = basisInfo.noOfBasisFuncs;

  const ergo_real twoTimesPiToPow5half = 2 * pitopow52;// = 2 * pow(pi, 2.5);
  ergo_real* summedIntegralList = bufferStructPtr->summedIntegralList;
  ergo_real* primitiveIntegralList = bufferStructPtr->primitiveIntegralList;

  int nChunks_1 = distr_org_struct_1.chunkCount;
  int nChunks_2 = distr_org_struct_2.chunkCount;
  const chunk_struct* chunkList_1 = &distr_org_struct_1.chunkList[0];
  const chunk_struct* chunkList_2 = &distr_org_struct_2.chunkList[0];
  const cluster_struct* clusterList_1 = &distr_org_struct_1.clusterList[0];
  const cluster_struct* clusterList_2 = &distr_org_struct_2.clusterList[0];
  const distr_group_struct* groupList_1 = &distr_org_struct_1.groupList[0];
  const distr_group_struct* groupList_2 = &distr_org_struct_2.groupList[0];
  const basis_func_pair_struct* basisFuncPairList_1 = &distr_org_struct_1.basisFuncPairList[0];
  const basis_func_pair_struct* basisFuncPairList_2 = &distr_org_struct_2.basisFuncPairList[0];
  const int* basisFuncListForChunks_1 = &distr_org_struct_1.basisFuncListForChunks[0];
  const int* basisFuncListForChunks_map_1 = &distr_org_struct_1.basisFuncListForChunks_map[0];
  const int* basisFuncListForChunks_2 = &distr_org_struct_2.basisFuncListForChunks[0];
  const int* basisFuncListForChunks_map_2 = &distr_org_struct_2.basisFuncListForChunks_map[0];
  const int* basisFuncList_1 = &distr_org_struct_1.basisFuncList[0];
  int basisFuncList_1_count = distr_org_struct_1.basisFuncListCount;
  const int* basisFuncList_2 = &distr_org_struct_2.basisFuncList[0];
  int basisFuncList_2_count = distr_org_struct_2.basisFuncListCount;

  const i_j_val_struct* spMatElementList_1 = &distr_org_struct_1.spMatElementList[0];
  const int* spMatCountList_1 = &distr_org_struct_1.spMatCountList[0];
  const int* spMatIdxList_1 = &distr_org_struct_1.spMatIdxList[0];
  const i_j_val_struct* spMatElementList_2 = &distr_org_struct_2.spMatElementList[0];
  const int* spMatCountList_2 = &distr_org_struct_2.spMatCountList[0];
  const int* spMatIdxList_2 = &distr_org_struct_2.spMatIdxList[0];

  // Set up "partial box-box density matrix"
  int nnn1 = basisFuncList_1_count;
  int nnn2 = basisFuncList_2_count;

  ergo_real* partial_dmat_1 = bufferStructPtr->partial_dmat_1;
  ergo_real* partial_dmat_2 = bufferStructPtr->partial_dmat_2;

  ergo_real maxabsdmatelement_boxbox = 0;
  for(int i1 = 0; i1 < nnn1; i1++)
    for(int i2 = 0; i2 < nnn2; i2++) {
      int a = basisFuncList_1[i1];
      int b = basisFuncList_2[i2];
      if(dens)
	partial_dmat_1[i1*nnn2+i2] = dens[a*n+b];
      else
	partial_dmat_1[i1*nnn2+i2] = ergo_CSR_get_element_inline(dens_CSR, a, b);
      ergo_real absval = std::fabs(partial_dmat_1[i1*nnn2+i2]);
      if(absval > maxabsdmatelement_boxbox)
	maxabsdmatelement_boxbox = absval;
    }

  if(symmetryFlag == 0) {
    for(int i1 = 0; i1 < nnn1; i1++)
      for(int i2 = 0; i2 < nnn2; i2++) {
	int a = basisFuncList_1[i1];
	int b = basisFuncList_2[i2];
	if(dens)
	  partial_dmat_2[i1*nnn2+i2] = dens[b*n+a];
	else
	  partial_dmat_2[i1*nnn2+i2] = ergo_CSR_get_element_inline(dens_CSR, b, a);
	ergo_real absval = std::fabs(partial_dmat_2[i1*nnn2+i2]);
	if(absval > maxabsdmatelement_boxbox)
	  maxabsdmatelement_boxbox = absval;
      }
  }

  ergo_real* partial_K_1 = bufferStructPtr->partial_K_1;
  ergo_real* partial_K_2 = bufferStructPtr->partial_K_2;

  for(int i1 = 0; i1 < nnn1; i1++)
    for(int i2 = 0; i2 < nnn2; i2++) {
      partial_K_1[i1*nnn2+i2] = 0;
      if(symmetryFlag == 0)
	partial_K_2[i1*nnn2+i2] = 0;
    }

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


	  // Set up "local density matrix" for this pair of chunks.
	  int nn1 = chunkList_1[chunk_i].basisFuncForChunkCount;
	  int nn2 = chunkList_2[chunk_j].basisFuncForChunkCount;
	  ergo_real local_dmat_1[nn1][nn2];
	  ergo_real local_dmat_2[nn1][nn2];
	  ergo_real maxabsdmatelement = 0;
	  for(int i1 = 0; i1 < nn1; i1++)
	    for(int i2 = 0; i2 < nn2; i2++)
	      {
		int a =   basisFuncListForChunks_1[chunkList_1[chunk_i].basisFuncForChunksIndex+i1];
		int b =   basisFuncListForChunks_2[chunkList_2[chunk_j].basisFuncForChunksIndex+i2];
		int a2 =  basisFuncListForChunks_map_1[chunkList_1[chunk_i].basisFuncForChunksIndex+i1];
		int b2 =  basisFuncListForChunks_map_2[chunkList_2[chunk_j].basisFuncForChunksIndex+i2];
		if(dens) {
		  local_dmat_1[i1][i2] = dens[a*n+b];
		  if(symmetryFlag == 0)
		    local_dmat_2[i1][i2] = dens[b*n+a];
		}
		else {
		  local_dmat_1[i1][i2] = partial_dmat_1[a2*nnn2+b2];
		  if(symmetryFlag == 0)
		    local_dmat_2[i1][i2] = partial_dmat_2[a2*nnn2+b2];
		}
		ergo_real absval = std::fabs(local_dmat_1[i1][i2]);
		if(absval > maxabsdmatelement)
		  maxabsdmatelement = absval;
		if(symmetryFlag == 0) {
		  ergo_real absval = std::fabs(local_dmat_2[i1][i2]);
		  if(absval > maxabsdmatelement)
		    maxabsdmatelement = absval;
		}
	      }

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

		  if(use_multipole_screening_for_clusters == 1)
		    {
		      // Try multipole screening
		      int maxDegree = 2;
		      ergo_real maxAbsContributionFromMultipole = mm_limits_get_max_abs_mm_contrib(maxDegree,
												   clusterList_1[cluster_i].multipoleEuclideanNormList,
												   maxDegree,
												   clusterList_2[cluster_j].multipoleEuclideanNormList,
												   boxDistance);
		      if(maxAbsContributionFromMultipole * maxabsdmatelement < threshold)
			continue;
		    } // END IF try multipole screening
		  
		  int group_i_start = clusterList_1[cluster_i].groupStartIndex;
		  int group_i_end = group_i_start + clusterList_1[cluster_i].noOfGroups;
		  int group_j_start = clusterList_2[cluster_j].groupStartIndex;
		  int group_j_end = group_j_start + clusterList_2[cluster_j].noOfGroups;

		  int n1max = clusterList_1[cluster_i].nmax;
		  int n2max = clusterList_2[cluster_j].nmax;
		  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1max];
		  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2max];

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
			  // Try Cauchy-Schwartz screening
			  if(groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabsdmatelement < threshold)
			    continue;

			  ergo_real dx = groupList_2[group_j].centerCoords[0] - groupList_1[group_i].centerCoords[0];
			  ergo_real dy = groupList_2[group_j].centerCoords[1] - groupList_1[group_i].centerCoords[1];
			  ergo_real dz = groupList_2[group_j].centerCoords[2] - groupList_1[group_i].centerCoords[2];

			  // Check if multipole screening can be used
			  ergo_real distance = std::sqrt(dx*dx+dy*dy+dz*dz);
			  if(distance > groupList_1[group_i].maxExtentGroup + groupList_2[group_j].maxExtentGroup)
			    {
			      // Try multipole screening
			      int maxDegree = 2;
			      ergo_real maxAbsContributionFromMultipole = mm_limits_get_max_abs_mm_contrib(maxDegree,
													   groupList_1[group_i].multipoleEuclideanNormList,
													   maxDegree,
													   groupList_2[group_j].multipoleEuclideanNormList,
													   distance);
			      if(maxAbsContributionFromMultipole * maxabsdmatelement < threshold)
				continue;
			    } // END IF try multipole screening

			  // now we can do all integrals needed for this pair of groups
			  // now we have dx dy dz alpha0 alpha1 n1max n2max. Get all integrals for this case.
			  get_related_integrals_hermite(integralInfo,
							CAM_params,
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

		int a =     basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_1;
		int b =     basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_2;
		int c =     basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_1;
		int d =     basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_2;
		ergo_real integralValueCurr = summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_2];

		int a_mod = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_1_mod;
		int b_mod = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_2_mod;
		int c_mod = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_1_mod;
		int d_mod = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_2_mod;

		int a_mod2 = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_inbox_1;
		int b_mod2 = basisFuncPairList_1[chunkList_1[chunk_i].basisFuncPairListIndex+idx_1].index_inbox_2;
		int c_mod2 = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_inbox_1;
		int d_mod2 = basisFuncPairList_2[chunkList_2[chunk_j].basisFuncPairListIndex+idx_2].index_inbox_2;


		if(a == c && b == d)
		  integralValueCurr *= 2;

		if(std::fabs(integralValueCurr)*maxabsdmatelement < threshold)
		  continue;

		ergo_real dens_ac = local_dmat_1[a_mod][c_mod];
		ergo_real dens_ad = local_dmat_1[a_mod][d_mod];
		ergo_real dens_bc = local_dmat_1[b_mod][c_mod];
		ergo_real dens_bd = local_dmat_1[b_mod][d_mod];

		ergo_real dens_ca = dens_ac;
		ergo_real dens_da = dens_ad;
		ergo_real dens_cb = dens_bc;
		ergo_real dens_db = dens_bd;
		
		if(symmetryFlag == 0) {
		  dens_ca = local_dmat_2[a_mod][c_mod];
		  dens_da = local_dmat_2[a_mod][d_mod];
		  dens_cb = local_dmat_2[b_mod][c_mod];
		  dens_db = local_dmat_2[b_mod][d_mod];
		}





		if(symmetryFlag && !K) {
		  
		  if(a != b && c != d && a != c && a != d && b != c && b != d) {
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
		  }
		  else if(a == b && c != d && a != c && a != d && b != c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
		  }
		  else if(a != b && c == d && a != c && a != d && b != c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
		  }
		  else if(a != b && c != d && a == c && a != d && b != c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr * 2.0;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
		  }
		  else if(a != b && c != d && a != c && a == d && b != c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr * 2.0;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
		  }
		  else if(a != b && c != d && a != c && a != d && b == c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr * 2.0;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
		  }
		  else if(a != b && c != d && a != c && a != d && b != c && b == d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr * 2.0;
		  }
		  else if(a != b && c != d && a == c && b == d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
		  }
		  else if(a == b && c == d && a != c && a != d && b != c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		  }
		  else if(a == b && c == d && a == c && a == d && b == c && b == d) { // OK
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		  }
		  else if(a == b && c != d && a == c && a != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr*2.0;
		  }
		  else if(a == b && c != d && a != c && a == d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr*2.0;
		    partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
		  }

		  else if(a != b && c == d && a == c && b != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr*2.0;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
		  }
		  else if(a != b && c == d && b == c && a != d) { // OK
		    partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
		    partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr*2.0;
		  }
		  else {
		    return -1;
		  }



		}
		else if(a != b && c != d && a != c && a != d && b != c && b != d)
		  {
		    if(symmetryFlag)
		      {
			partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
			partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
			partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
			partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;

			if(K)
			  {
			    if(d >= a)
			      K[a*n+d] += -0.5 * dens_bc * integralValueCurr;
			    else
			      K[d*n+a] += -0.5 * dens_bc * integralValueCurr;
			    if(c >= a)
			      K[a*n+c] += -0.5 * dens_bd * integralValueCurr;
			    else
			      K[c*n+a] += -0.5 * dens_bd * integralValueCurr;
			    if(c >= b)
			      K[b*n+c] += -0.5 * dens_ad * integralValueCurr;
			    else
			      K[c*n+b] += -0.5 * dens_ad * integralValueCurr;
			    if(d >= b)
			      K[b*n+d] += -0.5 * dens_ac * integralValueCurr;
			    else
			      K[d*n+b] += -0.5 * dens_ac * integralValueCurr;
			  }
		      }
		    else
		      {
			partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
			partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
			partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
			partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;

			partial_K_2[a_mod2*nnn2+d_mod2] += -0.5 * dens_cb * integralValueCurr;
			partial_K_2[a_mod2*nnn2+c_mod2] += -0.5 * dens_db * integralValueCurr;
			partial_K_2[b_mod2*nnn2+c_mod2] += -0.5 * dens_da * integralValueCurr;
			partial_K_2[b_mod2*nnn2+d_mod2] += -0.5 * dens_ca * integralValueCurr;

			if(K)
			  {
			    K[a*n+d] += -0.5 * dens[b*n+c] * integralValueCurr;
			    K[a*n+c] += -0.5 * dens[b*n+d] * integralValueCurr;
			    K[b*n+d] += -0.5 * dens[a*n+c] * integralValueCurr;
			    K[b*n+c] += -0.5 * dens[a*n+d] * integralValueCurr;
			    K[c*n+b] += -0.5 * dens[d*n+a] * integralValueCurr;
			    K[d*n+b] += -0.5 * dens[c*n+a] * integralValueCurr;
			    K[c*n+a] += -0.5 * dens[d*n+b] * integralValueCurr;
			    K[d*n+a] += -0.5 * dens[c*n+b] * integralValueCurr;
			  }
		      }
		  }
		else
		  {		    
		    abcd_struct list[8];
		    
		    /* determine unique configurations */
		    set_abcd_list_item_macro(0, a, b, c, d, dens_bc, a_mod2, d_mod2);
		    set_abcd_list_item_macro(1, a, b, d, c, dens_bd, a_mod2, c_mod2);
		    set_abcd_list_item_macro(2, b, a, c, d, dens_ac, b_mod2, d_mod2);
		    set_abcd_list_item_macro(3, b, a, d, c, dens_ad, b_mod2, c_mod2);

		    set_abcd_list_item_macro(4, c, d, a, b, dens_da, b_mod2, c_mod2);
		    set_abcd_list_item_macro(5, d, c, a, b, dens_ca, b_mod2, d_mod2);
		    set_abcd_list_item_macro(6, c, d, b, a, dens_db, a_mod2, c_mod2);
		    set_abcd_list_item_macro(7, d, c, b, a, dens_cb, a_mod2, d_mod2);

		    int ccc = 0;
	  
		    for(int ii = 0; ii < 8; ii++)
		      {
			abcd_struct* abcd = &list[ii];
			int aa, bb, cc, dd;

			/* check if this is a new unique configuration */
			int unique = 1;
			for(int jj = 0; jj < ii; jj++)
			  {
			    if(abcd->a == list[jj].a && 
			       abcd->b == list[jj].b && 
			       abcd->c == list[jj].c && 
			       abcd->d == list[jj].d)
			      unique = 0;
			  }
			if(unique == 0)
			  continue;
			/* now we know that this configuration is unique. */
			aa = abcd->a;
			bb = abcd->b;
			cc = abcd->c;
			dd = abcd->d;

			ccc++;

			if(symmetryFlag)
			  {
			    if(dd >= aa)
			      {
				partial_K_1[abcd->idx1*nnn2+abcd->idx2] += -0.5 * abcd->densValue * integralValueCurr;
				if(K)
				  K[aa*n+dd] += -0.5 * abcd->densValue * integralValueCurr;
			      }
			  }
			else
			  {
			    if(ii <= 3)
			      partial_K_1[abcd->idx1*nnn2+abcd->idx2] += -0.5 * abcd->densValue * integralValueCurr;
			    else
			      partial_K_2[abcd->idx1*nnn2+abcd->idx2] += -0.5 * abcd->densValue * integralValueCurr;
			    if(K)
			      K[aa*n+dd] += -0.5 * dens[bb*n+cc] * integralValueCurr;
			  }
			
		      } /* END FOR ii go through 8 configurations */
		  }

	      } // END FOR idx_1 idx_2
	} // END FOR chunk_j
    } // END FOR chunk_i

  if(K_CSR_shared)
    {
      if(K_CSR_shared->n)
	{
	  pthread_mutex_lock(&K_CSR_shared_access_mutex);
	  // Now move results from partial_K to K.
	  for(int i1 = 0; i1 < nnn1; i1++)
	    for(int i2 = 0; i2 < nnn2; i2++)
	      {
		int a = basisFuncList_1[i1];
		int b = basisFuncList_2[i2];
		ergo_CSR_add_to_element(K_CSR_shared, 
					a,
					b,
					partial_K_1[i1*nnn2+i2]);
		if(symmetryFlag == 0) {
		  ergo_CSR_add_to_element(K_CSR_shared, 
					  b,
					  a,
					  partial_K_2[i1*nnn2+i2]);
		}
	      }
	  pthread_mutex_unlock(&K_CSR_shared_access_mutex);
	}
    }

  return 0;
}





static int
create_joblist_exchange_for_two_boxes_recursive(const BasisInfoStruct & basisInfo,
						const IntegralInfo & integralInfo,
						int maxNoOfMonomials,
						ergo_real threshold,
						const box_struct* boxList,
						int numberOfLevels,
						const csr_matrix_struct* dmatLimitMatrixCSRList,
						const int* basisFuncGroupCounterList,
						int currLevel,
						int boxIndex_1,
						int boxIndex_2,
						job_list_entry_K_struct* jobList_K,
						int maxNoOfJobs
						)
{
  // Check if this pair of boxes can be skipped.
  int noOfRelevantBasisFuncGroups_1 = boxList[boxIndex_1].noOfRelevantBasisFuncGroups;
  int noOfRelevantBasisFuncGroups_2 = boxList[boxIndex_2].noOfRelevantBasisFuncGroups;
  const csr_matrix_struct* dmatLimitMatrixCSR = &dmatLimitMatrixCSRList[currLevel];

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
  
  ergo_real maxDistanceOutsideBox_1 = boxList[boxIndex_1].distrListForK.org.maxDistanceOutsideBox;
  ergo_real maxDistanceOutsideBox_2 = boxList[boxIndex_2].distrListForK.org.maxDistanceOutsideBox;
  
  int useMultipole = 0;
  if(boxIndex_1 != boxIndex_2 && distance > maxDistanceOutsideBox_1 + maxDistanceOutsideBox_2)
    useMultipole = 1;
  
  ergo_real maxValue_CauschySchwartz = 0;
  ergo_real maxValue_multipole = 0;
  for(int i = 0; i < noOfRelevantBasisFuncGroups_1; i++)
    for(int j = 0; j < noOfRelevantBasisFuncGroups_2; j++)
      {
	ergo_real size_1 = boxList[boxIndex_1].basisFuncGroupInfoList[i].max_CS_factor;
	int index_1 = boxList[boxIndex_1].basisFuncGroupInfoList[i].basisFuncGroupIndex;
	ergo_real size_2 = boxList[boxIndex_2].basisFuncGroupInfoList[j].max_CS_factor;
	int index_2 = boxList[boxIndex_2].basisFuncGroupInfoList[j].basisFuncGroupIndex;
	if(index_1 < 0 || index_2 < 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive: (index_1 < 0 || index_2 < 0)");
	    return -1;
	  }
	ergo_real maxDensElement = ergo_CSR_get_element(dmatLimitMatrixCSR, index_1, index_2);
	ergo_real currMax = size_1 * size_2 * maxDensElement;
	if(currMax > maxValue_CauschySchwartz)
	  maxValue_CauschySchwartz = currMax;
	if(useMultipole == 1)
	  {
	    int degreeNeeded_1 = boxList[boxIndex_1].basisFuncGroupInfoList[i].maxMultipoleDegree;
	    int degreeNeeded_2 = boxList[boxIndex_2].basisFuncGroupInfoList[j].maxMultipoleDegree;
	    ergo_real maxAbsContributionFromMultipole = mm_limits_get_max_abs_mm_contrib(degreeNeeded_1,
											 boxList[boxIndex_1].basisFuncGroupInfoList[i].maxMomentVectorNormList,
											 degreeNeeded_2,
											 boxList[boxIndex_2].basisFuncGroupInfoList[j].maxMomentVectorNormList,
											 distance);
	    ergo_real currMaxFromMultipole = maxAbsContributionFromMultipole * maxDensElement;
	    if(currMaxFromMultipole > maxValue_multipole)
	      maxValue_multipole = currMaxFromMultipole;
	  } // END IF useMultipole
      } // END FOR i j
  
  if(useMultipole == 1 && maxValue_multipole < threshold)
    return 0;
  if(maxValue_CauschySchwartz < threshold)
    return 0;
  if(currLevel == numberOfLevels-1)
    {
      // We are at the level of smallest boxes. Add job to job list.
      if(maxNoOfJobs <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive: (maxNoOfJobs <= 0)");
	  return -1;
	}
      if(jobList_K != NULL)
	{
	  jobList_K[0].boxIndex_1 = boxIndex_1;
	  jobList_K[0].boxIndex_2 = boxIndex_2;
	  jobList_K[0].useMultipole = useMultipole;
	  jobList_K[0].distance = distance;
	}
      return 1;
    }
  // Go to next level. Do interaction between all pairs of children of the two boxes.
  int noOfChildren_1 = boxList[boxIndex_1].basicBox.noOfChildBoxes;
  int noOfChildren_2 = boxList[boxIndex_2].basicBox.noOfChildBoxes;
  int jobCount = 0;
  for(int i = 0; i < noOfChildren_1; i++)
    {
      int start_j = 0;
      if(boxIndex_1 == boxIndex_2)
	start_j = i;
      for(int j = start_j; j < noOfChildren_2; j++)
	{
	  int childIndex_1 = boxList[boxIndex_1].basicBox.firstChildBoxIndex + i;
	  int childIndex_2 = boxList[boxIndex_2].basicBox.firstChildBoxIndex + j;
	  job_list_entry_K_struct* jobList_K_mod = NULL;
	  if(jobList_K != NULL)
	    jobList_K_mod = &jobList_K[jobCount];
	  int noOfJobs = create_joblist_exchange_for_two_boxes_recursive(basisInfo,
									 integralInfo,
									 maxNoOfMonomials,
									 threshold,
									 boxList,
									 numberOfLevels,
									 dmatLimitMatrixCSRList,
									 basisFuncGroupCounterList,
									 currLevel + 1,
									 childIndex_1,
									 childIndex_2,
									 jobList_K_mod,
									 maxNoOfJobs - jobCount
									 );
	  if(noOfJobs < 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive for child boxes");
	      return -1;
	    }
	  jobCount += noOfJobs;
	} // END FOR j
    } // END FOR i
  return jobCount;
}






struct K_joblist_thread_struct
{
  pthread_t thread;
  const BasisInfoStruct & basisInfo;
  const IntegralInfo* integralInfo;
  const JK::ExchWeights & CAM_params;
  ergo_real* K;
  csr_matrix_struct* K_CSR_shared;
  const ergo_real* dens;
  const csr_matrix_struct* densCSR;
  int maxNoOfMonomials;
  int basisFuncListCount_max;
  ergo_real threshold;
  const box_struct* boxList;
  const job_list_entry_K_struct* jobList_K;
  int noOfJobs_K_total;
  int thread_ID;
  int noOfThreads;
  int resultCode;
  int symmetryFlag;
  K_joblist_thread_struct(const BasisInfoStruct & basisInfoIn,
			  const JK::ExchWeights & CAM_paramsIn) :
    basisInfo(basisInfoIn), CAM_params(CAM_paramsIn) { }
};


static void*
execute_joblist_K_thread_func(void* arg)
{
  K_joblist_thread_struct* params = (K_joblist_thread_struct*)arg;
  try {
  const box_struct* boxList = params->boxList;
  ergo_real* K = params->K;
  int threadID = params->thread_ID;
  int noOfThreads = params->noOfThreads;

  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(*params->integralInfo, params->maxNoOfMonomials, params->basisFuncListCount_max, &bufferStruct);

  for(int jobIndex = 0; jobIndex < params->noOfJobs_K_total; jobIndex++)
    {
      if(jobIndex % noOfThreads != threadID)
	continue;
      int self = 0;
      int boxIndex_1 = params->jobList_K[jobIndex].boxIndex_1;
      int boxIndex_2 = params->jobList_K[jobIndex].boxIndex_2;
      if(boxIndex_1 == boxIndex_2)
	self = 1;
      if(get_K_contribs_from_2_interacting_boxes(params->basisInfo,
						 *params->integralInfo,
						 params->CAM_params,
						 params->maxNoOfMonomials,
						 K,
						 params->K_CSR_shared,
						 params->dens,
						 params->densCSR,
						 params->symmetryFlag,						 
						 boxList[boxIndex_1].distrListForK.org,
						 boxList[boxIndex_2].distrListForK.org,
						 self,
						 params->threshold,
						 &bufferStruct,
						 params->jobList_K[jobIndex].useMultipole,
						 params->jobList_K[jobIndex].distance
						 ) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_K_contribs_from_2_interacting_boxes");
	  params->resultCode = -1;
	  return NULL;
	}
    } // END FOR jobIndex

  if(params->symmetryFlag && K)
    {
      // Fill the other triangle of K
      int n = params->basisInfo.noOfBasisFuncs;
      for(int i = 0; i < n; i++)
	for(int j = 0; j < i; j++)
	  K[i*n+j] = K[j*n+i];
    }

  free_buffers_needed_by_integral_code(&bufferStruct);

  params->resultCode = 0;
  }
  catch(char const* e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "char* exception caught in execute_joblist_K_thread_func: '%s'", e);    
    do_output_time(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
    params->resultCode = -1;
  }
  catch (std::exception & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "std::exception caught in execute_joblist_K_thread_func: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
  }
  return NULL;
}



static int
execute_joblist_K_threaded(int noOfThreads,
			   const ergo_real* dens,
			   csr_matrix_struct* densCSR,
			   const BasisInfoStruct & basisInfo,
			   const IntegralInfo & integralInfo,
			   const JK::ExchWeights & CAM_params,
			   int maxNoOfMonomials,
			   int basisFuncListCount_max,
			   const box_struct* boxList,
			   const job_list_entry_K_struct* jobList_K,
			   int noOfJobs_K,
			   ergo_real threshold,
			   ergo_real* K,
			   csr_matrix_struct* K_CSR,
			   int symmetryFlag
			   )
{
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "execute_joblist_K_threaded, noOfThreads = %2i, basisFuncListCount_max = %5i", noOfThreads, basisFuncListCount_max);

  if((!dens && K) || (dens && !K))
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in execute_joblist_K_threaded: CSR used for one of D, K but not for the other!?!.");
      return -1;
    }

  K_joblist_thread_struct* threadParamsList[noOfThreads];
  
  // Set common parameters for all threads
  for(int i = 0; i < noOfThreads; i++)
    {
      threadParamsList[i] = new K_joblist_thread_struct(basisInfo, CAM_params);
      threadParamsList[i]->dens = dens;
      threadParamsList[i]->densCSR = densCSR;
      threadParamsList[i]->integralInfo = &integralInfo;
      threadParamsList[i]->maxNoOfMonomials = maxNoOfMonomials;
      threadParamsList[i]->basisFuncListCount_max = basisFuncListCount_max;
      threadParamsList[i]->boxList = boxList;
      threadParamsList[i]->jobList_K = jobList_K;
      threadParamsList[i]->noOfJobs_K_total = noOfJobs_K;
      threadParamsList[i]->noOfThreads = noOfThreads;
      threadParamsList[i]->resultCode = -1; // initialize to error code
      threadParamsList[i]->threshold = threshold;
      threadParamsList[i]->symmetryFlag = symmetryFlag;
      threadParamsList[i]->K = NULL; // will be set later if needed
      threadParamsList[i]->K_CSR_shared = K_CSR;
    } // END FOR i
  
  // Set result K pointer for thread 0
  // For the full matrix case, thread 0 uses the original K pointer.
  threadParamsList[0]->K = K;

  // Set result K pointer for other threads
  int n = basisInfo.noOfBasisFuncs;
  for(int i = 1; i < noOfThreads; i++)
    {
      if(K)
	{
	  threadParamsList[i]->K = new ergo_real[n*n];
	  memset(threadParamsList[i]->K, 0, n*n*sizeof(ergo_real));
	}
      else
	{
	  // CSR case; do nothing. In this case all threads share the
	  // same K_CSR pointer, coordinating access using a mutex.
	}
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
			execute_joblist_K_thread_func, 
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
  time_t workStartTime;
  time(&workStartTime);

  /* wait for threads to finish */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(pthread_join(threadParamsList[i]->thread, NULL) != 0)
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
    } /* END FOR i */

  time_t workEndTime;
  time(&workEndTime);
  int secondsTaken = workEndTime - workStartTime;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "all %i threads have finished, took %8i wall s.", noOfThreads, secondsTaken);
  
  /* now all threads have finished, check for errors */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(threadParamsList[i]->resultCode != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_K_thread_func"
		    " for thread %i", i);
	  return -1;
	}
    } /* END FOR i */
  

  // add contributions from other threads
  for(int i = 1; i < noOfThreads; i++)
    {
      if(K)
       {
	  for(int j = 0; j < n*n; j++)
	    K[j] += threadParamsList[i]->K[j];
	}
      else
	{
	  // CSR case. Do nothing here since all threads have already
	  // placed their results in the same shared K_CSR matrix.
	}
    }
  
  // Free extra K buffers used by threads.
  // Note that this loop must start with 1, not 0.
  for(int i = 1; i < noOfThreads; i++)
    {
      if(K)
	{
	  delete [] threadParamsList[i]->K;
	}
      else
	{
	  // CSR case. Do nothing here.
	}
    }

  for(int i = 0; i < noOfThreads; i++)
    delete threadParamsList[i];

  return 0;
}


static int
execute_joblist_K_serial(const ergo_real* dens,
                        csr_matrix_struct* densCSR,
			 const BasisInfoStruct & basisInfo,
			 const IntegralInfo & integralInfo,
			 const JK::ExchWeights & CAM_params,
			 int maxNoOfMonomials,
			 int basisFuncListCount_max,
			 const box_struct* boxList,
			 const job_list_entry_K_struct* jobList_K,
			 int noOfJobs_K,
			 ergo_real threshold,
			 ergo_real* K,
			 csr_matrix_struct* K_CSR,
			 int symmetryFlag)
{
  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(integralInfo, maxNoOfMonomials, basisFuncListCount_max, &bufferStruct);
  for(int jobIndex = 0; jobIndex < noOfJobs_K; jobIndex++)
    {
      int self = 0;
      int boxIndex_1 = jobList_K[jobIndex].boxIndex_1;
      int boxIndex_2 = jobList_K[jobIndex].boxIndex_2;
      if(boxIndex_1 == boxIndex_2)
	self = 1;

      if(get_K_contribs_from_2_interacting_boxes(basisInfo,
						 integralInfo,
						 CAM_params,
						 maxNoOfMonomials,
						 K,
                                                 K_CSR,
						 dens,
						 densCSR,
						 symmetryFlag,
						 boxList[boxIndex_1].distrListForK.org,
						 boxList[boxIndex_2].distrListForK.org,						
						 self,
						 threshold,
						 &bufferStruct,
						 jobList_K[jobIndex].useMultipole,
						 jobList_K[jobIndex].distance
						 ) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_K_contribs_from_2_interacting_boxes");
	  return -1;
	}
    } // END FOR jobIndex

  if(symmetryFlag && K)
    {
      // Fill the other triangle of K
      int n = basisInfo.noOfBasisFuncs;
      for(int i = 0; i < n; i++)
	for(int j = 0; j < i; j++)
	  K[i*n+j] = K[j*n+i];
    }

  free_buffers_needed_by_integral_code(&bufferStruct);

  return 0;
}





typedef struct
{
  int i1;
  int i2;
} basisFuncGroupPairStruct;


static int
compare_basisFuncGroupPairs(const void* p1, const void* p2)
{
  basisFuncGroupPairStruct* pair_1 = (basisFuncGroupPairStruct*)p1;
  basisFuncGroupPairStruct* pair_2 = (basisFuncGroupPairStruct*)p2;
  if(pair_1->i1 > pair_2->i1)
    return 1;
  if(pair_1->i1 < pair_2->i1)
    return -1;
  if(pair_1->i2 > pair_2->i2)
    return 1;
  if(pair_1->i2 < pair_2->i2)
    return -1;
  return 0;
}






static int
get_basisFuncGroupInfoList_size(int distrCountTot,
				const DistributionSpecStructLabeled* distrList,
				int numberOfLevels,
				const int* levelStartIndexList,
				const int* levelCounterList,
				const box_struct* boxList,
				int** basisFuncGroupList)
{
  int basisFuncGroupInfoList_count = 0;
  std::vector<basisFuncGroupPairStruct> pairList(2*distrCountTot);
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
    {
      int pairCount = 0;
      for(int i = levelStartIndexList[levelNumber]; i < levelStartIndexList[levelNumber] + levelCounterList[levelNumber]; i++)
	{
	  // go through all distrs of this box, and update basisFuncGroupInfoList accordingly.
	  int distrStartIndex = boxList[i].basicBox.firstItemIndex;
	  int distrCountCurrBox = boxList[i].basicBox.noOfItems;
	  for(int j = distrStartIndex; j < distrStartIndex + distrCountCurrBox; j++)
	    {
	      const DistributionSpecStructLabeled & currDistr = distrList[j];
	      int basisFuncGroup_1 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_1];
	      int basisFuncGroup_2 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_2];
	      pairList[pairCount].i1 = i;
	      pairList[pairCount].i2 = basisFuncGroup_1;
	      pairCount++;
	      pairList[pairCount].i1 = i;
	      pairList[pairCount].i2 = basisFuncGroup_2;
	      pairCount++;
	    } // END FOR j
	} // END FOR i
      // sort pairList
      qsort(&pairList[0], pairCount, sizeof(basisFuncGroupPairStruct), compare_basisFuncGroupPairs);
      int nn = 0;
      int i = 0;
      while(i < pairCount)
	{
	  // now i should point to a new i1
	  int i1 = pairList[i].i1;
	  int j = i;
	  while(j < pairCount && pairList[j].i1 == i1)
	    {
	      nn++;
	      int i2 = pairList[j].i2;
	      // now skip until another i2 is found.
	      while(j < pairCount && pairList[j].i1 == i1 && pairList[j].i2 == i2)
		j++;
	    }
	  i = j;
	}
      basisFuncGroupInfoList_count += nn;
    } // END FOR levelNumber
  return basisFuncGroupInfoList_count;
}




typedef struct
{
  int i1;
  int i2;
  ergo_real x;
} dmatElementStruct;

static int 
compare_dmatElements(const void* p1, const void* p2)
{
  dmatElementStruct* e1 = (dmatElementStruct*)p1;
  dmatElementStruct* e2 = (dmatElementStruct*)p2;
  if(e1->i1 > e2->i1)
    return 1;
  if(e1->i1 < e2->i1)
    return -1;
  if(e1->i2 > e2->i2)
    return 1;
  if(e1->i2 < e2->i2)
    return -1;
  return 0;
}




static int
create_reduced_vector(int nvalues,
		      const std::vector<dmatElementStruct> & dmatElementList,
		      std::vector<dmatElementStruct> & resultVector) {
  resultVector.resize(nvalues);
  int i = 0;
  int nvalues2 = 0;
  int curr_i1 = dmatElementList[0].i1;
  int curr_i2 = dmatElementList[0].i2;
  ergo_real curr_maxAbs = std::fabs(dmatElementList[0].x);
  while(i < nvalues) {
    i++;
    int closeCurr = 0;
    if(i == nvalues)
      closeCurr = 1;
    else {
      // now we know it is safe to access element i
      int i1 = dmatElementList[i].i1;
      int i2 = dmatElementList[i].i2;
      if(i1 != curr_i1 || i2 != curr_i2)
	closeCurr = 1;
      else {
	// now we know this i is just a continuation of the current batch
	ergo_real absx = std::fabs(dmatElementList[i].x);
	if(absx > curr_maxAbs)
	  curr_maxAbs = absx;
      }
    }
    if(closeCurr) {
      resultVector[nvalues2].i1 = curr_i1;
      resultVector[nvalues2].i2 = curr_i2;
      resultVector[nvalues2].x  = curr_maxAbs;
      nvalues2++;
      if(i < nvalues) {
	// Now we know it is safe to access element i. Start new batch.
	curr_i1 = dmatElementList[i].i1;
	curr_i2 = dmatElementList[i].i2;
	curr_maxAbs = std::fabs(dmatElementList[i].x);
      }
    }
  }
  resultVector.resize(nvalues2);
  return nvalues2;
} /* End create_reduced_vector */



static int 
getDmatLimitMatrixCSRList(csr_matrix_struct* dmatLimitMatrixCSRList, 
			  int numberOfLevels,
			  const csr_matrix_struct* densCSR,
			  const int* const* basisFuncGroupList,
			  const int* basisFuncGroupCounterList)
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "getDmatLimitMatrixCSRList start.");
  memset(dmatLimitMatrixCSRList, 0, numberOfLevels*sizeof(csr_matrix_struct));
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++) {
    // Populate dmatElementList with info for this level, one row at a time.
    int nn = basisFuncGroupCounterList[levelNumber];
    std::vector< std::vector<dmatElementStruct> > dmatElementListList(densCSR->n);
    for(int dmatrow = 0; dmatrow < densCSR->n; dmatrow++) {
      int nValuesCurrRow = ergo_CSR_get_nvalues_singlerow(densCSR, dmatrow);
      if(nValuesCurrRow == 0)
	continue;
      std::vector<int> colind(nValuesCurrRow);
      std::vector<ergo_real> values(nValuesCurrRow);
      if(ergo_CSR_get_values_singlerow(densCSR,
				       dmatrow,
				       &colind[0], 
				       &values[0],
				       nValuesCurrRow) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_get_values_singlerow.");
	return -1;
      }
      std::vector<dmatElementStruct> dmatElementList(nValuesCurrRow);
      for(int i = 0; i < nValuesCurrRow; i++) {
	int grrow = basisFuncGroupList[levelNumber][dmatrow];
	int grcol = basisFuncGroupList[levelNumber][colind[i]];
	if(grrow < grcol) {
	  dmatElementList[i].i1 = grrow;
	  dmatElementList[i].i2 = grcol;
	}
	else {
	  dmatElementList[i].i1 = grcol;
	  dmatElementList[i].i2 = grrow;
	}
	dmatElementList[i].x  = values[i];
      }
      // sort list to gather equal i1 i2 pairs together
      qsort(&dmatElementList[0], nValuesCurrRow, sizeof(dmatElementStruct), compare_dmatElements);
      // Create reduced vector.
      std::vector<dmatElementStruct> reducedVector;
      int nReduced = create_reduced_vector(nValuesCurrRow, dmatElementList, reducedVector);
      // Store result for this row in dmatElementListList.
      dmatElementListList[dmatrow].resize(nReduced);
      for(int i = 0; i < nReduced; i++)
	dmatElementListList[dmatrow][i] = reducedVector[i];
    }
    // OK, all rows done. Now create a single long list of all rows.
    int nTot = 0;
    for(int row = 0; row < densCSR->n; row++)
      nTot += dmatElementListList[row].size();
    std::vector<dmatElementStruct> dmatElementList(nTot);
    int count = 0;
    for(int row = 0; row < densCSR->n; row++)
      for(int i = 0; i < (int)dmatElementListList[row].size(); i++)
	dmatElementList[count++] = dmatElementListList[row][i];
    if(count != nTot) {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in getDmatLimitMatrixCSRList: (count != nTot).");
      return -1;
    }
    // OK, single long list created.
    // sort list to gather equal i1 i2 pairs together
    qsort(&dmatElementList[0], nTot, sizeof(dmatElementStruct), compare_dmatElements);
    // Create reduced vector.
    std::vector<dmatElementStruct> reducedVector;
    int nReduced = create_reduced_vector(nTot, dmatElementList, reducedVector);
    // Create CSR matrix for this level.
    std::vector<int> rowind2(nReduced);
    std::vector<int> colind2(nReduced);
    std::vector<ergo_real> values2(nReduced);
    for(int i = 0; i < nReduced; i++) {
      rowind2[i] = reducedVector[i].i1;
      colind2[i] = reducedVector[i].i2;
      values2[i] = reducedVector[i].x;
    }
    // Create CSR
    csr_matrix_struct* currCSR = &dmatLimitMatrixCSRList[levelNumber];
    if(ergo_CSR_create(currCSR, 
		       1,
		       nn,
		       nReduced,
		       &rowind2[0],
		       &colind2[0]) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_create for dmatLimitMatrixCSRList.");
	return -1;
      }
    for(int i = 0; i < nReduced; i++) {
      ergo_CSR_add_to_element(currCSR, 
			      rowind2[i],
			      colind2[i],
			      values2[i]);
    }
  }
  return 0;
}






/*
NOTE: This function adds its result to K.
This means that if only K is wanted, it must be set to zero before calling this function.
*/
int
compute_K_by_boxes(const BasisInfoStruct & basisInfo,
		   const IntegralInfo & integralInfo,
		   const JK::ExchWeights & CAM_params_in,
		   const JK::Params& J_K_params,
		   ergo_real* K,
                   csr_matrix_struct* K_CSR,
		   const ergo_real* dens,
		   csr_matrix_struct* densCSR,
		   int symmetryFlag)
{
  Util::TimeMeter timeMeterTot;

  Util::TimeMeter timeMeterIntMatLimits;
  ergo_real maxDistance = getSafeMaxDistance(basisInfo);
  mm_limits_init(maxDistance);
  timeMeterIntMatLimits.print(LOG_AREA_INTEGRALS, "mm_limits_init");

  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "entering compute_K_by_boxes, no of basis funcs = %5i, threshold_K = %7.3g, exchange_box_size = %6.2f", 
	    n, (double)J_K_params.threshold_K, (double)J_K_params.exchange_box_size);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "beginning of compute_K_by_boxes");

  const JK::ExchWeights CAM_params(CAM_params_in);


  Util::TimeMeter timeMeterDistrList;

  ergo_real maxDensityMatrixElement;
  if(dens)
    maxDensityMatrixElement = get_max_abs_vector_element(n*n, dens);
  else
    maxDensityMatrixElement = ergo_CSR_get_max_abs_element(densCSR);


  // get largest limiting factor
  ergo_real maxLimitingFactor = 0;
  if(get_list_of_labeled_distrs_maxLimitingFactor(basisInfo,
						  integralInfo,
						  J_K_params.threshold_K,
						  &maxLimitingFactor,
						  maxDensityMatrixElement) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs_maxLimitingFactor");
      return -1;
    }

  // Get number of distributions
  int distrCountTot = get_list_of_labeled_distrs(basisInfo,
						 integralInfo,
						 J_K_params.threshold_K,
						 NULL,
						 0,
						 maxLimitingFactor,
						 NULL,
						 maxDensityMatrixElement);
  if(distrCountTot == 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_K_by_boxes: (distrCountTot == 0), skipping.");
      return 0;
    }
  if(distrCountTot <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_K_by_boxes: (distrCountTot <= 0)");
      return -1;
    }

  std::vector<DistributionSpecStructLabeled> distrList(distrCountTot);

  // create list of product primitives, with labels
  int distrCountTemp = get_list_of_labeled_distrs(basisInfo,
						  integralInfo,
						  J_K_params.threshold_K,
						  &distrList[0],
						  distrCountTot,
						  maxLimitingFactor,
						  NULL,
						  maxDensityMatrixElement);
  if(distrCountTemp != distrCountTot)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_K_by_boxes:(distrCountTemp != distrCountTot)");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating list of primitive distributions");


  // compute extent for all distrs
  Util::TimeMeter timeMeterComputeExtentForAllDistrs;
  compute_extent_for_list_of_distributions(distrCountTot, 
					   &distrList[0], 
					   J_K_params.threshold_K,
					   maxLimitingFactor,
					   maxDensityMatrixElement);
  timeMeterComputeExtentForAllDistrs.print(LOG_AREA_INTEGRALS, "Compute extent for all distrs");
  
  
  // get maximum number of monomials
  int maxNoOfMonomials = 0;
  for(int i = 0; i < distrCountTot; i++)
    {
      int degree = 0;
      for(int j = 0; j < 3; j++)
	degree += distrList[i].distr.monomialInts[j];
      int noOfMonomials = integralInfo.monomial_info.no_of_monomials_list[degree];
      if(noOfMonomials > maxNoOfMonomials)
	maxNoOfMonomials = noOfMonomials;
    } // END FOR ABcount

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Creating list of distributions done, distrCountTot = %9i", distrCountTot);
  timeMeterDistrList.print(LOG_AREA_INTEGRALS, "Creating list of distributions");




  //
  // This is where we start to worry about the box system
  //

  Util::TimeMeter timeMeterBoxes;

  BoxSystem boxSystem;
  if(create_box_system_and_reorder_distrs(distrCountTot,
					  &distrList[0],
					  J_K_params.exchange_box_size,
					  boxSystem) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system_and_reorder_distrs");
      return -1;
    }

  // Create new list of boxes (more advanced boxes this time)
  std::vector<box_struct> boxList(boxSystem.totNoOfBoxes);
  // TODO: need to clear contents of boxList here?  
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


  // Set up basisFuncGroups for all levels
  // Create another box system, this time with basis functions as items

  std::vector<box_item_struct> itemListBasisFuncs(n);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < 3; j++)
	itemListBasisFuncs[i].centerCoords[j] = basisInfo.basisFuncList[i].centerCoords[j];
      itemListBasisFuncs[i].originalIndex = i;
    } // END FOR i

  
  const ergo_real maxToplevelBoxSizeBasisFuncs = J_K_params.exchange_box_size;
  
  BoxSystem boxSystemBasisFuncs;

  if(boxSystemBasisFuncs.create_box_system(&itemListBasisFuncs[0],
					   n,
					   maxToplevelBoxSizeBasisFuncs * 0.5) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after creating second box system");

  if(boxSystemBasisFuncs.noOfLevels < boxSystem.noOfLevels)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (boxSystemBasisFuncs.noOfLevels < boxSystem.noOfLevels)");
      return -1;
    }
  int noOfLevelsBasisFuncs = boxSystemBasisFuncs.noOfLevels;
  int noOfLevelsDiff = noOfLevelsBasisFuncs - numberOfLevels;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "noOfLevelsBasisFuncs = %i, noOfLevelsDiff = %i", noOfLevelsBasisFuncs, noOfLevelsDiff);
  
  int* basisFuncGroupList[numberOfLevels];
  int basisFuncGroupCounterList[numberOfLevels];
  for(int i = 0; i < numberOfLevels; i++)
    basisFuncGroupList[i] = new int[n];

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating basisFuncGroupList");

  int maxNoOfBasisFuncGroupsPerLevel = 0;
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
    {
      // Set up basisFuncGroup list for this level.
      int noOfBoxesCurrLevel = boxSystemBasisFuncs.levelList[levelNumber+noOfLevelsDiff].noOfBoxes;
      int startIndex         = boxSystemBasisFuncs.levelList[levelNumber+noOfLevelsDiff].startIndexInBoxList;
      for(int i = startIndex; i < startIndex + noOfBoxesCurrLevel; i++)
	{
	  // assign basis funcs of this box to basisFuncGroup i
	  int firstItemIndex = boxSystemBasisFuncs.boxList[i].firstItemIndex;
	  for(int j = firstItemIndex; j < firstItemIndex + boxSystemBasisFuncs.boxList[i].noOfItems; j++)
	    {
	      int basisFuncIndex = itemListBasisFuncs[j].originalIndex;
	      basisFuncGroupList[levelNumber][basisFuncIndex] = i - startIndex;
	    } // END FOR j
	} // END FOR i
      basisFuncGroupCounterList[levelNumber] = noOfBoxesCurrLevel;
      if(noOfBoxesCurrLevel > maxNoOfBasisFuncGroupsPerLevel)
	maxNoOfBasisFuncGroupsPerLevel = noOfBoxesCurrLevel;
    } // END FOR levelNumber
  
  // OK, basisFuncGroups done.


  // OK, boxes created.


  timeMeterBoxes.print(LOG_AREA_INTEGRALS, "Creating boxes etc");



  Util::TimeMeter timeMeterGetMultipoleNormVectors;

  // Create list of multipole norm vectors, for later use.
  std::vector<ergo_real> multipoleNormVectorList((MAX_MULTIPOLE_DEGREE_BASIC+1)*distrCountTot);
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating multipoleNormVectorList");

  std::vector<int> multipoleDegreeList(distrCountTot);
  for(int j = 0; j < distrCountTot; j++)
    {
      DistributionSpecStructLabeled* currDistr = &distrList[j];
      ergo_real* multipoleNormVectorList_curr = &multipoleNormVectorList[j*(MAX_MULTIPOLE_DEGREE_BASIC+1)];
      multipole_struct_small multipole;
      if(compute_multipole_moments(integralInfo, &currDistr->distr, &multipole) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_multipole_moments");
	  return -1;
	}
      multipoleDegreeList[j] = multipole.degree;
      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	multipoleNormVectorList_curr[l] = 0;
      for(int l = 0; l <= multipole.degree; l++)
	{
	  int startIndex = l*l;
	  int endIndex = (l+1)*(l+1);
	  ergo_real sum = 0;
	  for(int A = startIndex; A < endIndex; A++)
	    sum += multipole.momentList[A]*multipole.momentList[A];
	  ergo_real subNorm = std::sqrt(sum);
	  multipoleNormVectorList_curr[l] = subNorm;
	}
    } // END FOR j
  
  timeMeterGetMultipoleNormVectors.print(LOG_AREA_INTEGRALS, "getting multipoleNormVectorList");


  // For each box at each level, store information about the largest size distr associated with each basisFuncGroup.

  Util::TimeMeter timeMeterGetLimitsAllLevels;


  // Predict size of basisFuncGroupInfoList
  int basisFuncGroupInfoList_count_predicted = get_basisFuncGroupInfoList_size(distrCountTot,
									       &distrList[0],
									       numberOfLevels,
									       levelStartIndexList,
									       levelCounterList,
									       &boxList[0],
									       basisFuncGroupList);
  if(basisFuncGroupInfoList_count_predicted <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basisFuncGroupInfoList_size");
      return -1;
    }

  std::vector<basis_func_group_info_for_box> basisFuncGroupInfoList(basisFuncGroupInfoList_count_predicted);
  for(int i = 0; i < basisFuncGroupInfoList_count_predicted; i++)
    {
      // TODO: OK to set to -1 here? It should be OK...?
      basisFuncGroupInfoList[i].basisFuncGroupIndex = -1;
      basisFuncGroupInfoList[i].max_CS_factor = 0;
      basisFuncGroupInfoList[i].maxMultipoleDegree = -1;
      memset(basisFuncGroupInfoList[i].maxMomentVectorNormList, 0, (MAX_MULTIPOLE_DEGREE_BASIC+1)*sizeof(ergo_real));
    }
 
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating basisFuncGroupInfoList");

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "starting loop to setup basisFuncGroupInfoList, with added index checks.");
  
  int basisFuncGroupInfoListTotCount = 0;
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
    {
      // Set up basisFuncGroupInfoList for each box at level.
      
      for(int boxIndex = levelStartIndexList[levelNumber]; boxIndex < levelStartIndexList[levelNumber] + levelCounterList[levelNumber]; boxIndex++)
	{
	  if(boxIndex < 0 || boxIndex >= boxSystem.totNoOfBoxes)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (boxIndex < 0 || boxIndex >= boxSystem.totNoOfBoxes)");
	      return -1;
	    }
	  box_struct & currBox = boxList[boxIndex];
	  // Set basisFuncGroupInfoList pointer for this box to point
	  // to current position in previously allocated list.
	  currBox.basisFuncGroupInfoList = &basisFuncGroupInfoList[basisFuncGroupInfoListTotCount];
	  // Define maxCount for this box, so we can check that we do
	  // not try to access anything outside what is allocated.
	  int maxCount = basisFuncGroupInfoList_count_predicted - basisFuncGroupInfoListTotCount;
	  // go through all distrs of this box, and update basisFuncGroupInfoList accordingly.
	  int count = 0;
	  int distrStartIndex = currBox.basicBox.firstItemIndex;
	  int distrCountCurrBox = currBox.basicBox.noOfItems;
	  for(int j = distrStartIndex; j < distrStartIndex + distrCountCurrBox; j++)
	    {
	      if(j < 0 || j >= distrCountTot)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (j < 0 || j >= distrCountTot)");
		  return -1;
		}
	      const DistributionSpecStructLabeled & currDistr = distrList[j];
	      const ergo_real* multipoleNormVectorList_curr = &multipoleNormVectorList[j*(MAX_MULTIPOLE_DEGREE_BASIC+1)];
	      int multipoleDegree_curr = multipoleDegreeList[j];
	      int basisFuncGroup_1 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_1];
	      int basisFuncGroup_2 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_2];
	      ergo_real CS_factor = currDistr.limitingFactor;

	      // check if basisFuncGroup_1 and/or basisFuncGroup_2 is already present
	      int foundIndex_1 = -1;
	      int foundIndex_2 = -1;
	      if(count > maxCount)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (count > maxCount)");
		  return -1;
		}
	      for(int k = 0; k < count; k++)
		{
		  if(currBox.basisFuncGroupInfoList[k].basisFuncGroupIndex == basisFuncGroup_1)
		    foundIndex_1 = k;
		  if(currBox.basisFuncGroupInfoList[k].basisFuncGroupIndex == basisFuncGroup_2)
		    foundIndex_2 = k;
		}
	      if(foundIndex_1 >= 0)
		{
		  // check if max_CS_factor needs updating
		  if(CS_factor > currBox.basisFuncGroupInfoList[foundIndex_1].max_CS_factor)
		    currBox.basisFuncGroupInfoList[foundIndex_1].max_CS_factor = CS_factor;
		  // modfy maxMomentVectorNormList if needed.
		  if(multipoleDegree_curr > currBox.basisFuncGroupInfoList[foundIndex_1].maxMultipoleDegree)
		    currBox.basisFuncGroupInfoList[foundIndex_1].maxMultipoleDegree = multipoleDegree_curr;
		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    {
		      if(multipoleNormVectorList_curr[l] > currBox.basisFuncGroupInfoList[foundIndex_1].maxMomentVectorNormList[l])
			currBox.basisFuncGroupInfoList[foundIndex_1].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
		    }
		}
	      else
		{
		  // add new entry for basisFuncGroup_1
		  if(count >= maxCount)
		    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (count >= maxCount)");
		      return -1;
		    }
		  currBox.basisFuncGroupInfoList[count].basisFuncGroupIndex = basisFuncGroup_1;
		  currBox.basisFuncGroupInfoList[count].max_CS_factor = CS_factor;
		  currBox.basisFuncGroupInfoList[count].maxMultipoleDegree = multipoleDegree_curr;
		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    currBox.basisFuncGroupInfoList[count].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
		  count++;
		}
	      if(basisFuncGroup_2 != basisFuncGroup_1)
		{
		  if(foundIndex_2 >= 0)
		    {
		      // check if maxSize needs updating
		      if(CS_factor > currBox.basisFuncGroupInfoList[foundIndex_2].max_CS_factor)
			currBox.basisFuncGroupInfoList[foundIndex_2].max_CS_factor = CS_factor;
		      // modfy maxMomentVectorNormList if needed.
		      if(multipoleDegree_curr > currBox.basisFuncGroupInfoList[foundIndex_2].maxMultipoleDegree)
			currBox.basisFuncGroupInfoList[foundIndex_2].maxMultipoleDegree = multipoleDegree_curr;
		      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
			{
			  if(multipoleNormVectorList_curr[l] > currBox.basisFuncGroupInfoList[foundIndex_2].maxMomentVectorNormList[l])
			    currBox.basisFuncGroupInfoList[foundIndex_2].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
			}
		    }
		  else
		    {
		      // add new entry for basisFuncGroup_2
		      if(count >= maxCount)
			{
			  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (count >= maxCount)");
			  return -1;
			}
		      currBox.basisFuncGroupInfoList[count].basisFuncGroupIndex = basisFuncGroup_2;
		      currBox.basisFuncGroupInfoList[count].max_CS_factor = CS_factor;
		      currBox.basisFuncGroupInfoList[count].maxMultipoleDegree = multipoleDegree_curr;
		      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
			currBox.basisFuncGroupInfoList[count].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
		      count++;
		    }
		}
	      // OK, distr j done
	    } // END FOR j
	  // we are done for this box
	  currBox.noOfRelevantBasisFuncGroups = count;
	  basisFuncGroupInfoListTotCount += count;
	} // END FOR boxIndex
    } // END FOR levelNumber

  if(basisFuncGroupInfoListTotCount != basisFuncGroupInfoList_count_predicted)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (basisFuncGroupInfoListTotCount != basisFuncGroupInfoList_count_predicted)");
      return -1;
    }

  // OK, basisFuncGroup info done for all boxes.

  timeMeterGetLimitsAllLevels.print(LOG_AREA_INTEGRALS, "GetLimitsAllLevels");







  int noOfBoxesTopLevel = levelCounterList[numberOfLevels-1];
  box_struct* boxListTopLevel = &boxList[levelStartIndexList[numberOfLevels-1]];

  int basisFuncListCount_max = 0;
  
  // Now call organize_distributions for each top-level box
  Util::TimeMeter timeMeterKorg;
  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      DistributionSpecStructLabeled* distrListCurrBox = &distrList[boxListTopLevel[i].basicBox.firstItemIndex];
      int distrCountCurrBox = boxListTopLevel[i].basicBox.noOfItems;
      if(organize_distributions(integralInfo,
				distrListCurrBox,
				distrCountCurrBox,
				&boxListTopLevel[i].distrListForK.org,
				boxListTopLevel[i].basicBox.centerCoords,
				boxListTopLevel[i].basicBox.width) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in organize_distributions for box %i", i);
	  return -1;
	}
      if(boxListTopLevel[i].distrListForK.org.basisFuncListCount > basisFuncListCount_max)
	basisFuncListCount_max = boxListTopLevel[i].distrListForK.org.basisFuncListCount;
    } // END FOR i
  timeMeterKorg.print(LOG_AREA_INTEGRALS, "K organize_distributions for all boxes");





  // Now go through the other levels, getting info for parent boxes

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
	  // We want to get maxDistanceOutsideBox for parent box (use largest value found among the children).
	  ergo_real maxDistanceOutsideBox = 0;	  
	  for(int childIndex = 0; childIndex < noOfChildren; childIndex++)
	    {
	      int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
	      box_struct* childBox = &boxList[childIndexInBoxList];	      
	      if(childBox->distrListForK.org.maxDistanceOutsideBox > maxDistanceOutsideBox)
		maxDistanceOutsideBox = childBox->distrListForK.org.maxDistanceOutsideBox;
	    } // END FOR childIndex
	  currBox->distrListForK.org.maxDistanceOutsideBox = maxDistanceOutsideBox;
	} // END FOR boxIndex
    } // END FOR levelNumber












  // Generate multipole limits for each group, and for each box.
  //// Also get maxExtentGroup for each group.
  
  Util::TimeMeter timeMeterGetMultipoleLimits;

  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      const chunk_struct* chunkList = &boxListTopLevel[i].distrListForK.org.chunkList[0];
      cluster_struct* clusterList = &boxListTopLevel[i].distrListForK.org.clusterList[0];
      distr_group_struct* groupList = &boxListTopLevel[i].distrListForK.org.groupList[0];
      const minimal_distr_struct* minimalDistrList = &boxListTopLevel[i].distrListForK.org.minimalDistrList[0];
      int chunkCount = boxListTopLevel[i].distrListForK.org.chunkCount;

      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	boxListTopLevel[i].multipoleEuclideanNormList[l] = 0;

      boxListTopLevel[i].largestCSfactor = 0;
      
      for(int chunkIndex = 0; chunkIndex < chunkCount; chunkIndex++)
	{
	  int clusterCount = chunkList[chunkIndex].noOfClusters;
	  int cluster_start = chunkList[chunkIndex].clusterStartIndex;
	  for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++)
	    {
	      int group_start = clusterList[clusterIndex].groupStartIndex;
	      int group_end = group_start + clusterList[clusterIndex].noOfGroups;
	      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		clusterList[clusterIndex].multipoleEuclideanNormList[l] = 0;
	      for(int groupIndex = group_start; groupIndex < group_end; groupIndex++)
		{
		  distr_group_struct* currGroup = &groupList[groupIndex];

		  ergo_real maxMomentVectorNormForDistrsList[MAX_MULTIPOLE_DEGREE_BASIC+1];
		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    maxMomentVectorNormForDistrsList[l] = 0;
		 
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

		      // modfy maxMomentVectorNormForDistrsList if needed.		      
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
		      
		    } // END FOR distrIndex

		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    currGroup->multipoleEuclideanNormList[l] = maxMomentVectorNormForDistrsList[l];

		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    {
		      if(currGroup->multipoleEuclideanNormList[l] > clusterList[clusterIndex].multipoleEuclideanNormList[l])
			clusterList[clusterIndex].multipoleEuclideanNormList[l] = currGroup->multipoleEuclideanNormList[l];
		    }

		  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
		    {
		      if(maxMomentVectorNormForDistrsList[l] > boxListTopLevel[i].multipoleEuclideanNormList[l])
			boxListTopLevel[i].multipoleEuclideanNormList[l] = maxMomentVectorNormForDistrsList[l];
		    }

		  if(currGroup->maxLimitingFactorGroup > boxListTopLevel[i].largestCSfactor)
		    boxListTopLevel[i].largestCSfactor = currGroup->maxLimitingFactorGroup;
		  
		} // END FOR groupIndex
	    } // END FOR clusterIndex
	} // END FOR chunkIndex
    } // END FOR i

  timeMeterGetMultipoleLimits.print(LOG_AREA_INTEGRALS, "Generate multipole limits for each group and for each box");






  // OK, basisFuncGroup info done for all boxes.





  Util::TimeMeter timeMeterGetDensityMatrixLimitMatrixList;

  // prepare densityMatrixLimit matrix for each level.

  csr_matrix_struct dmatLimitMatrixCSRList[numberOfLevels];
  if(dens)
    {
      ergo_real* densityMatrixLimitMatrixList[numberOfLevels];
      for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
	{
	  int nn = basisFuncGroupCounterList[levelNumber];
	  densityMatrixLimitMatrixList[levelNumber] = new ergo_real[nn*nn];
	  memset(densityMatrixLimitMatrixList[levelNumber], 0, nn*nn*sizeof(ergo_real));
	}
      for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++)
	  {
	    ergo_real currabs;
	    if(dens)
	      currabs = std::fabs(dens[i*n+j]);
	    else
	      currabs = std::fabs(ergo_CSR_get_element(densCSR, i, j));
	    for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
	      {
		int nn = basisFuncGroupCounterList[levelNumber];
		int gr_i = basisFuncGroupList[levelNumber][i];
		int gr_j = basisFuncGroupList[levelNumber][j];
		if(currabs > densityMatrixLimitMatrixList[levelNumber][gr_i*nn+gr_j])
		  densityMatrixLimitMatrixList[levelNumber][gr_i*nn+gr_j] = currabs;
	      }
	  }
      // Convert from full to CSR form
      for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
	{
	  // Create full CSR structure without symmetry
	  int nn = basisFuncGroupCounterList[levelNumber];
	  int nvalues = nn*nn;
	  std::vector<int> rowind(nvalues);
	  std::vector<int> colind(nvalues);
	  for(int i = 0; i < nn; i++)
	    for(int j = 0; j < nn; j++)
	      {
		rowind[i*nn+j] = i;
		colind[i*nn+j] = j;
	      }
	  csr_matrix_struct* currCSR = &dmatLimitMatrixCSRList[levelNumber];
	  if(ergo_CSR_create(currCSR, 
			     0,
			     nn,
			     nvalues,
			     &rowind[0],
			     &colind[0]) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_create for dmatLimitMatrixCSRList.");
	      return -1;
	    }
	  for(int i = 0; i < nn; i++)
	    for(int j = 0; j < nn; j++)
	      ergo_CSR_add_to_element(currCSR, i, j, densityMatrixLimitMatrixList[levelNumber][i*nn+j]);
	}
      for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
	delete [] densityMatrixLimitMatrixList[levelNumber];
    }
  else
    {
      output_current_memory_usage(LOG_AREA_INTEGRALS, "before calling getDmatLimitMatrixCSRList");
      if(getDmatLimitMatrixCSRList(dmatLimitMatrixCSRList, 
				   numberOfLevels, 
				   densCSR, 
				   basisFuncGroupList, 
				   basisFuncGroupCounterList) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in getDmatLimitMatrixCSRList.");
	  return -1;
	}
    }

  // OK, densityMatrixLimitMatrixList done.

  timeMeterGetDensityMatrixLimitMatrixList.print(LOG_AREA_INTEGRALS, "getting densityMatrixLimitMatrixList");
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after doing densityMatrixLimitMatrixList");






  // Crete job-list for K
  Util::TimeMeter timeMeterKjoblist;

  // compute number of jobs before allocating list.
  int noOfJobs_K_firstCount = create_joblist_exchange_for_two_boxes_recursive(basisInfo,
									      integralInfo,
									      maxNoOfMonomials,
									      J_K_params.threshold_K,
									      &boxList[0],
									      numberOfLevels,
									      dmatLimitMatrixCSRList,
									      basisFuncGroupCounterList,
									      0,
									      0, 
									      0,
									      NULL,
									      HUGE_INTEGER_NUMBER
									      );
  if(noOfJobs_K_firstCount < 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive");
      return -1;
    }

  std::vector<job_list_entry_K_struct> jobList_K(noOfJobs_K_firstCount);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating jobList_K");

  int noOfJobs_K = create_joblist_exchange_for_two_boxes_recursive(basisInfo,
								   integralInfo,
								   maxNoOfMonomials,
								   J_K_params.threshold_K,
								   &boxList[0],
								   numberOfLevels,
								   dmatLimitMatrixCSRList,
								   basisFuncGroupCounterList,
								   0,
								   0, 
								   0,
								   &jobList_K[0],
								   noOfJobs_K_firstCount
								   );
  if(noOfJobs_K < 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "job list for K created, %8i jobs", noOfJobs_K);
  timeMeterKjoblist.print(LOG_AREA_INTEGRALS, "creating job list for K");


  // Execute job-list for K

  Util::TimeMeter timeMeterK;

  int noOfThreads = J_K_params.noOfThreads_K;
  if(noOfThreads <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (noOfThreads <= 0)");
      return -1;
    }
  if(noOfThreads == 1)
    {
      // no threading requested
      if(execute_joblist_K_serial(dens,
				  densCSR,
				  basisInfo,
				  integralInfo,
				  CAM_params,
				  maxNoOfMonomials,
				  basisFuncListCount_max,
				  &boxList[0],
				  &jobList_K[0],
				  noOfJobs_K,
				  J_K_params.threshold_K,
				  K,
				  K_CSR,
				  symmetryFlag) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_K_serial");
	  return -1;
	}
    }
  else
    {
      if(execute_joblist_K_threaded(noOfThreads,
				    dens,
				    densCSR,
				    basisInfo,
				    integralInfo,
				    CAM_params,
				    maxNoOfMonomials,
				    basisFuncListCount_max,
				    &boxList[0],
				    &jobList_K[0],
				    noOfJobs_K,
				    J_K_params.threshold_K,
				    K,
				    K_CSR,
				    symmetryFlag) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_K_threaded");
	  return -1;
	}
    }

  timeMeterK.print(LOG_AREA_INTEGRALS, "Executing job list for K");
  
  for(int i = 0; i < numberOfLevels; i++)
    delete [] basisFuncGroupList[i];

  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
    ergo_CSR_destroy(&dmatLimitMatrixCSRList[levelNumber]);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after freeing stuff at end of compute_K_by_boxes");

  timeMeterTot.print(LOG_AREA_INTEGRALS, "compute_K_by_boxes");

  return 0;
}



