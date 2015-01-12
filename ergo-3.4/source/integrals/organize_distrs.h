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

#ifndef ORGANIZE_DISTRS_HEADER
#define ORGANIZE_DISTRS_HEADER

#include "output.h"
#include "multipole.h"
#include "simple_sparse_mat.h"

#include <vector>


typedef struct
{
  int startIndex;
  int distrCount;
  int nmax;
  ergo_real centerCoords[3];
  ergo_real exponent;
  ergo_real maxSizeGroup;
  ergo_real maxExtentGroup;
  ergo_real maxLimitingFactorGroup;
  ergo_real maxAbsDmatElementGroup;
  multipole_struct_small* multipolePtr;
  ergo_real multipoleEuclideanNormList[MAX_MULTIPOLE_DEGREE_BASIC+1];
} distr_group_struct;

typedef struct
{
  int basisFuncPairIndex;
  int monomialIndex;
  ergo_real coeff;
} minimal_distr_struct;

typedef struct
{
  int nmax;
  ergo_real exponent;
  int groupStartIndex;
  int noOfGroups;
  ergo_real maxLimitingFactorForCluster;
  ergo_real multipoleEuclideanNormList[MAX_MULTIPOLE_DEGREE_BASIC+1];
} cluster_struct;

typedef struct
{
  int index_1;
  int index_2;
  int index_1_mod;
  int index_2_mod;
  int index_inbox_1;
  int index_inbox_2;
  int pairIndex;
  ergo_real dmatElement;
} basis_func_pair_struct;

#ifndef BASIS_FUNC_POLY_MAX_DEGREE
#error The constant BASIS_FUNC_POLY_MAX_DEGREE must be defined.
#endif
#if BASIS_FUNC_POLY_MAX_DEGREE<6
#define MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK 1000
#else
#define MAX_NO_OF_BASIS_FUNC_PAIRS_PER_CHUNK 10000
#endif

typedef struct
{
  int clusterStartIndex;
  int noOfClusters;
  int noOfBasisFuncPairs;
  int basisFuncPairListIndex;
  int basisFuncForChunksIndex;
  int basisFuncForChunkCount;
  int global_debug_id;
} chunk_struct;


struct distr_org_struct {
  std::vector<minimal_distr_struct> minimalDistrList;
  std::vector<distr_group_struct> groupList;
  std::vector<cluster_struct> clusterList;
  std::vector<chunk_struct> chunkList;
  std::vector<basis_func_pair_struct> basisFuncPairList;
  std::vector<int> basisFuncListForChunks;
  std::vector<int> basisFuncListForChunks_map;
  std::vector<int> basisFuncList;
  std::vector<i_j_val_struct> spMatElementList;
  std::vector<int> spMatCountList;
  std::vector<int> spMatIdxList;
  int minimalDistrCount;
  int groupCount;
  int clusterCount;
  int chunkCount;
  int basisFuncPairCount;
  int basisFuncForChunksCount;
  int basisFuncListCount;
  ergo_real maxExtent;
  ergo_real maxDistanceOutsideBox;
  distr_org_struct():
    minimalDistrCount(0), 
    groupCount(0), 
    clusterCount(0), 
    chunkCount(0), 
    basisFuncPairCount(0), 
    basisFuncForChunksCount(0), 
    basisFuncListCount(0), 
    maxExtent(0), 
    maxDistanceOutsideBox(0)
  {}
};




int
organize_distributions(const IntegralInfo & integralInfo,
		       DistributionSpecStructLabeled* distrList_in, 
		       int distrCount, 
		       distr_org_struct* result,
		       const ergo_real* boxCenterCoords,
		       ergo_real boxWidth);

#endif
