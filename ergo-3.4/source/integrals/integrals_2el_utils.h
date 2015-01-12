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

#ifndef INTEGRALS_2EL_UTILS_HEADER
#define INTEGRALS_2EL_UTILS_HEADER


#include "organize_distrs.h"
#include "box_system.h"


struct distr_list_description_struct {
  distr_org_struct org;
  ergo_real totCharge;
  ergo_real multipolePoint[3];
  multipole_struct_large multipole;
  ergo_real maxMomentVectorNormForDistrsList[MAX_MULTIPOLE_DEGREE_BASIC+1];
  distr_list_description_struct();
};


#define MAX_NO_OF_BRANCHES 10


typedef struct
{
  int basisFuncGroupIndex;
  ergo_real max_CS_factor;
  ergo_real maxMomentVectorNormList[MAX_MULTIPOLE_DEGREE_BASIC+1];
  int maxMultipoleDegree;
} basis_func_group_info_for_box;


struct box_struct {
  box_struct_basic basicBox;
  distr_list_description_struct distrListForK;
  distr_list_description_struct branchList[MAX_NO_OF_BRANCHES];
  ergo_real multipolePoint[3];
  int branchIndexList[MAX_NO_OF_BRANCHES];
  int branchCountList[MAX_NO_OF_BRANCHES];
  int noOfBasisFuncs;
  int noOfRelevantBasisFuncGroups;
  basis_func_group_info_for_box* basisFuncGroupInfoList;
  ergo_real largestCSfactor; // Largest Cauchy-Schwartz factor in box
  ergo_real multipoleEuclideanNormList[MAX_MULTIPOLE_DEGREE_BASIC+1];
  box_struct();
};


typedef struct
{
  ergo_real* summedIntegralList;
  ergo_real* primitiveIntegralList;
  ergo_real* primitiveIntegralList_work;
  ergo_real* partial_dmat_1;
  ergo_real* partial_dmat_2; // used in non-symmetric case
  ergo_real* partial_K_1;
  ergo_real* partial_K_2; // used in non-symmetric case
} JK_contribs_buffer_struct;




ergo_real get_max_abs_vector_element(int n, const ergo_real* vector);

void
allocate_buffers_needed_by_integral_code(const IntegralInfo & integralInfo, 
					 int maxNoOfMonomials,
					 int basisFuncListCount_max,
					 JK_contribs_buffer_struct* bufferStruct);

void
free_buffers_needed_by_integral_code(JK_contribs_buffer_struct* bufferStruct);

int
get_related_integrals_h(const IntegralInfo & integralInfo,
			const JK::ExchWeights & CAM_params,
			int n1max, int noOfMonomials_1,
			int n2max, int noOfMonomials_2,
			ergo_real dx0, 
			ergo_real dx1, 
			ergo_real dx2, 
			ergo_real alpha1, 
			ergo_real alpha2,
			ergo_real alpha0,
			ergo_real* primitiveIntegralList,
			ergo_real* primitiveIntegralList_work,
			ergo_real resultPreFactor);

void
compute_extent_for_list_of_distributions(int n,
					 DistributionSpecStructLabeled* distrList,
					 ergo_real threshold,
					 ergo_real maxLimitingFactor,
					 ergo_real maxabsDmatelement);

int
get_list_of_labeled_distrs_maxLimitingFactor(const BasisInfoStruct & basisInfo,
					     const IntegralInfo & integralInfo,
					     ergo_real threshold,
					     ergo_real* resultMaxLimitingFactor,
					     ergo_real maxDensityMatrixElement);

int
get_list_of_labeled_distrs(const BasisInfoStruct & basisInfo,
			   const IntegralInfo & integralInfo,
			   ergo_real threshold,
			   DistributionSpecStructLabeled* resultList,
			   int maxCountDistrs,
			   ergo_real maxLimitingFactor,
			   const ergo_real* dens,
			   ergo_real maxDensityMatrixElement);

int
create_box_system_and_reorder_distrs(int distrCount,
				     DistributionSpecStructLabeled* distrList,
				     ergo_real toplevelBoxSize,
				     BoxSystem & boxSystem);



#endif
