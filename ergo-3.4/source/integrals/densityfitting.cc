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

/* ergo - linearly scaling DFT program.
   Copyright(c)2005 by ....
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "densityfitting.h"
#include "output.h"
#include "memorymanag.h"
#include "integrals_2el_single.h"
#include "solve_lin_eq_syst.h"
#include "utilities.h"
#include "integrals_general.h"
#include "integrals_2el.h"
#include "pi.h"
#include "boysfunction.h"

#include "mat_gblas.h"


#if 0
static void print_distr(DistributionSpecStruct* psi, const char* s)
{
  printf("print_distr for '%s':\n", s);
  printf("centerCoords: %22.11f %22.11f %22.11f\n", 
	 (double)psi->centerCoords[0],
	 (double)psi->centerCoords[1],
	 (double)psi->centerCoords[2]);
}
#endif



static ergo_real
do_2center_integral(const IntegralInfo* integralInfo, 
		    const BasisInfoStruct & basisInfo, 
		    int alpha, int beta)
{
  int i, j;
  ergo_real sum, currIntegral;

  const JK::ExchWeights CAM_params_not_used;

  const BasisFuncStruct* basisFuncAlpha = &basisInfo.basisFuncList[alpha];
  const BasisFuncStruct* basisFuncBeta  = &basisInfo.basisFuncList[beta];

  sum = 0;
  for(i = 0; i < basisFuncAlpha->noOfSimplePrimitives; i++)
    {
      DistributionSpecStruct* prim_psi1 = &basisInfo.simplePrimitiveList[basisFuncAlpha->simplePrimitiveIndex + i];
      for(j = 0; j < basisFuncBeta->noOfSimplePrimitives; j++)
	{
	  DistributionSpecStruct* prim_psi2 = &basisInfo.simplePrimitiveList[basisFuncBeta->simplePrimitiveIndex + j];
	  currIntegral = do_2e_integral_using_symb_info(CAM_params_not_used, prim_psi1, prim_psi2, *integralInfo);
	  sum += currIntegral;
	} /* END FOR j */
    } /* END FOR i */
  
  return sum;
}



#if 0
static ergo_real
do_3center_integral(const IntegralInfo* integralInfo, 
		    const BasisInfoStruct & basisInfoMain, 
		    const BasisInfoStruct & basisInfoDensFit, 
		    int alpha, int a, int b)
{
  int i, j;
  ergo_real sum, currIntegral;
  const int maxCount = 1000;
  DistributionSpecStruct list_psi1[maxCount];

  /* form product of basisfuncs mu and nu, store product in psi1 */
  int n_psi1 = get_product_simple_primitives(basisInfoMain, a,
					     basisInfoMain, b,
					     list_psi1,
					     maxCount,
					     0);
  if(n_psi1 <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives\n");
      exit(0);
      return 0;
    }

  BasisFuncStruct* basisFuncAlpha = &basisInfoDensFit.basisFuncList[alpha];

  const JK::ExchWeights CAM_params_not_used;

  sum = 0;
  for(i = 0; i < n_psi1; i++)
    {
      DistributionSpecStruct* prim_psi1 = &list_psi1[i];
      for(j = 0; j < basisFuncAlpha->noOfSimplePrimitives; j++)
	{
	  DistributionSpecStruct* prim_psi2 = &basisInfoDensFit.simplePrimitiveList[basisFuncAlpha->simplePrimitiveIndex + j];

	  //printf("before calling do_2e_integral_using_symb_info:\n");
	  //print_distr(prim_psi1, "prim_psi1");
	  //print_distr(prim_psi2, "prim_psi2");

	  currIntegral = do_2e_integral_using_symb_info(CAM_params_not_used, prim_psi1, prim_psi2, integralInfo);
	  sum += currIntegral;
	} /* END FOR j */
    } /* END FOR i */

  return sum;
}
#endif





#define MAX_NO_OF_INTEGRALS_PER_SHELL_COMB 1000
#define MAX_NO_OF_EXPPAIRS_PER_SHELL_COMB 200
#define MAX_NO_OF_XYZ_XYZ_ENTRIES 4000
#define MAX_NO_OF_ITERMLIST_ENTRIES 10000





#if 0
static void
create_subpolys_abc(const IntegralInfo* integralInfo,
		    ShellSpecStruct* shellA,
		    ShellSpecStruct* shellB,
		    ShellSpecStruct* shellC,
		    shell_pair_struct* shellPairAB,
		    int nABmax,
		    int nCmax,
		    int Nmax,
		    ergo_real* a0amg_comb_value_list_all,
		    ergo_real* alpha0ValueList)
{
  ergo_real exponentList_ab[888];
  int count_ab, count_cd;
  int i1, i2, k;
  int count;
  int expabidx, expcdidx;
  ergo_real alpha1, alpha2, alphasum, alphaproduct, alpha0, alpham, g, resultPreFactor;
  ergo_real factor, oneoveralpham, alpha0inv;
  ergo_real gtopowj[2*MAX_N1_N2_P1];
  int noOfNeeded_a0amgcombs;
  ergo_real alphaminvtopowp[2*MAX_N1_N2_P1];
  ergo_real alpha0topow_shifted[2*MAX_N1_N2_P1];
  a0amg_comb_struct* a0amg_comb_list;
  const ergo_real twoTimesPiToPow5half = 2 * std::pow(pi, 2.5);
  int pow_alpha0_min = integralInfo->n1_n2_case_info_list[nABmax][nCmax].pow_alpha0_min;
  int pow_alpha0_max = integralInfo->n1_n2_case_info_list[nABmax][nCmax].pow_alpha0_max;
  int pow_alpha0_count = pow_alpha0_max - pow_alpha0_min + 1;

  count = 0;
  for(i1 = 0; i1 < shellA->noOfContr; i1++)
    for(i2 = 0; i2 < shellB->noOfContr; i2++)
      {
	exponentList_ab[count] = shellA->exponentList[i1] + shellB->exponentList[i2];
	count++;
      }
  count_ab = count;
  count_cd = shellC->noOfContr;

  noOfNeeded_a0amgcombs = integralInfo->n1_n2_case_info_list[nABmax][nCmax].noOfNeeded_a0amgcombs;
  a0amg_comb_list = integralInfo->n1_n2_case_info_list[nABmax][nCmax].a0amg_comb_list;
	      
  for(expabidx = 0; expabidx < count_ab; expabidx++)
    for(expcdidx = 0; expcdidx < count_cd; expcdidx++)
      {
	ergo_real* a0amg_comb_value_list = &a0amg_comb_value_list_all[expabidx*MAX_NO_OF_A0AMG_COMBINATIONS*MAX_NO_OF_CONTR_GAUSSIANS+expcdidx*MAX_NO_OF_A0AMG_COMBINATIONS];
	
	alpha1 = exponentList_ab[expabidx];
	alpha2 = shellC->exponentList[expcdidx];

#if 0
	get_a0amg_comb_value_list_and_alpha0(integralInfo, alpha1, alpha2, nABmax, nCDmax,
					     a0amg_comb_value_list, &alpha0ValueList[expabidx*prep->maxNoOfExponentPairs+expcdidx]);
#else


	alphasum = alpha1 + alpha2;
	alphaproduct = alpha1 * alpha2;
	alpha0 = alphaproduct / alphasum;
		    
	alpha0ValueList[expabidx*MAX_NO_OF_CONTR_GAUSSIANS+expcdidx] = alpha0;
		    
	alpham = 0.5 * alphasum;
	g = (alpha2 - alpha1) / alphasum;
	resultPreFactor = twoTimesPiToPow5half / (alphaproduct*sqrt(alphasum));
		    
	gtopowj[0] = 1;
	gtopowj[1] = g;
		    
	oneoveralpham = 1.0 / alpham;
	factor = 1;
	for(k = 0; k <= Nmax/2; k++)
	  {
	    alphaminvtopowp[k] = factor;
	    factor *= oneoveralpham;
	  }
		    
	alpha0inv = 1.0 / alpha0;
	factor = 1;
	for(k = 0; k < -pow_alpha0_min; k++)
	  factor *= alpha0inv;
	for(k = 0; k < pow_alpha0_count; k++)
	  {
	    alpha0topow_shifted[k] = factor;
	    factor *= alpha0;
	  }
	
	for(k = 0; k < noOfNeeded_a0amgcombs; k++)
	  {
	    ergo_real factor = 1;
	    factor *= gtopowj        [ a0amg_comb_list[k].pow_g];
	    factor *= alphaminvtopowp[-a0amg_comb_list[k].pow_am];
	    factor *= alpha0topow_shifted[-pow_alpha0_min+a0amg_comb_list[k].pow_a0];
	    a0amg_comb_value_list[k] = factor * resultPreFactor;
	  } // END FOR k
#endif		    
      } /* END FOR expabidx expcdidx */
}

#endif






typedef struct
{
  int a;
  int b;
  int alpha;
} abalpha_struct;




static int 
compute_gamma_or_J_shelldriven(const BasisInfoStruct & basisInfoMain,
			       const BasisInfoStruct & basisInfoDensFit,
			       const IntegralInfo* integralInfo,
			       ergo_real* gamma,
			       ergo_real* J,
			       ergo_real* dens,
			       ergo_real* c_vector,
			       ergo_real threshold)
{
  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: compute_gamma_or_J_shelldriven used old integral code, needs to be rewritten. Right now not implemented.");
  return -1; 
 
#if 0
  int n = basisInfoMain->noOfBasisFuncs;
  int nDensFit = basisInfoDensFit.noOfBasisFuncs;
  int A, B, C; /* A B C are indeces of shells */
  int ABcount;
  int noOfMonomialsAB;
  int noOfExponentPairs_ab;
  ergo_real integralValueList[MAX_NO_OF_INTEGRALS_PER_SHELL_COMB];
  abalpha_struct integralList[MAX_NO_OF_INTEGRALS_PER_SHELL_COMB];

  int shell_ID_A;
  int shell_ID_B;
  int shell_ID_C;
  int shell_ID_A_prev;
  int shell_ID_B_prev;
  int shell_ID_C_prev;
  
  ix1_ix2_case_struct* ix1_ix2_case_AB;

  ShellSpecStruct* shellA;
  ShellSpecStruct* shellB;
  ShellSpecStruct* shellC;

  ergo_real* xyz_xyz_list_direct_ab_all;

  ergo_real* a0amg_comb_value_list_all;
  ergo_real* alpha0ValueList;
  
  itermliststruct* itermlist_ab;
  
  itermlist_value_list_struct* itermlist_ab_value_list_all_2;


  ergo_real* primitiveIntegralValueList = ergo_new(2000000,ergo_real);
  ergo_real* work_subpolyValueList3     = ergo_new(2000000,ergo_real);
  ergo_real* work_dxdydzhValueList      = ergo_new(2000000,ergo_real);
  simplified_W_poly_term_struct* work_common_termlist_for_W_list_localized =
    ergo_new(2000000,simplified_W_poly_term_struct);
  simplified_W_poly_struct* work_W_list_localized = 
    ergo_new(2000000,simplified_W_poly_struct);
  dxdydzh_term_struct* work_dxdydzh_term_list_localized = 
    ergo_new(2000000,dxdydzh_term_struct);

  char gamma_or_J_string[88];
  if(gamma != NULL)
    {
      if(J != NULL)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: both gamma and J givem");
	  return -1;
	}
      if(dens == NULL)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (dens == NULL) when computing gamma");
	  return -1;
	}
      strcpy(gamma_or_J_string, "gamma");
    }
  else
    {
      if(J == NULL)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: neither gamma nor J given");
	  return -1;
	}
      if(c_vector == NULL)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (c_vector == NULL) when computing J");
	  return -1;
	}
      strcpy(gamma_or_J_string, "J");
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "compute_gamma_or_J_shelldriven, computing %s, "
	    "threshold = %4.1g", 
	    gamma_or_J_string, (double)threshold);
  
  if(J != NULL)
    memset(J, 0, n*n*sizeof(ergo_real));
  if(gamma != NULL)
    memset(gamma, 0, nDensFit*sizeof(ergo_real));


  // allocate memory

  a0amg_comb_value_list_all =
    ergo_new(MAX_NO_OF_A0AMG_COMBINATIONS*prep->maxNoOfExponentPairs*
             MAX_NO_OF_CONTR_GAUSSIANS, ergo_real);
  alpha0ValueList =
    ergo_new(prep->maxNoOfExponentPairs*MAX_NO_OF_CONTR_GAUSSIANS, ergo_real);
  xyz_xyz_list_direct_ab_all =
    ergo_new(MAX_NO_OF_EXPPAIRS_PER_SHELL_COMB*MAX_NO_OF_XYZ_XYZ_ENTRIES,
             ergo_real);
  itermlist_ab_value_list_all_2 = 
    ergo_new(MAX_NO_OF_EXPPAIRS_PER_SHELL_COMB,itermlist_value_list_struct);

  shell_ID_A_prev = -1;
  shell_ID_B_prev = -1;
  shell_ID_C_prev = -1;

  for(ABcount = 0; ABcount < prep->noOfShellPairs; ABcount++)
    {

      shell_pair_struct* shellPairAB = &prep->shellPairList[ABcount];

      A = shellPairAB->shell_idx_1;
      B = shellPairAB->shell_idx_2;

      shellA = &basisInfoMain->shellList[A];
      shellB = &basisInfoMain->shellList[B];

      shell_ID_A = shellA->shell_ID;
      shell_ID_B = shellB->shell_ID;

      int nABmax = shellA->shellType + shellB->shellType;

      nABmax = shellA->shellType + shellB->shellType;

      noOfMonomialsAB = integralInfo->monomial_info.no_of_monomials_list[nABmax];
      noOfExponentPairs_ab = shellPairAB->noOfExponentPairs;
      
      ix1_ix2_case_AB = &integralInfo->ix1_ix2_list[shellA->shellType][shellB->shellType];
      
      compute_iterm_value_list(prep,
			       shellA->shellType, shellB->shellType, noOfExponentPairs_ab, shellPairAB,
			       xyz_xyz_list_direct_ab_all,
			       ix1_ix2_case_AB,
			       itermlist_ab_value_list_all_2);

      ergo_real screeningValue = 1;
      if(gamma != NULL)
	{
	  // we are computing gamma. Find largest absolute density matrix element for this AB shellpair.
	  int a, b;
	  ergo_real maxabs = 0;
	  for(a = 0; a < shellA->noOfBasisFuncs; a++)
	    for(b = 0; b < shellB->noOfBasisFuncs; b++)
	      {
		int a2 = shellA->startIndexInMatrix + a;
		int b2 = shellB->startIndexInMatrix + b;
		ergo_real currElementAbs = std::fabs(dens[a2*n+b2]);
		if(currElementAbs > maxabs)
		  maxabs = currElementAbs;
	      } // END FOR a b
	  screeningValue = maxabs;
	}      

      int noOfExponentPairs_ab = shellPairAB->noOfExponentPairs;

      for(C = 0; C < basisInfoDensFit.noOfShells; C++)
	{
	  shellC = &basisInfoDensFit.shellList[C];

	  shell_ID_C = shellC->shell_ID;

#if 1
	  if(shell_ID_A != shell_ID_A_prev || shell_ID_B != shell_ID_B_prev || shell_ID_C != shell_ID_C_prev)
	    {
	      setup_localized_dxdydzh_and_W_lists(integralInfo,
						  nABmax, shellC->shellType,
						  work_common_termlist_for_W_list_localized,
						  work_W_list_localized,
						  work_dxdydzh_term_list_localized);
	      int Nmax = nABmax + shellC->shellType;
	      create_subpolys_abc(integralInfo,
				  prep,
				  shellA,
				  shellB,
				  shellC,
				  shellPairAB,
				  nABmax,
				  shellC->shellType,
				  Nmax,
				  a0amg_comb_value_list_all,
				  alpha0ValueList);

	      //create_subpolys(integralInfo, prep, 
	      //	      shellA, shellB, shellC, shellD, shellPairAB, shellPairCD, nABmax, nCDmax, Nmax, a0amg_comb_value_list_all, alpha0ValueList);
	    } /* END IF compute sub-polys */
#endif	  

#if 1
	  shell_ID_A_prev = shell_ID_A;
	  shell_ID_B_prev = shell_ID_B;
	  shell_ID_C_prev = shell_ID_C;
#endif

	  

	  if(J != NULL)
	    {
	      // we are computing J. Find largest absolute c_vector element for this C shell.
	      int c;
	      ergo_real maxabs = 0;
	      for(c = 0; c < shellC->noOfBasisFuncs; c++)
		{
		  int alpha = shellC->startIndexInMatrix + c;
		  ergo_real currElementAbs = std::fabs(c_vector[alpha]);
		  if(currElementAbs > maxabs)
		    maxabs = currElementAbs;
		} // END FOR c
	      screeningValue = maxabs;
	    }

	  //printf("C = %i, shellC->noOfContr (direct) = %i\n", C, shellC->noOfContr);

	  // set up integralList
	  int a, b, c;
	  int noOfIntegrals = 0;
	  for(a = 0; a < shellA->noOfBasisFuncs; a++)
	    for(b = 0; b < shellB->noOfBasisFuncs; b++)
	      for(c = 0; c < shellC->noOfBasisFuncs; c++)
		{
		  integralList[noOfIntegrals].a = a;
		  integralList[noOfIntegrals].b = b;
		  integralList[noOfIntegrals].alpha = c;
		  integralValueList[noOfIntegrals] = 0;
		  noOfIntegrals++;
		  if(noOfIntegrals >= MAX_NO_OF_INTEGRALS_PER_SHELL_COMB)
		    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (noOfIntegrals >= MAX_NO_OF_INTEGRALS_PER_SHELL_COMB)");
		      return -1;
		    }
		} // END FOR a b c


	  int noOfMonomials_C = integralInfo->monomial_info.no_of_monomials_list[shellC->shellType];

	  int expabidx;
	  for(expabidx = 0; expabidx < noOfExponentPairs_ab; expabidx++)
	    {
	      exponent_pair_struct * currExpPair_ab;
	      currExpPair_ab = &prep->exponentPairList[shellPairAB->exponentPairIndex+expabidx];

#if 1
	      if(currExpPair_ab->sizeOfProduct * screeningValue < threshold)
		break;
#endif

	      int productOrgIdx_ab = currExpPair_ab->productOrgIdx;

	      ergo_real coeff_ab = currExpPair_ab->coeff_12;
	      ergo_real* centerCoords_ab = currExpPair_ab->centerCoords_12;

	      ergo_real* centerCoords_C = shellC->centerCoords;

	      int primC;
	      for(primC = 0; primC < shellC->noOfContr; primC++)
		{
		  ergo_real sizeOfPrimitive_C = shellC->sizeList[primC];
		  if(sizeOfPrimitive_C * currExpPair_ab->sizeOfProduct * screeningValue < threshold)
		      continue;
		  ergo_real coeff_C = shellC->coeffList[primC];
		  
		  // compute all primitive integrals needed for this case
#if 0
		  ergo_real exponent_C = shellC->exponentList[primC];
		  ergo_real alpha0 = alpha0ValueList[productOrgIdx_ab*prep->maxNoOfExponentPairs+primC];
		  ergo_real* a0amg_comb_value_list = &a0amg_comb_value_list_all[productOrgIdx_ab*MAX_NO_OF_A0AMG_COMBINATIONS*prep->maxNoOfExponentPairs+primC*MAX_NO_OF_A0AMG_COMBINATIONS];
		  get_primitive_integral_values(
						integralInfo,
						centerCoords_ab,
						centerCoords_C,
						currExpPair_ab->exponent_12,
						exponent_C,
						nABmax, shellC->shellType,
						primitiveIntegralValueList,
						work_subpolyValueList3,
						work_dxdydzhValueList,
						work_common_termlist_for_W_list_localized,
						work_W_list_localized,
						a0amg_comb_value_list,
						work_dxdydzh_term_list_localized,
						alpha0
						);
#else

		  ergo_real alpha0 = alpha0ValueList[productOrgIdx_ab*MAX_NO_OF_CONTR_GAUSSIANS+primC];
		  ergo_real* a0amg_comb_value_list = &a0amg_comb_value_list_all[productOrgIdx_ab*MAX_NO_OF_A0AMG_COMBINATIONS*MAX_NO_OF_CONTR_GAUSSIANS+primC*MAX_NO_OF_A0AMG_COMBINATIONS];
		  int n1 = nABmax;
		  int n2 = shellC->shellType;
		  ergo_real* centerCoords_1 = centerCoords_ab;
		  ergo_real* centerCoords_2 = centerCoords_C;
		  ergo_real deltaxtopow_list[3][2*MAX_N1_N2_P1];
		  ergo_real dx0, dx1, dx2, R_12_squared, factor0, factor1, factor2, arg;
		  ergo_real BoysFunctionList[2*MAX_N1_N2_P1];
		  int iloop, i, k;
		  int nNeeded = integralInfo->n1_n2_case_info_list[n1][n2].noOfSubpolysNeeded;
		  int Nmax = n1 + n2;
		  int noOfdxdydzh_needed = integralInfo->n1_n2_case_info_list[n1][n2].count_dxdydzh;
		  ergo_real factor;
		  int noOfMonomialsAB = integralInfo->monomial_info.no_of_monomials_list[n1];
		  int noOfMonomialsCD = integralInfo->monomial_info.no_of_monomials_list[n2];		  

		  subpoly_struct* subpolyList = integralInfo->n1_n2_case_info_list[n1][n2].subpolyList;

		  for(iloop = 0; iloop < nNeeded; iloop++)
		    {
		      /* compute this subPoly */
		      subpoly_struct* subPoly = &subpolyList[iloop];	      
		      ergo_real subPolySum = 0;
		      int kk;
		      for(kk = 0; kk < subPoly->noOfTerms; kk++)
			subPolySum += subPoly->termList[kk].coeff * a0amg_comb_value_list[subPoly->termList[kk].a0amg_comb_index];
		      work_subpolyValueList3[iloop] = subPolySum;
		    } /* END FOR iloop */
	  
		  dx0 = centerCoords_2[0] - centerCoords_1[0];
		  dx1 = centerCoords_2[1] - centerCoords_1[1];
		  dx2 = centerCoords_2[2] - centerCoords_1[2];
		  R_12_squared = dx0*dx0 + dx1*dx1 + dx2*dx2;

		  factor0 = dx0;
		  factor1 = dx1;
		  factor2 = dx2;
		  deltaxtopow_list[0][0] = 1;
		  deltaxtopow_list[0][1] = dx0;
		  deltaxtopow_list[1][0] = 1;
		  deltaxtopow_list[1][1] = dx1;
		  deltaxtopow_list[2][0] = 1;
		  deltaxtopow_list[2][1] = dx2;
		  for(k = 2; k <= Nmax; k++)
		    {
		      factor0 *= dx0;
		      factor1 *= dx1;
		      factor2 *= dx2;
		      deltaxtopow_list[0][k] = factor0;
		      deltaxtopow_list[1][k] = factor1;
		      deltaxtopow_list[2][k] = factor2;
		    }

		  /* Compute all Boys function values needed */
		  /* Use downward recursion to get Boys function values */
		  arg =  alpha0 * R_12_squared;
		  BoysFunctionList[Nmax] = BoysFunction(Nmax, arg);
		  for(i = Nmax-1; i >= 0; i--)
		    BoysFunctionList[i] = (2*arg*BoysFunctionList[i+1] + std::exp(-arg)) / (2*i+1);

		  /* prepare list of all dxdydzh values needed for this case */
		  for(iloop = 0; iloop < noOfdxdydzh_needed; iloop++)
		    {
		      dxdydzh_term_struct* term;
		      /* compute this combination of dx dy dz h */
		      term = &work_dxdydzh_term_list_localized[iloop];
		      factor = 
			deltaxtopow_list[0][term->pows[0]] * 
			deltaxtopow_list[1][term->pows[1]] *
			deltaxtopow_list[2][term->pows[2]] *
			BoysFunctionList   [term->pows[3]];
		      work_dxdydzhValueList[iloop] = factor;
		    } /* END FOR i */

		  prep_products_of_Ix6_polys(noOfMonomialsAB, noOfMonomialsCD, integralInfo->W_list_3, 
					     work_subpolyValueList3, work_dxdydzhValueList, work_W_list_localized, primitiveIntegralValueList);


#endif

		  
		  int integralNo;
		  for(integralNo = 0; integralNo < noOfIntegrals; integralNo++)
		    {
		      ergo_real sum = 0;
		      int a2 = shellA->startIndexInMatrix + integralList[integralNo].a;
		      int b2 = shellB->startIndexInMatrix + integralList[integralNo].b;
		      int c2 = shellC->startIndexInMatrix + integralList[integralNo].alpha;
		      int basfuncpolyidA = basisInfoMain->basisFuncPolyIDStructList[a2].id;
		      int basfuncpolyidB = basisInfoMain->basisFuncPolyIDStructList[b2].id;
		      int basfuncpolyidC = basisInfoDensFit.basisFuncPolyIDStructList[c2].id;
		      int n_ab = prep->itermlist_list[basfuncpolyidA][basfuncpolyidB].n;
		      itermlist_ab = prep->itermlist_list[basfuncpolyidA][basfuncpolyidB].itermlist;

		      basis_func_poly_struct* poly_C = &integralInfo->basis_func_poly_list[basfuncpolyidC];

		      int basfuncidsAB = basfuncpolyidA*MAX_NO_OF_BASIS_FUNC_POLYS+basfuncpolyidB;
		      ergo_real* itermlist_ab_val_list = itermlist_ab_value_list_all_2[expabidx].valueList[basfuncidsAB];


		      int mab;
		      for(mab = 0; mab < n_ab; mab++)
			{
			  int monomialID_ab = itermlist_ab[mab].monomialIndex;
			  int mc;
			  for(mc = 0; mc < poly_C->noOfTerms; mc++)
			    {
			      int monomialID_C = poly_C->termList[mc].monomialID;
			      sum += 
				coeff_ab * 
				coeff_C *
				itermlist_ab_val_list[mab] *
				poly_C->termList[mc].coeff *
				primitiveIntegralValueList[monomialID_ab*noOfMonomials_C + monomialID_C];
			    } // END FOR mc

			} // END FOR mab
		      integralValueList[integralNo] += sum;
		    } // END FOR integralNo

		} // END FOR primC

	    } /* END FOR exp ab */

	  // integrals computed. now add to gamma or J.
	  if(gamma != NULL)
	    {
	      int integralNo;
	      for(integralNo = 0; integralNo < noOfIntegrals; integralNo++)
		{
		  int a2 = shellA->startIndexInMatrix + integralList[integralNo].a;
		  int b2 = shellB->startIndexInMatrix + integralList[integralNo].b;
		  int c2 = shellC->startIndexInMatrix + integralList[integralNo].alpha;
		  if(A == B)
		    //if(a2 == b2)
		    gamma[c2] += integralValueList[integralNo] * dens[a2*n+b2];
		  else
		    gamma[c2] += 2 * integralValueList[integralNo] * dens[a2*n+b2];
		} // END FOR integralNo
	    } // END IF add to gamma
	  if(J != NULL)
	    {
	      int integralNo;
	      for(integralNo = 0; integralNo < noOfIntegrals; integralNo++)
		{
		  int a2 = shellA->startIndexInMatrix + integralList[integralNo].a;
		  int b2 = shellB->startIndexInMatrix + integralList[integralNo].b;
		  int c2 = shellC->startIndexInMatrix + integralList[integralNo].alpha;
		  if(A == B)
		    {
		      J[a2*n+b2] += integralValueList[integralNo] * c_vector[c2];
		    }
		  else
		    {
		      J[a2*n+b2] += integralValueList[integralNo] * c_vector[c2];
		      J[b2*n+a2] += integralValueList[integralNo] * c_vector[c2];
		    }
		} // END FOR integralNo
	    } // END IF add to J

	} /* END FOR C */
    } /* END FOR ABcount */


  ergo_free(primitiveIntegralValueList);
  ergo_free(work_subpolyValueList3);
  ergo_free(work_dxdydzhValueList);
  ergo_free(work_common_termlist_for_W_list_localized);
  ergo_free(work_W_list_localized);
  ergo_free(work_dxdydzh_term_list_localized);
  ergo_free(a0amg_comb_value_list_all);
  ergo_free(alpha0ValueList);
  ergo_free(xyz_xyz_list_direct_ab_all);
  ergo_free(itermlist_ab_value_list_all_2);

  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "compute_gamma_or_J_shelldriven ending OK.");

  return 0;

#endif
}










int 
densfit_compute_gamma(const IntegralInfo* integralInfo,
		      const BasisInfoStruct & basisInfoMain,
		      const BasisInfoStruct & basisInfoDensFit,
		      ergo_real* densityMatrix,
		      ergo_real* result_gamma,
		      ergo_real threshold)
{
  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "entering densfit_compute_gamma, %i aux. b.f.s",
            basisInfoDensFit.noOfBasisFuncs);

  if(basisInfoDensFit.noOfBasisFuncs <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in densfit_compute_gamma: (basisInfoDensFit.noOfBasisFuncs <= 0)");
      return -1;
    }

#if 0
  int n = basisInfoMain->noOfBasisFuncs;
  // construct vector gamma, one element at a time
  int alpha;
  for(alpha = 0; alpha < basisInfoDensFit.noOfBasisFuncs; alpha++)
    {
      ergo_real sum = 0;
      int a, b;
      for(a = 0; a < n; a++)
	for(b = 0; b < n; b++)
	  {
	    ergo_real integral = do_3center_integral(integralInfo, basisInfoMain, basisInfoDensFit, alpha, a, b);
	    //printf("integral = %22.11f\n", integral);
	    sum += integral * densityMatrix[a*n+b];
	  } // END FOR a b
      result_gamma[alpha] = sum;
    } // END FOR alpha
#if 0
  printf("gamma (simple):\n");
  int ii;
  for(ii = 0; ii < basisInfoDensFit.noOfBasisFuncs; ii++)
    printf("%22.11f\n", result_gamma[ii]);
  exit(0);
#endif
#else
  if(compute_gamma_or_J_shelldriven(basisInfoMain,
				    basisInfoDensFit,
				    integralInfo,
				    result_gamma,
				    NULL,
				    densityMatrix,
				    NULL,
				    threshold) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_gamma_or_J_shelldriven when computing gamma");
      return -1;
    }
#endif

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "densfit_compute_gamma ending OK.");
  timeMeter.print(LOG_AREA_INTEGRALS, "densfit_compute_gamma");

  return 0;
}



int
densfit_compute_alpha_beta_matrix_inverse(const IntegralInfo* integralInfo,
					  const BasisInfoStruct & basisInfoDensFit,
					  ergo_real* result_U_inverse)
{
  Util::TimeMeter timeMeter;
  int alpha, beta;
  int n = basisInfoDensFit.noOfBasisFuncs;
  
  // construct one triangle of matrix ( alpha | beta )
  for(alpha = 0; alpha < n; alpha++)
    for(beta = 0; beta <= alpha; beta++)
      {
	result_U_inverse[alpha*n+beta] = do_2center_integral(integralInfo, basisInfoDensFit, alpha, beta);
      } // END FOR alpha beta
  // fill the other triangle
  for(alpha = 0; alpha < n; alpha++)
    for(beta = alpha+1; beta < n; beta++)
      {
	result_U_inverse[alpha*n+beta] = result_U_inverse[beta*n+alpha];
      }
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "Matrix (alpha|beta) computed, inverting...");
  
  //  char uplo = 'U';
  int info = -1;
  template_lapack_potf2("U", &n, result_U_inverse, &n, &info);
  if( info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS,  "%s: error in potf2, info=%d", __FUNCTION__, info);
      return -1;
    }

  // set other triangle to zero
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	if(i < j)
	  result_U_inverse[i*n+j] = 0;
      }

  // Get Uinv
  //  char diag = 'N';
  template_lapack_trtri("U", "N", &n, result_U_inverse, &n, &info);
  
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in trtri");
      return -1;
    }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "densfit_compute_alpha_beta_matrix_inverse ending OK.");
  timeMeter.print(LOG_AREA_INTEGRALS, "densfit_compute_alpha_beta_matrix_inverse");

  return 0;
}

/* AlphaBetaMemSzLimit: largest dimension of the (alpha|beta) array
 * that is kept in RAM. */
//static const size_t ALPHA_BETA_MEM_SZ_LIMIT = 0x8000;
static const size_t ALPHA_BETA_MEM_SZ_LIMIT = 40000;
DensfitData*
densfit_init(const IntegralInfo* integralInfo,
             const BasisInfoStruct & basisInfoDensFit)
{
  int alpha, beta;
  int n = basisInfoDensFit.noOfBasisFuncs;
  ergo_real *U = ergo_new(n*n, ergo_real);
  DensfitData *d = ergo_new(1, DensfitData);

  // construct one triangle of matrix ( alpha | beta )
  for(alpha = 0; alpha < n; alpha++)
    for(beta = 0; beta <= alpha; beta++)
      U[alpha*n+beta] =
        do_2center_integral(integralInfo, basisInfoDensFit, alpha, beta);
 /* 2GB limit on storing stuff in memory */
  d->using_file = n*sizeof(ergo_real)> ALPHA_BETA_MEM_SZ_LIMIT;
  if(d->using_file) {
    char fname[256];
    snprintf(fname, sizeof(fname), "/scratch/alphabeta.%d", getpid());
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "(alpha|beta) too large  - using a temp file");
    d->f = fopen(fname, "w+b");
    if(!d->f) goto err;
    unlink(fname); /* in case of crashes... */
    for(alpha = 0; alpha < n; alpha++)
      if(fwrite(U + alpha*n, sizeof(ergo_real), alpha+1, d->f) != 
	 (size_t)alpha+1)
	goto err;
    ergo_free(U);
  } else d->ptr = U;
  return d;
 err:
  ergo_free(U);
  ergo_free(d);
  return NULL;
}

void
densfit_destroy(DensfitData *d)
{
  if(d->using_file) fclose(d->f);
  else ergo_free(d->ptr);
  ergo_free(d);
}

/* compute C vector by direct linear equation solving */
int
densfit_compute_c_vector(const IntegralInfo* integralInfo,
			 const BasisInfoStruct & basisInfoDensFit,
			 DensfitData* df_data,
			 ergo_real* gamma,
			 ergo_real* result_c_vector)
{
  //  static int ONEI = 1;
  int alpha;
  int n = basisInfoDensFit.noOfBasisFuncs, info;
  ergo_real *U = ergo_new(n*n, ergo_real);

  memcpy(result_c_vector, gamma, n*sizeof(ergo_real));
  if(df_data->using_file) {
    rewind(df_data->f);
    for(alpha = 0; alpha < n; alpha++) {
      info = fread(U + alpha*n, sizeof(ergo_real), alpha+1, df_data->f);
      if( info != alpha+1) {
	do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "densfit_compute_c_vector_lr: IO failed, Expected %d read %d",
		  alpha+1, info);
	return -1;
      }
    }
  } else
    memcpy(U, df_data->ptr, n*n*sizeof(ergo_real));

  
  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in densfit_compute_c_vector: posv not implemented.");
  return -1;
  //info=dposv_wrapper('U', n, 1, U, n, result_c_vector, n);

  ergo_free(U);
  if( info != 0) {
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "densfit_compute_c_vector_lr failed (info=%d)", info);
    return -1;
  } else {
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "densfit_compute_c_vector_lr ending OK.");
    return 0;
  }
}


int
densfit_compute_J(const IntegralInfo* integralInfo,
		  const BasisInfoStruct & basisInfoMain,
		  const BasisInfoStruct & basisInfoDensFit,
		  ergo_real* c_vector,
		  ergo_real* result_J,
		  ergo_real threshold)
{
  Util::TimeMeter timeMeter;

#if 0
  int a, b, n, alpha;
  n = basisInfoMain->noOfBasisFuncs;
  for(a = 0; a < n; a++)
    for(b = 0; b < n; b++)
      {
	ergo_real sum = 0;
	for(alpha = 0; alpha < basisInfoDensFit.noOfBasisFuncs; alpha++)
	  {
	    sum += do_3center_integral(integralInfo, basisInfoMain, basisInfoDensFit, alpha, a, b) * c_vector[alpha];
	  } // END FOR alpha
	result_J[a*n+b] = sum;
      } // END FOR a b
#else
  if(compute_gamma_or_J_shelldriven(basisInfoMain,
				    basisInfoDensFit,
				    integralInfo,
				    NULL,
				    result_J,
				    NULL,
				    c_vector,
				    threshold) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_gamma_or_J_shelldriven when computing J");
      return -1;
    }
#endif

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "densfit_compute_J ending OK.");
  timeMeter.print(LOG_AREA_INTEGRALS, "densfit_compute_J");

  return 0;
}

