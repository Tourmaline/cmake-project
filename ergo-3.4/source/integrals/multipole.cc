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

/** @file multipole.cc

    \brief Code for computing multipole moments, and multipole
    interaction and translation matrices.

    @author: Elias Rudberg <em>responsible</em>. 
*/

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "multipole.h"
#include "output.h"
#include "integrals_general.h"


typedef struct {
  int l;
  int m;
} l_m_struct;

static int global_initialized_flag = 0;
static ergo_real global_lm_factor_list[MAX_MULTIPOLE_DEGREE+1][MAX_MULTIPOLE_DEGREE+1];
static l_m_struct global_l_m_list[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];

static ergo_real
slow_factorial(int n) {
  if(n == 0)
    return 1;
  return n * slow_factorial(n-1);
}

static ergo_real
get_lm_factor_slow(int l, int m) {
  return 1.0 / std::sqrt(slow_factorial(l-m)*slow_factorial(l+m));
}

static ergo_real
get_lm_factor(int l, int m) {
  if(m >= 0)
    return global_lm_factor_list[l][m];
  else
    return global_lm_factor_list[l][-m];
}

static void
get_l_m_from_index(int index, int* result_l, int* result_m) {
  int l = 0;
  int count = 0;
  while(count < index) {
    l++;
    count += 2 * l + 1;
  }
  int m = l + index - count;
  *result_l = l;
  *result_m = m;
}


int 
init_multipole_code() {
  if(global_initialized_flag)
    return 0;
  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE; l++)
    for(int m = 0; m <= l; m++)
      global_lm_factor_list[l][m] = get_lm_factor_slow(l, m);
  global_initialized_flag = 1;
  // initialize global_l_m_list
  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
    get_l_m_from_index(A, &global_l_m_list[A].l, &global_l_m_list[A].m);
  return 0;
}


int
compute_multipole_moments(const IntegralInfo& integralInfo, 
			  const DistributionSpecStruct* distr,
			  multipole_struct_small* result)
{
  int l = 0;
  int funcCountCurr_l = 0;

  int distrDegree = 0;
  for(int k = 0; k < 3; k++)
    distrDegree += distr->monomialInts[k];

  result->noOfMoments = 0;

  memset(result, 0, sizeof(multipole_struct_small));

  for(int i = 0; i < MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC; i++) {
    // get polynomial for scaled solid harmonic function
    if(i >= integralInfo.no_of_basis_func_polys) 
      throw "Error in compute_multipole_moments: (i >= integralInfo.no_of_basis_func_polys).";
    const basis_func_poly_struct* curr = &integralInfo.basis_func_poly_list[i];

    // do one term at a time
    ergo_real sum = 0;
    int savedDegree = -1;
    for(int j = 0; j < curr->noOfTerms; j++) {
      DistributionSpecStruct prim = *distr;
      int termDegree = 0;
      for(int k = 0; k < 3; k++) {
	prim.monomialInts[k] += curr->termList[j].monomialInts[k];
	termDegree += curr->termList[j].monomialInts[k];
      }
      // check degree
      if(j > 0) {
	if(termDegree != savedDegree)
	  throw "Error in compute_multipole_moments: (termDegree != savedDegree).";
    }
    savedDegree = termDegree;
    prim.coeff *= curr->termList[j].coeff;
    sum += compute_integral_of_simple_prim(&prim);
  } // END FOR j

  if(savedDegree <= distrDegree) {
    result->noOfMoments++;
    result->momentList[i] = curr->scaledSolidHarmonicPrefactor * sum;
  }
  else
    {
      // now the computed moment should be zero
      if(std::fabs(sum) > 1e-3) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: computed moment not zero when (savedDegree > distrDegree)");
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "computed moment : %22.11f", (double)sum);
	exit(0);
      }
      result->momentList[i] = 0;
    }
      
  funcCountCurr_l++;
  if(funcCountCurr_l == 2*l+1) {
    l++;
    funcCountCurr_l = 0;
  }
} // END FOR i

if(distrDegree > MAX_MULTIPOLE_DEGREE_BASIC)
  {
    // This should not happen, really, but for testing purposes it is nice to be able to set
    // MAX_MULTIPOLE_DEGREE_BASIC to be lower than needed to describe product distributions correctly.
    // The accuracy will then be bad, but the program should still work.
    result->degree = MAX_MULTIPOLE_DEGREE_BASIC;
  }
 else
   result->degree = distrDegree;

for(int k = 0; k < 3; k++)
  result->centerCoords[k] = distr->centerCoords[k];
  
return 0;
}


MMTranslator::MMTranslator()
{
  buffer_W_cc = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  buffer_W_cs = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  buffer_W_sc = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  buffer_W_ss = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];

  if(global_initialized_flag == 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS,
                "error in get_translation_matrix: (global_initialized_flag == 0)");
      exit(0);
    }
}

MMTranslator::~MMTranslator()
{
  delete []buffer_W_cc;
  delete []buffer_W_cs;
  delete []buffer_W_sc;
  delete []buffer_W_ss;
}

int
MMTranslator::getTranslationMatrix(ergo_real dx,
                                   ergo_real dy,
                                   ergo_real dz,
                                   int l_1,
                                   int l_2,
                                   ergo_real* result_W) const
{
  ergo_real r2 = dx*dx + dy*dy + dz*dz;

  // generate values of all needed "scaled regular solid harmonics"
  int largest_l_needed = MAX_MULTIPOLE_DEGREE;

  int noOf_l_values = largest_l_needed + 1;
  int L = noOf_l_values;
  ergo_real R_c[L][L];
  ergo_real R_s[L][L];
  
  R_c[0][0] = 1;
  R_s[0][0] = 0;

  // generate all R_c and R_s with l = m
  for(int l = 0; l < L-1; l++) {
    R_c[l+1][l+1] = -(dx * R_c[l][l] - dy * R_s[l][l]) / (2*l+2);
    R_s[l+1][l+1] = -(dy * R_c[l][l] + dx * R_s[l][l]) / (2*l+2);
  }
  // generate all R_c and R_s with l > m
  for(int l = 0; l < L-1; l++) {
    for(int m = 0; m <= l; m++) {
      ergo_real R_c_lmin1m = 0;
      ergo_real R_s_lmin1m = 0;
      if(l > 0 && m < l) {
	R_c_lmin1m = R_c[l-1][m];
	R_s_lmin1m = R_s[l-1][m];
      }
      R_c[l+1][m] = ((2*l+1)*dz*R_c[l][m] - r2*R_c_lmin1m) / ((l + m + 1) * (l - m + 1));
      R_s[l+1][m] = ((2*l+1)*dz*R_s[l][m] - r2*R_s_lmin1m) / ((l + m + 1) * (l - m + 1));
    }
  }

  ergo_real m1topowlist[MAX_MULTIPOLE_DEGREE+1];
  m1topowlist[0] = 1;
  for(int k = 1; k <= MAX_MULTIPOLE_DEGREE; k++)
    m1topowlist[k] = m1topowlist[k-1] * -1;
  ergo_real onehalftopowlist[2];
  onehalftopowlist[0] = 1.0;
  onehalftopowlist[1] = 0.5;

  // Use R_c and R_s values to generate translation matrix
  ergo_real (*W_cc)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_W_cc;
  ergo_real (*W_cs)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_W_cs;
  ergo_real (*W_sc)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_W_sc;
  ergo_real (*W_ss)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_W_ss;

  for(int l = 0; l <= l_1; l++)
    for(int j = 0; j <= l_2; j++)
      for(int m = 0; m <= l; m++)
	for(int k = 0; k <= j; k++) {
	  ergo_real R_c_lmjmmk = 0;
	  ergo_real R_s_lmjmmk = 0;
	  ergo_real R_c_lmjmpk = 0;
	  ergo_real R_s_lmjmpk = 0;
	  int lmj = l - j;
	  int mmk = m - k;
	  int mpk = m + k;
	  if(lmj >= 0) {
	    if(mmk >= -lmj && mmk <= lmj) {
	      if(mmk >= 0) {
		R_c_lmjmmk = R_c[lmj][mmk];
		R_s_lmjmmk = R_s[lmj][mmk];
	      }
	      else {
		R_c_lmjmmk =  m1topowlist[-mmk] * R_c[lmj][-mmk];
		R_s_lmjmmk = -m1topowlist[-mmk] * R_s[lmj][-mmk];
	      }
	    }
	    if(mpk >= -lmj && mpk <= lmj) {
	      if(mpk >= 0) {
		R_c_lmjmpk = R_c[lmj][mpk];
		R_s_lmjmpk = R_s[lmj][mpk];
	      }
	      else {
		R_c_lmjmpk =  m1topowlist[-mpk] * R_c[lmj][-mpk];
		R_s_lmjmpk = -m1topowlist[-mpk] * R_s[lmj][-mpk];
	      }
	    }
	  }

	  int dk0 = 0;
	  if(k == 0)
	    dk0 = 1;

	  W_cc[l][m][j][k] = onehalftopowlist[dk0] * ( R_c_lmjmmk + m1topowlist[k] * R_c_lmjmpk);
	  W_cs[l][m][j][k] =                         (-R_s_lmjmmk + m1topowlist[k] * R_s_lmjmpk);
	  W_sc[l][m][j][k] = onehalftopowlist[dk0] * ( R_s_lmjmmk + m1topowlist[k] * R_s_lmjmpk);
	  W_ss[l][m][j][k] =                         ( R_c_lmjmmk - m1topowlist[k] * R_c_lmjmpk);
	} // END FOR l j m k

  int noOfMoments_1 = (l_1+1)*(l_1+1);
  int noOfMoments_2 = (l_2+1)*(l_2+1);

  for(int A = 0; A < noOfMoments_1; A++)
    for(int B = 0; B < noOfMoments_2; B++)
      {
	int l = global_l_m_list[A].l;
	int m = global_l_m_list[A].m;
	int j = global_l_m_list[B].l;
	int k = global_l_m_list[B].m;

	result_W[A*noOfMoments_2+B] = 0;

	if(m >= 0 && k >= 0)
	  result_W[A*noOfMoments_2+B] = W_cc[l][m][j][k];
	if(m >= 0 && k < 0)
	  result_W[A*noOfMoments_2+B] = W_cs[l][m][j][-k];
	if(m < 0 && k >= 0)
	  result_W[A*noOfMoments_2+B] = W_sc[l][-m][j][k];
	if(m < 0 && k < 0)
	  result_W[A*noOfMoments_2+B] = W_ss[l][-m][j][-k];

	result_W[A*noOfMoments_2+B] *= get_lm_factor(j, k) / get_lm_factor(l, m);
      }

  return 0;
}


MMInteractor::MMInteractor()
{
  buffer_T_cc = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  buffer_T_cs = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  buffer_T_sc = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  buffer_T_ss = new ergo_real[MMDP1*MMDP1*MMDP1*MMDP1];
  if(global_initialized_flag == 0)
    throw "Error in MMInteractor::MMInteractor: (global_initialized_flag == 0).";
}

MMInteractor::~MMInteractor()
{
  delete [] buffer_T_cc;
  delete [] buffer_T_cs;
  delete [] buffer_T_sc;
  delete [] buffer_T_ss;
}

int
MMInteractor::getInteractionMatrix(ergo_real dx,
                                   ergo_real dy,
                                   ergo_real dz,
                                   int l_1,
                                   int l_2,
                                   ergo_real* result_T)
{
  ergo_real r2 = dx*dx + dy*dy + dz*dz;
  ergo_real oneOverR2 = (ergo_real)1 / r2;
  ergo_real oneOverR  = std::sqrt(oneOverR2);

  // generate values of all needed "scaled irregular solid harmonics"
  //int largest_l_needed = distrMultipole->degree + otherMultipole->degree;
  int largest_l_needed = (2*MAX_MULTIPOLE_DEGREE);

  int noOf_l_values = largest_l_needed + 1;
  int L = noOf_l_values;
  ergo_real I_c[L][L];
  ergo_real I_s[L][L];
  
  I_c[0][0] = oneOverR;
  I_s[0][0] = 0;

  int LL = l_1 + l_2;

  // generate all I_c and I_s with l = m
  for(int l = 0; l < LL; l++) {
    I_c[l+1][l+1] = -(2*l + 1) * oneOverR2 * (dx * I_c[l][l] - dy * I_s[l][l]);
    I_s[l+1][l+1] = -(2*l + 1) * oneOverR2 * (dy * I_c[l][l] + dx * I_s[l][l]);
  }
  // generate all I_c and I_s with l > m
  for(int l = 0; l < LL; l++) {
    for(int m = 0; m <= l; m++) {
      ergo_real I_c_lmin1m = 0;
      ergo_real I_s_lmin1m = 0;
      if(l > 0 && m < l) {
	I_c_lmin1m = I_c[l-1][m];
	I_s_lmin1m = I_s[l-1][m];
      }
      I_c[l+1][m] = oneOverR2 * ((2*l+1)*dz*I_c[l][m] - (l*l - m*m)*I_c_lmin1m);
      I_s[l+1][m] = oneOverR2 * ((2*l+1)*dz*I_s[l][m] - (l*l - m*m)*I_s_lmin1m);
    }
  }  

  ergo_real m1topowlist[MAX_MULTIPOLE_DEGREE+1];
  m1topowlist[0] = 1;
  for(int k = 1; k <= MAX_MULTIPOLE_DEGREE; k++)
    m1topowlist[k] = m1topowlist[k-1] * -1;
  ergo_real onehalftopowlist[3];
  onehalftopowlist[0] = 1.0;
  onehalftopowlist[1] = 0.5;
  onehalftopowlist[2] = 0.25;

  // Use I_c and I_s values to generate Interaction matrix
  const int MMDP1 = MAX_MULTIPOLE_DEGREE+1;
  ergo_real (*T_cc)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_T_cc;
  ergo_real (*T_cs)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_T_cs;
  ergo_real (*T_sc)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_T_sc;
  ergo_real (*T_ss)[MMDP1][MMDP1][MMDP1] = (ergo_real(*)[MMDP1][MMDP1][MMDP1])buffer_T_ss;

  for(int l = 0; l <= l_1; l++)
    for(int j = 0; j <= l_2; j++)
      for(int m = 0; m <= l; m++)
	for(int k = 0; k <= j; k++) {
	  ergo_real I_c_lpjmmk;
	  ergo_real I_s_lpjmmk;
	  if(m >= k) {
	    I_c_lpjmmk = I_c[l+j][m-k];
	    I_s_lpjmmk = I_s[l+j][m-k];
	  }
	  else {
	    I_c_lpjmmk =  m1topowlist[k-m] * I_c[l+j][k-m];
	    I_s_lpjmmk = -m1topowlist[k-m] * I_s[l+j][k-m];
	  }	
	  int dm0 = 0;
	  int dk0 = 0;
	  if(m == 0)
	    dm0 = 1;
	  if(k == 0)
	    dk0 = 1;

	  ergo_real commonPreFactor = m1topowlist[j] * onehalftopowlist[dm0 + dk0] * 2;
	  ergo_real m1topowk = m1topowlist[k];
	  ergo_real I_c_lpjmpk = I_c[l+j][m+k];
	  ergo_real I_s_lpjmpk = I_s[l+j][m+k];

	  T_cc[l][m][j][k] = commonPreFactor * ( I_c_lpjmpk + m1topowk * I_c_lpjmmk);
	  T_cs[l][m][j][k] = commonPreFactor * ( I_s_lpjmpk - m1topowk * I_s_lpjmmk);
	  T_sc[l][m][j][k] = commonPreFactor * ( I_s_lpjmpk + m1topowk * I_s_lpjmmk);
	  T_ss[l][m][j][k] = commonPreFactor * (-I_c_lpjmpk + m1topowk * I_c_lpjmmk);

	} // END FOR l j m k

  int noOfMoments_1 = (l_1+1)*(l_1+1);
  int noOfMoments_2 = (l_2+1)*(l_2+1);

  for(int A = 0; A < noOfMoments_1; A++)
    for(int B = 0; B < noOfMoments_2; B++) {
      int l = global_l_m_list[A].l;
      int m = global_l_m_list[A].m;
      int j = global_l_m_list[B].l;
      int k = global_l_m_list[B].m;

      result_T[A*noOfMoments_2+B] = 0;

      if(m >= 0 && k >= 0)
	result_T[A*noOfMoments_2+B] = T_cc[l][m][j][k];
      if(m >= 0 && k < 0)
	result_T[A*noOfMoments_2+B] = T_cs[l][m][j][-k];
      if(m < 0 && k >= 0)
	result_T[A*noOfMoments_2+B] = T_sc[l][-m][j][k];
      if(m < 0 && k < 0)
	result_T[A*noOfMoments_2+B] = T_ss[l][-m][j][-k];

      result_T[A*noOfMoments_2+B] *= get_lm_factor(l, m) * get_lm_factor(j, k);
    }

  return 0;
}


int
setup_multipole_maxAbsMomentList(multipole_struct_large* multipole) {
  ergo_real largestAbsMomentSoFar = 0;
  for(int degree = MAX_MULTIPOLE_DEGREE; degree >= 0; degree--) {
    int startIndex = degree*degree;
    int endIndex = (degree+1)*(degree+1);
    ergo_real sumOfSquares = 0;
    for(int i = startIndex; i < endIndex; i++) {
      ergo_real absMoment = std::fabs(multipole->momentList[i]);
      if(absMoment > largestAbsMomentSoFar)
	largestAbsMomentSoFar = absMoment;
      sumOfSquares += absMoment * absMoment;
    }
    multipole->maxAbsMomentList[degree] = largestAbsMomentSoFar;
    multipole->euclideanNormList[degree] = std::sqrt(sumOfSquares);
  } // END FOR degree
  return 0;
}



