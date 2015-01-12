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

/** @file ci.cc

    \brief Configuration Interaction (CI) code.

    @author: Elias Rudberg <em>responsible</em>. 
*/
#include <cstdlib>
#include <cstdio>
#include <memory.h>
#include <assert.h>
#include <vector>

#include "ci.h"
#include "output.h"
#include "utilities.h"
#include "integrals_2el_explicit.h"
#include "densfromf_full.h"

#include "../matrix/mat_gblas.h"


const int MAX_AOS = 30;
const int MAX_SOS = 2 * MAX_AOS;
const int MAX_ELECTRONS = 40;

const int SPIN_A = 1;
const int SPIN_B = 2;

typedef struct
{
  ergo_real x[MAX_AOS][MAX_AOS][MAX_AOS][MAX_AOS];
} four_idx_AO_struct;

typedef struct
{
  ergo_real x[MAX_SOS][MAX_SOS][MAX_SOS][MAX_SOS];
} four_idx_SO_struct;

typedef struct
{
  ergo_real x[MAX_SOS][MAX_SOS];
} two_idx_SO_struct;

typedef struct
{
  ergo_real coeffs[MAX_AOS];
  int spin;
} SO_struct;

typedef struct
{
  char SO_list[MAX_ELECTRONS];
  int startIndex; // at next level
  int count; // at next level
} SlaterDet_struct;

typedef struct
{
  int a;
  int b;
  int nDiff;
  char SOs_a[2];
  char SOs_b[2];
  char SOs_a_pos[2];
  char SOs_b_pos[2];
} SlaterDet_pair_struct;


#if 0
static void printPair(SlaterDet_pair_struct* p)
{
  printf("a b nDiff SOs_a_0 SOs_a_1 SOs_b_0 SOs_b_1 SOs_a_pos_0 SOs_a_pos_1 SOs_b_pos_0 SOs_b_pos_1 : %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i\n",
	 p->a, p->b, p->nDiff, p->SOs_a[0], p->SOs_a[1], p->SOs_b[0], p->SOs_b[1], p->SOs_a_pos[0], p->SOs_a_pos[1], p->SOs_b_pos[0], p->SOs_b_pos[1]);
}
#endif


static ergo_real get_vector_norm(int n, const ergo_real* v)
{
  ergo_real sqSum = 0;
  for(int i = 0; i < n; i++)
    sqSum += v[i] * v[i];
  return std::sqrt(sqSum);
}

static void normalize_vector(int n, ergo_real* v)
{
  ergo_real factor = 1.0 / get_vector_norm(n, v);
  for(int i = 0; i < n; i++)
    v[i] *= factor;
}

#if 0
static void print_large_coeffs(int n, const ergo_real* coeffList)
{
  const ergo_real limit = 0.01;
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "printing all coeffs larger than %6.3f", (double)limit);
  const int smax = 888;
  char s[smax];
  s[0] = 0;
  for(int i = 0; i < n; i++)
    {
      if(coeffList[i] > limit)
	{
	  char stmp[88];
	  sprintf(stmp, "%7.4f (%9i) ", (double)coeffList[i], i);
	  strcat(s, stmp);
	}
      if(strlen(s) > 99)
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_CI, s);
	  s[0] = 0;
	}
    }
  if(strlen(s) > 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, s);
      s[0] = 0;
    }
}
#endif


void get_1el_energy_and_gradient(int nSOs,
				 int nEl,
				 int nSlaterDets, 
				 const SlaterDet_struct* SlaterDetList, 
				 int nSlaterDetPairs,
				 const SlaterDet_pair_struct* SlaterDet_pair_list,
				 const int* pairCountList,
				 int noOfTrialVectors,
				 ergo_real* energy_list,
				 ergo_real** coeffListList,
				 const two_idx_SO_struct* h_SO, 
				 ergo_real** resultGradient_list)
{
  for(int k = 0; k < noOfTrialVectors; k++)
    {
      energy_list[k] = 0;
      for(int i = 0; i < nSlaterDets; i++)
	resultGradient_list[k][i] = 0;
    }
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {      

      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  // Do nothing here!
	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{

	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];

	  int count, savedCount;
	  int signa = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signa *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;

	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];

	  int signb = signa;
	  // annihilate q from b
	  char SO_list_mod_b_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[b].SO_list[i] == q)
		savedCount = count;
	      else
		{
		  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signb *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
		  
	  int equal = 1;
	  for(int i = 0; i < nEl-1; i++)
	    {
	      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		equal = 0;
	    }
		  
	  if(equal)
	    {
	      // OK, we have a contribution
	      ergo_real x = signb;
	      x *= 2; // To account for a > b case
	      ergo_real h_pq = h_SO->x[p][q];
	      for(int k = 0; k < noOfTrialVectors; k++)
		{
		  energy_list[k] += x * h_pq * coeffListList[k][a] * coeffListList[k][b];
		  resultGradient_list[k][a] += x * h_pq * coeffListList[k][b];
		  resultGradient_list[k][b] += x * h_pq * coeffListList[k][a];
		}
	    }
	}
      else
	{
	  
	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  for(int p = 0; p < nSOs; p++)
	    {
	      int count, savedCount;
	      int signa = 1;
	      // annihilate p from a
	      char SO_list_mod_a_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(int i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[a].SO_list[i] == p)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		      count++;
		      if(savedCount < 0)
			signa *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;
	      for(int q = 0; q < nSOs; q++)
		{
		  int signb = signa;
		  // annihilate q from b
		  char SO_list_mod_b_1[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(int i = 0; i < nEl; i++)
		    {
		      if(SlaterDetList[b].SO_list[i] == q)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
			  count++;
			  if(savedCount < 0)
			    signb *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  
		  int equal = 1;
		  for(int i = 0; i < nEl-1; i++)
		    {
		      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
			equal = 0;
		    }
		  
		  if(equal)
		    {
		      // OK, we have a contribution
		      ergo_real x = signb;
		      ergo_real h_pq = h_SO->x[p][q];
		      for(int k = 0; k < noOfTrialVectors; k++)
			{
			  energy_list[k] += x * h_pq * coeffListList[k][a] * coeffListList[k][b];
			  resultGradient_list[k][a] += x * h_pq * coeffListList[k][b];
			  resultGradient_list[k][b] += x * h_pq * coeffListList[k][a];
			}
		    }
		} // END FOR q
	    } // END FOR p
	} // END ELSE
    } // END FOR pairIdx
}




void get_1el_contribs_to_mult_or_dmat(int nSOs,
				      int nEl,
				      int nSlaterDets, 
				      const SlaterDet_struct* SlaterDetList, 
				      const SlaterDet_pair_struct* SlaterDetPair,
				      const two_idx_SO_struct* h_SO, 
				      const ergo_real* sourceVector,
				      ergo_real* resultVector, // if result of matrix-vector mult is requested
				      two_idx_SO_struct* resultdmat // if dmat is requested
				      )
{
  if(SlaterDetPair->nDiff == 2)
    {
      // Do nothing here!
    }
  else if(SlaterDetPair->nDiff == 1)
    {

      int a = SlaterDetPair->a;
      int b = SlaterDetPair->b;

      int p = SlaterDetPair->SOs_a[0];

      int count, savedCount;
      int signa = 1;
      // annihilate p from a
      char SO_list_mod_a_1[MAX_ELECTRONS];
      count = 0;
      savedCount = -1;
      for(int i = 0; i < nEl; i++)
	{
	  if(SlaterDetList[a].SO_list[i] == p)
	    savedCount = count;
	  else
	    {
	      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
	      count++;
	      if(savedCount < 0)
		signa *= -1;
	    }
	}
      if(savedCount < 0)
	return;

      int q = SlaterDetPair->SOs_b[0];

      int signb = signa;
      // annihilate q from b
      char SO_list_mod_b_1[MAX_ELECTRONS];
      count = 0;
      savedCount = -1;
      for(int i = 0; i < nEl; i++)
	{
	  if(SlaterDetList[b].SO_list[i] == q)
	    savedCount = count;
	  else
	    {
	      SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
	      count++;
	      if(savedCount < 0)
		signb *= -1;
	    }
	}
      if(savedCount < 0)
	return;
		  
      int equal = 1;
      for(int i = 0; i < nEl-1; i++)
	{
	  if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
	    equal = 0;
	}
		  
      if(equal)
	{
	  // OK, we have a contribution
	  ergo_real x = signb;
	  ergo_real h_pq = h_SO->x[p][q];
	  
	  ergo_real resultMatrix_ab = x * h_pq;
	  ergo_real resultMatrix_ba = x * h_pq;

	  if(resultVector)
	    {
	      resultVector[a] += resultMatrix_ab * sourceVector[b];
	      resultVector[b] += resultMatrix_ba * sourceVector[a];
	    }
	  else
	    {
	      resultdmat->x[p][q] += x * sourceVector[a] * sourceVector[b];
	      resultdmat->x[q][p] += x * sourceVector[a] * sourceVector[b];
	    }
	}
    }
  else
    {
	  
      int a = SlaterDetPair->a;
      int b = SlaterDetPair->b;

      for(int p = 0; p < nSOs; p++)
	{
	  int count, savedCount;
	  int signa = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signa *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
	  for(int q = 0; q < nSOs; q++)
	    {
	      int signb = signa;
	      // annihilate q from b
	      char SO_list_mod_b_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(int i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[b].SO_list[i] == q)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		      count++;
		      if(savedCount < 0)
			signb *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;
		  
	      int equal = 1;
	      for(int i = 0; i < nEl-1; i++)
		{
		  if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		    equal = 0;
		}
		  
	      if(equal)
		{
		  // OK, we have a contribution
		  ergo_real x = signb;
		  ergo_real h_pq = h_SO->x[p][q];

		  ergo_real resultMatrix_ab = x * h_pq;

		  if(resultVector)
		    resultVector[a] += resultMatrix_ab * sourceVector[b];
		  else
		    resultdmat->x[p][q] += x * sourceVector[a] * sourceVector[b];
		}
	    } // END FOR q
	} // END FOR p
    } // END ELSE
}






void get_1el_contribs(int nSOs,
		      int nEl,
		      int nSlaterDets, 
		      const SlaterDet_struct* SlaterDetList, 
		      int nSlaterDetPairs,
		      const SlaterDet_pair_struct* SlaterDet_pair_list,
		      const int* pairCountList,
		      const two_idx_SO_struct* h_SO, 
		      ergo_real* resultMatrix)
{
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {            
      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  // Do nothing here!
	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{

	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];

	  int count, savedCount;
	  int signa = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signa *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;

	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];

	  int signb = signa;
	  // annihilate q from b
	  char SO_list_mod_b_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[b].SO_list[i] == q)
		savedCount = count;
	      else
		{
		  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signb *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
		  
	  int equal = 1;
	  for(int i = 0; i < nEl-1; i++)
	    {
	      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		equal = 0;
	    }
		  
	  if(equal)
	    {
	      // OK, we have a contribution
	      ergo_real x = signb;
	      //x *= 2; // To account for a > b case
	      ergo_real h_pq = h_SO->x[p][q];

	      //printf("1el contrib to a b = %i %i : %22.11f\n", a, b, x * h_pq);

	      resultMatrix[a*nSlaterDets+b] += x * h_pq;
	      resultMatrix[b*nSlaterDets+a] += x * h_pq;
	    }
	}
      else
	{
	  
	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  for(int p = 0; p < nSOs; p++)
	    {
	      int count, savedCount;
	      int signa = 1;
	      // annihilate p from a
	      char SO_list_mod_a_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(int i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[a].SO_list[i] == p)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		      count++;
		      if(savedCount < 0)
			signa *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;
	      for(int q = 0; q < nSOs; q++)
		{
		  int signb = signa;
		  // annihilate q from b
		  char SO_list_mod_b_1[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(int i = 0; i < nEl; i++)
		    {
		      if(SlaterDetList[b].SO_list[i] == q)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
			  count++;
			  if(savedCount < 0)
			    signb *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  
		  int equal = 1;
		  for(int i = 0; i < nEl-1; i++)
		    {
		      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
			equal = 0;
		    }
		  
		  if(equal)
		    {
		      // OK, we have a contribution
		      ergo_real x = signb;
		      ergo_real h_pq = h_SO->x[p][q];

		      //printf("1el contrib to a b = %i %i : %22.11f\n", a, b, x * h_pq);

		      resultMatrix[a*nSlaterDets+b] += x * h_pq;
		    }
		} // END FOR q
	    } // END FOR p
	} // END ELSE
    } // END FOR pairIdx
}






typedef struct
{
  int p;
  int q;
  int r;
  int s;
  int sign;
} contrib_debug_struct;



void get_2el_energy_and_gradient(int nSOs,
				 int nEl,
				 int nSlaterDets, 
				 const SlaterDet_struct* SlaterDetList, 
				 int nSlaterDetPairs,
				 const SlaterDet_pair_struct* SlaterDet_pair_list,
				 const int* pairCountList,
				 int noOfTrialVectors,
				 ergo_real* energy_list,
				 ergo_real** coeffListList,
				 const four_idx_SO_struct* g_SO, 
				 ergo_real** resultGradient_list)
{
  for(int k = 0; k < noOfTrialVectors; k++)
    {
      energy_list[k] = 0;
      for(int i = 0; i < nSlaterDets; i++)
	resultGradient_list[k][i] = 0;
    }
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {      
      int a = SlaterDet_pair_list[pairIdx].a;
      int b = SlaterDet_pair_list[pairIdx].b;

      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];
	  
	  int r = SlaterDet_pair_list[pairIdx].SOs_a[1];
	  int s = SlaterDet_pair_list[pairIdx].SOs_b[1];
	  int pos_r = SlaterDet_pair_list[pairIdx].SOs_a_pos[1];
	  int pos_s = SlaterDet_pair_list[pairIdx].SOs_b_pos[1];

	  int sign_a = 1;
	  if(pos_p%2 == 1)
	    sign_a *= -1;
	  if(pos_r%2 == 1)
	    sign_a *= -1;
	  if(r > p)
	    sign_a *= -1;
	  
	  int sign_b = 1;
	  if(pos_q%2 == 1)
	    sign_b *= -1;
	  if(pos_s%2 == 1)
	    sign_b *= -1;
	  if(s > q)
	    sign_a *= -1;

	  // OK, we have a contribution
	  ergo_real x = sign_a * sign_b;
	  x *= 2; // To account for a > b case

	  for(int k = 0; k < noOfTrialVectors; k++)
	    {
	      // pqrs
	      ergo_real g_pqrs = g_SO->x[p][q][r][s];
	      energy_list[k]            += 0.5 * x * g_pqrs * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] += 0.5 * x * g_pqrs * coeffListList[k][b];
	      resultGradient_list[k][b] += 0.5 * x * g_pqrs * coeffListList[k][a];

	      // rspq
	      ergo_real g_rspq = g_SO->x[r][s][p][q];
	      energy_list[k]            += 0.5 * x * g_rspq * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] += 0.5 * x * g_rspq * coeffListList[k][b];
	      resultGradient_list[k][b] += 0.5 * x * g_rspq * coeffListList[k][a];

	      // rqps
	      ergo_real g_rqps = g_SO->x[r][q][p][s];
	      energy_list[k]            -= 0.5 * x * g_rqps * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] -= 0.5 * x * g_rqps * coeffListList[k][b];
	      resultGradient_list[k][b] -= 0.5 * x * g_rqps * coeffListList[k][a];

	      // psrq
	      ergo_real g_psrq = g_SO->x[p][s][r][q];
	      energy_list[k]            -= 0.5 * x * g_psrq * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] -= 0.5 * x * g_psrq * coeffListList[k][b];
	      resultGradient_list[k][b] -= 0.5 * x * g_psrq * coeffListList[k][a];
	    } // END FOR k
	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];

	  for(int rr = 0; rr < nEl; rr++)
	    {
	      int r = SlaterDetList[a].SO_list[rr];
	      if(r == p)
		continue;

	      int sign_a = 1;
	      if(pos_p%2 == 1)
		sign_a *= -1;		
	      if(rr%2 == 1)
		sign_a *= -1;
	      if(r > p)
		sign_a *= -1;
	      
	      for(int ss = 0; ss < nEl; ss++)
		{
		  int s = SlaterDetList[b].SO_list[ss];
		  if(s == q)
		    continue;
		  if(s != r)
		    continue;

		  int sign_b = 1;
		  if(pos_q%2 == 1)
		    sign_b *= -1;
		  if(ss%2 == 1)
		    sign_b *= -1;
		  if(s > q)
		    sign_b *= -1;

		  // OK, we have a contribution
		  ergo_real x = sign_a * sign_b;
		  x *= 2; // To account for a > b case

		  for(int k = 0; k < noOfTrialVectors; k++)
		    {
		      // pqrs
		      ergo_real g_pqrs = g_SO->x[p][q][r][s];
		      energy_list[k]            += 0.5 * x * g_pqrs * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] += 0.5 * x * g_pqrs * coeffListList[k][b];
		      resultGradient_list[k][b] += 0.5 * x * g_pqrs * coeffListList[k][a];

		      // rspq
		      ergo_real g_rspq = g_SO->x[r][s][p][q];
		      energy_list[k]            += 0.5 * x * g_rspq * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] += 0.5 * x * g_rspq * coeffListList[k][b];
		      resultGradient_list[k][b] += 0.5 * x * g_rspq * coeffListList[k][a];

		      // rqps
		      ergo_real g_rqps = g_SO->x[r][q][p][s];
		      energy_list[k]            -= 0.5 * x * g_rqps * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] -= 0.5 * x * g_rqps * coeffListList[k][b];
		      resultGradient_list[k][b] -= 0.5 * x * g_rqps * coeffListList[k][a];

		      // psrq
		      ergo_real g_psrq = g_SO->x[p][s][r][q];
		      energy_list[k]            -= 0.5 * x * g_psrq * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] -= 0.5 * x * g_psrq * coeffListList[k][b];
		      resultGradient_list[k][b] -= 0.5 * x * g_psrq * coeffListList[k][a];
		    } // END FOR k

		} // END FOR ss
	    } // END FOR rr
	}
      else
	{
	  for(int pp = 0; pp < nEl; pp++)
	    {
	      int p = SlaterDetList[a].SO_list[pp];
	      for(int rr = 0; rr < nEl; rr++)
		{
		  int r = SlaterDetList[a].SO_list[rr];
		  int sign_a = 1;
		  if(pp%2 == 1)
		    sign_a *= -1;
		  if(rr%2 == 1)
		    sign_a *= -1;
		  if(pp < rr)
		    sign_a *= -1;		    
		  int qlist[2];
		  qlist[0] = pp;
		  qlist[1] = rr;
		  for(int qqq = 0; qqq < 2; qqq++)
		    {
		      int qq = qlist[qqq];
		      int q = SlaterDetList[b].SO_list[qq];
		      int slist[2];
		      slist[0] = pp;
		      slist[1] = rr;
		      for(int sss = 0; sss < 2; sss++)
			{
			  int ss = slist[sss];
			  int s = SlaterDetList[b].SO_list[ss];
			  if(s == q)
			    continue;
			  int sign_b = 1;
			  if(qq%2 == 1)
			    sign_b *= -1;
			  if(ss%2 == 1)
			    sign_b *= -1;
			  if(qq < ss)
			    sign_b *= -1;

			  // OK, we have a contribution
			  ergo_real x = sign_a * sign_b;
			  ergo_real g_pqrs = g_SO->x[p][q][r][s];

			  for(int k = 0; k < noOfTrialVectors; k++)
			    {
			      energy_list[k]            += 0.5 * x * g_pqrs * coeffListList[k][a] * coeffListList[k][b];
			      resultGradient_list[k][a] += 0.5 * x * g_pqrs * coeffListList[k][b];
			      resultGradient_list[k][b] += 0.5 * x * g_pqrs * coeffListList[k][a];
			    }
			  
			} // END FOR sss
		    } // END FOR qqq
		} // END FOR rr
	    } // END FOR pp

	} // END ELSE

    } // END FOR pair
  
}





void get_2el_contribs_to_mult_or_dmat(int nSOs,
				      int nEl,
				      int nSlaterDets, 
				      const SlaterDet_struct* SlaterDetList, 
				      const SlaterDet_pair_struct* SlaterDetPair,
				      const four_idx_SO_struct* g_SO, 
				      const ergo_real* sourceVector,
				      ergo_real* resultVector, // if result of matrix-vector mult is requested
				      four_idx_SO_struct* result_dmat_2el // if dmat is requested
				      )
{
  int a = SlaterDetPair->a;
  int b = SlaterDetPair->b;

  if(SlaterDetPair->nDiff == 2)
    {
      int p = SlaterDetPair->SOs_a[0];
      int q = SlaterDetPair->SOs_b[0];
      int pos_p = SlaterDetPair->SOs_a_pos[0];
      int pos_q = SlaterDetPair->SOs_b_pos[0];
      
      int r = SlaterDetPair->SOs_a[1];
      int s = SlaterDetPair->SOs_b[1];
      int pos_r = SlaterDetPair->SOs_a_pos[1];
      int pos_s = SlaterDetPair->SOs_b_pos[1];
      
      int sign_a = 1;
      if(pos_p%2 == 1)
	sign_a *= -1;
      if(pos_r%2 == 1)
	sign_a *= -1;
      if(r > p)
	sign_a *= -1;
      
      int sign_b = 1;
      if(pos_q%2 == 1)
	sign_b *= -1;
      if(pos_s%2 == 1)
	sign_b *= -1;
      if(s > q)
	sign_a *= -1;

      // OK, we have a contribution
      ergo_real x = sign_a * sign_b;

      ergo_real g_pqrs = g_SO->x[p][q][r][s];
      ergo_real g_rspq = g_SO->x[r][s][p][q];
      ergo_real g_rqps = g_SO->x[r][q][p][s];
      ergo_real g_psrq = g_SO->x[p][s][r][q];

      ergo_real resultMatrix_ab = 0;
      resultMatrix_ab += 0.5 * x * g_pqrs;
      resultMatrix_ab += 0.5 * x * g_rspq;
      resultMatrix_ab -= 0.5 * x * g_rqps;
      resultMatrix_ab -= 0.5 * x * g_psrq;
      ergo_real resultMatrix_ba = resultMatrix_ab;


      if(resultVector)
	{
	  resultVector[a] += resultMatrix_ab * sourceVector[b];
	  resultVector[b] += resultMatrix_ba * sourceVector[a];
	}
      else
	{
	  ergo_real value = x * sourceVector[a] * sourceVector[b];
	  result_dmat_2el->x[p][q][r][s] += value;
	  result_dmat_2el->x[q][p][s][r] += value;
	  result_dmat_2el->x[r][s][p][q] += value;
	  result_dmat_2el->x[s][r][q][p] += value;
	  result_dmat_2el->x[r][q][p][s] -= value;
	  result_dmat_2el->x[q][r][s][p] -= value;
	  result_dmat_2el->x[p][s][r][q] -= value;
	  result_dmat_2el->x[s][p][q][r] -= value;
	}
    }
  else if(SlaterDetPair->nDiff == 1)
    {
      int p = SlaterDetPair->SOs_a[0];
      int q = SlaterDetPair->SOs_b[0];
      int pos_p = SlaterDetPair->SOs_a_pos[0];
      int pos_q = SlaterDetPair->SOs_b_pos[0];
      
      for(int rr = 0; rr < nEl; rr++)
	{
	  int r = SlaterDetList[a].SO_list[rr];
	  if(r == p)
	    continue;

	  int sign_a = 1;
	  if(pos_p%2 == 1)
	    sign_a *= -1;		
	  if(rr%2 == 1)
	    sign_a *= -1;
	  if(r > p)
	    sign_a *= -1;
	  
	  for(int ss = 0; ss < nEl; ss++)
	    {
	      int s = SlaterDetList[b].SO_list[ss];
	      if(s == q)
		continue;
	      if(s != r)
		continue;

	      int sign_b = 1;
	      if(pos_q%2 == 1)
		sign_b *= -1;
	      if(ss%2 == 1)
		sign_b *= -1;
	      if(s > q)
		sign_b *= -1;

	      // OK, we have a contribution
	      ergo_real x = sign_a * sign_b;

	      ergo_real g_pqrs = g_SO->x[p][q][r][s];
	      ergo_real g_rspq = g_SO->x[r][s][p][q];
	      ergo_real g_rqps = g_SO->x[r][q][p][s];
	      ergo_real g_psrq = g_SO->x[p][s][r][q];

	      ergo_real resultMatrix_ab = 0;
	      resultMatrix_ab += 0.5 * x * g_pqrs;
	      resultMatrix_ab += 0.5 * x * g_rspq;
	      resultMatrix_ab -= 0.5 * x * g_rqps;
	      resultMatrix_ab -= 0.5 * x * g_psrq;
	      ergo_real resultMatrix_ba = resultMatrix_ab;

	      if(resultVector)
		{
		  resultVector[a] += resultMatrix_ab * sourceVector[b];
		  resultVector[b] += resultMatrix_ba * sourceVector[a];
		}
	      else
		{
		  ergo_real value = x * sourceVector[a] * sourceVector[b];
		  result_dmat_2el->x[p][q][r][s] += value;
		  result_dmat_2el->x[q][p][r][s] += value;
		  result_dmat_2el->x[r][s][p][q] += value;
		  result_dmat_2el->x[r][s][q][p] += value;
		  result_dmat_2el->x[r][q][p][s] -= value;
		  result_dmat_2el->x[q][r][s][p] -= value;
		  result_dmat_2el->x[p][s][r][q] -= value;
		  result_dmat_2el->x[s][p][q][r] -= value;
		}
	    } // END FOR ss
	} // END FOR rr
    }
  else
    {
      for(int pp = 0; pp < nEl; pp++)
	{
	  int p = SlaterDetList[a].SO_list[pp];
	  for(int rr = 0; rr < nEl; rr++)
	    {
	      int r = SlaterDetList[a].SO_list[rr];
	      int sign_a = 1;
	      if(pp%2 == 1)
		sign_a *= -1;
	      if(rr%2 == 1)
		sign_a *= -1;
	      if(pp < rr)
		sign_a *= -1;		    
	      int qlist[2];
	      qlist[0] = pp;
	      qlist[1] = rr;
	      for(int qqq = 0; qqq < 2; qqq++)
		{
		  int qq = qlist[qqq];
		  int q = SlaterDetList[b].SO_list[qq];
		  int slist[2];
		  slist[0] = pp;
		  slist[1] = rr;
		  for(int sss = 0; sss < 2; sss++)
		    {
		      int ss = slist[sss];
		      int s = SlaterDetList[b].SO_list[ss];
		      if(s == q)
			continue;
		      int sign_b = 1;
		      if(qq%2 == 1)
			sign_b *= -1;
		      if(ss%2 == 1)
			sign_b *= -1;
		      if(qq < ss)
			sign_b *= -1;

		      // OK, we have a contribution
		      ergo_real x = sign_a * sign_b;
		      ergo_real g_pqrs = g_SO->x[p][q][r][s];

		      ergo_real resultMatrix_ab = 0;
		      resultMatrix_ab += 0.5 * x * g_pqrs;

		      if(resultVector)
			{
			  resultVector[a] += resultMatrix_ab * sourceVector[b];
			}
		      else
			{
			  ergo_real value = x * sourceVector[a] * sourceVector[b];
			  result_dmat_2el->x[p][q][r][s] += value;
			}
		    } // END FOR sss
		} // END FOR qqq
	    } // END FOR rr
	} // END FOR pp
      
    } // END ELSE

}






void get_2el_contribs(int nSOs,
		      int nEl,
		      int nSlaterDets, 
		      const SlaterDet_struct* SlaterDetList, 
		      int nSlaterDetPairs,
		      const SlaterDet_pair_struct* SlaterDet_pair_list,
		      const int* pairCountList,
		      const four_idx_SO_struct* g_SO, 
		      ergo_real* resultMatrix)
{
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {      
      int a = SlaterDet_pair_list[pairIdx].a;
      int b = SlaterDet_pair_list[pairIdx].b;

      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];
	  
	  int r = SlaterDet_pair_list[pairIdx].SOs_a[1];
	  int s = SlaterDet_pair_list[pairIdx].SOs_b[1];
	  int pos_r = SlaterDet_pair_list[pairIdx].SOs_a_pos[1];
	  int pos_s = SlaterDet_pair_list[pairIdx].SOs_b_pos[1];

	  int sign_a = 1;
	  if(pos_p%2 == 1)
	    sign_a *= -1;
	  if(pos_r%2 == 1)
	    sign_a *= -1;
	  if(r > p)
	    sign_a *= -1;
	  
	  int sign_b = 1;
	  if(pos_q%2 == 1)
	    sign_b *= -1;
	  if(pos_s%2 == 1)
	    sign_b *= -1;
	  if(s > q)
	    sign_a *= -1;

	  // OK, we have a contribution
	  ergo_real x = sign_a * sign_b;
	  //x *= 2; // To account for a > b case

	  ergo_real g_pqrs = g_SO->x[p][q][r][s];
	  ergo_real g_rspq = g_SO->x[r][s][p][q];
	  ergo_real g_rqps = g_SO->x[r][q][p][s];
	  ergo_real g_psrq = g_SO->x[p][s][r][q];

	  //printf("2el contrib to a b = %i %i : %22.11f\n", a, b, 0.5 * x * g_pqrs);

	  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_pqrs;
	  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_rspq;
	  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_rqps;
	  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_psrq;

	  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_pqrs;
	  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_rspq;
	  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_rqps;
	  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_psrq;

	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];

	  for(int rr = 0; rr < nEl; rr++)
	    {
	      int r = SlaterDetList[a].SO_list[rr];
	      if(r == p)
		continue;

	      int sign_a = 1;
	      if(pos_p%2 == 1)
		sign_a *= -1;		
	      if(rr%2 == 1)
		sign_a *= -1;
	      if(r > p)
		sign_a *= -1;
	      
	      for(int ss = 0; ss < nEl; ss++)
		{
		  int s = SlaterDetList[b].SO_list[ss];
		  if(s == q)
		    continue;
		  if(s != r)
		    continue;

		  int sign_b = 1;
		  if(pos_q%2 == 1)
		    sign_b *= -1;
		  if(ss%2 == 1)
		    sign_b *= -1;
		  if(s > q)
		    sign_b *= -1;

		  // OK, we have a contribution
		  ergo_real x = sign_a * sign_b;
		  //x *= 2; // To account for a > b case

		  ergo_real g_pqrs = g_SO->x[p][q][r][s];
		  ergo_real g_rspq = g_SO->x[r][s][p][q];
		  ergo_real g_rqps = g_SO->x[r][q][p][s];
		  ergo_real g_psrq = g_SO->x[p][s][r][q];

		  //printf("2el contrib to a b = %i %i : %22.11f\n", a, b, 0.5 * x * g_pqrs);

		  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_pqrs;
		  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_rspq;
		  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_rqps;
		  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_psrq;

		  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_pqrs;
		  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_rspq;
		  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_rqps;
		  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_psrq;

		} // END FOR ss
	    } // END FOR rr
	}
      else
	{
	  for(int pp = 0; pp < nEl; pp++)
	    {
	      int p = SlaterDetList[a].SO_list[pp];
	      for(int rr = 0; rr < nEl; rr++)
		{
		  int r = SlaterDetList[a].SO_list[rr];
		  int sign_a = 1;
		  if(pp%2 == 1)
		    sign_a *= -1;
		  if(rr%2 == 1)
		    sign_a *= -1;
		  if(pp < rr)
		    sign_a *= -1;		    
		  int qlist[2];
		  qlist[0] = pp;
		  qlist[1] = rr;
		  for(int qqq = 0; qqq < 2; qqq++)
		    {
		      int qq = qlist[qqq];
		      int q = SlaterDetList[b].SO_list[qq];
		      int slist[2];
		      slist[0] = pp;
		      slist[1] = rr;
		      for(int sss = 0; sss < 2; sss++)
			{
			  int ss = slist[sss];
			  int s = SlaterDetList[b].SO_list[ss];
			  if(s == q)
			    continue;
			  int sign_b = 1;
			  if(qq%2 == 1)
			    sign_b *= -1;
			  if(ss%2 == 1)
			    sign_b *= -1;
			  if(qq < ss)
			    sign_b *= -1;

			  // OK, we have a contribution
			  ergo_real x = sign_a * sign_b;
			  ergo_real g_pqrs = g_SO->x[p][q][r][s];

			  //printf("2el contrib to a b = %i %i : %22.11f\n", a, b, 0.5 * x * g_pqrs);

			  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_pqrs;
			  //resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_pqrs;

			} // END FOR sss
		    } // END FOR qqq
		} // END FOR rr
	    } // END FOR pp

	} // END ELSE

    } // END FOR pair
  
}








int get_1e_density_matrix(int nSOs,
			  int nEl,
			  two_idx_SO_struct* D,
			  int nSlaterDets, 
			  const SlaterDet_struct* SlaterDetList, 
			  const ergo_real* coeffList)
{
  int p, q;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      {
	// Compute element pq
	ergo_real sum = 0;
	int a, b;
	for(a = 0; a < nSlaterDets; a++)
	  for(b = 0; b < nSlaterDets; b++)
	    {
	      int i, count, savedCount;
	      int sign = 1;
	      // annihilate p from a
	      char SO_list_mod_a_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[a].SO_list[i] == p)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		      count++;
		      if(savedCount < 0)
			sign *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;

	      // annihilate q from b
	      char SO_list_mod_b_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[b].SO_list[i] == q)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		      count++;
		      if(savedCount < 0)
			sign *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;

	      int equal = 1;
	      for(i = 0; i < nEl-1; i++)
		{
		  if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		    equal = 0;
		}
	      
	      if(equal)
		sum += sign * coeffList[a] * coeffList[b];
	    } // END FOR a b
	D->x[p][q] = sum;
      } // END FOR p q

  // Verify symmetry
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      {
	if(std::fabs(D->x[p][q] - D->x[q][p]) > 1e-9)
	  printf("ERROR 11!\n");
      }
  ergo_real sum = 0;
  for(p = 0; p < nSOs; p++)
    sum += D->x[p][p];
#if 0
  if(std::fabs(sum - nEl) > 1e-9)
    printf("Tr ERROR!\n");
#endif
  //printf("Tr(D) = %22.11f\n", sum);  

  return 0;
}


int get_2e_density_matrix(int nSOs,
			  int nEl,
			  four_idx_SO_struct* d,
			  int nSlaterDets, 
			  const SlaterDet_struct* SlaterDetList, 
			  const ergo_real* coeffList)
{
  int p, q, r, s;
  int a, b;
  for(p = 0; p < nSOs; p++)
    for(r = 0; r < nSOs; r++)
      for(a = 0; a < nSlaterDets; a++)
	{
	  int i, count, savedCount;
	  int sign1 = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    sign1 *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
	  // annihilate r from a
	  char SO_list_mod_a_2[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(i = 0; i < nEl-1; i++)
	    {
	      if(SO_list_mod_a_1[i] == r)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_2[count] = SO_list_mod_a_1[i];
		  count++;
		  if(savedCount < 0)
		    sign1 *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
	  
	  for(q = 0; q < nSOs; q++)
	    for(s = 0; s < nSOs; s++)
	      for(b = 0; b < nSlaterDets; b++)
		{
		  int sign = sign1;
		  
		  // annihilate q from b
		  char SO_list_mod_b_1[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(i = 0; i < nEl; i++)
		    {
		      if(SlaterDetList[b].SO_list[i] == q)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
			  count++;
			  if(savedCount < 0)
			    sign *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  // annihilate s from b
		  char SO_list_mod_b_2[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(i = 0; i < nEl-1; i++)
		    {
		      if(SO_list_mod_b_1[i] == s)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_2[count] = SO_list_mod_b_1[i];
			  count++;
			  if(savedCount < 0)
			    sign *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  
		  int equal = 1;
		  for(i = 0; i < nEl-2; i++)
		    {
		      if(SO_list_mod_a_2[i] != SO_list_mod_b_2[i])
			equal = 0;
		    }
		  
		  if(equal)
		    d->x[p][q][r][s] += sign * coeffList[a] * coeffList[b];
		} // END FOR q s b
	} // END FOR p r a

  //printf("d done\n");

  // Verify symmetry
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      for(r = 0; r < nSOs; r++)
	for(s = 0; s < nSOs; s++)
	  {
	    if(std::fabs(d->x[p][q][r][s] - d->x[r][s][p][q]) > 1e-9)
	      printf("ERROR 1!");
	    if(std::fabs(d->x[p][q][r][s] + d->x[r][q][p][s]) > 1e-9)
	      printf("ERROR 2!");
	    if(std::fabs(d->x[p][q][r][s] + d->x[p][s][r][q]) > 1e-9)
	      printf("ERROR 3!");
	    if(std::fabs(d->x[p][q][p][s]) > 1e-9)
	      printf("ERROR 4!");
	    if(std::fabs(d->x[p][q][r][q]) > 1e-9)
	      printf("ERROR 5!");
	    if(std::fabs(d->x[p][q][p][q]) > 1e-9)
	      printf("ERROR 6!");
	  }
  
#if 0
  p = 1; q = 2; r = 6; s = 7;
  printf("d[p][q][r][s] = %22.11f\n", d->x[p][q][r][s]);
  printf("d[r][q][p][s] = %22.11f\n", d->x[r][q][p][s]);
  printf("d[p][s][r][q] = %22.11f\n", d->x[p][s][r][q]);
  printf("d[r][s][p][q] = %22.11f\n", d->x[r][s][p][q]);
  printf("\n");
#endif

  return 0;
}






ergo_real get_CI_energy(int nSOs,
			int nEl,
			const four_idx_SO_struct* g_SO, 
			const two_idx_SO_struct* h_SO, 
			int nSlaterDets, 
			const SlaterDet_struct* SlaterDetList, 
			const ergo_real* coeffList,
			ergo_real nuclearEnergy)
{
  two_idx_SO_struct* D = new two_idx_SO_struct;
  four_idx_SO_struct* d = new four_idx_SO_struct;

  get_1e_density_matrix(nSOs,
			nEl,
			D, 
			nSlaterDets, 
			SlaterDetList, 
			coeffList);
  
  get_2e_density_matrix(nSOs,
			nEl,
			d, 
			nSlaterDets, 
			SlaterDetList, 
			coeffList);

  int p, q, r, s;
  
  ergo_real energy_1el = 0;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      energy_1el += D->x[p][q] * h_SO->x[p][q];
  
  ergo_real energy_2el = 0;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      for(r = 0; r < nSOs; r++)
	for(s = 0; s < nSOs; s++)
	  energy_2el += 0.5 * d->x[p][q][r][s] * g_SO->x[p][q][r][s];

  delete D;
  delete d;
  
  return energy_1el + energy_2el + nuclearEnergy;
}



void get_CI_energy_and_gradient(int nSOs,
				int nEl,
				const four_idx_SO_struct* g_SO, 
				const two_idx_SO_struct* h_SO, 
				int nSlaterDets, 
				const SlaterDet_struct* SlaterDetList, 
				int nSlaterDetPairs,
				const SlaterDet_pair_struct* SlaterDet_pair_list,
				const int* pairCountList,
				int noOfTrialVectors,
				ergo_real* energyList,
				ergo_real** coeffListList,
				ergo_real** gradientList,
				ergo_real nuclearEnergy)
{
  ergo_real* gradient_1el_list[88];
  for(int i = 0; i < noOfTrialVectors; i++)
    gradient_1el_list[i] = new ergo_real[nSlaterDets];
  ergo_real energy_1el_list[88];
  get_1el_energy_and_gradient(nSOs,
			      nEl,
			      nSlaterDets, 
			      SlaterDetList, 
			      nSlaterDetPairs,
			      SlaterDet_pair_list,
			      pairCountList,
			      noOfTrialVectors,
			      energy_1el_list,
			      coeffListList,
			      h_SO,
			      gradient_1el_list);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "1-electron part done");

  //printf("energy_1el_list[0] = %22.11f\n", energy_1el_list[0]);

  ergo_real* gradient_2el_list[88];
  for(int i = 0; i < noOfTrialVectors; i++)
    gradient_2el_list[i] = new ergo_real[nSlaterDets];
  ergo_real energy_2el_list[88];
  get_2el_energy_and_gradient(nSOs,
			      nEl,
			      nSlaterDets,
			      SlaterDetList,
			      nSlaterDetPairs,
			      SlaterDet_pair_list,
			      pairCountList,
			      noOfTrialVectors,
			      energy_2el_list,
			      coeffListList,
			      g_SO,
			      gradient_2el_list);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "2-electron part done");
  
  // Take into account that we want the derivative when the sum of squares is fixed.
  for(int k = 0; k < noOfTrialVectors; k++)
    for(int i = 0; i < nSlaterDets; i++)
      gradientList[k][i] = gradient_1el_list[k][i] + gradient_2el_list[k][i] - 2 * (energy_1el_list[k] + energy_2el_list[k]) * coeffListList[k][i];
    
  for(int k = 0; k < noOfTrialVectors; k++)
    {
      delete gradient_1el_list[k];
      delete gradient_2el_list[k];
    }
  
  for(int k = 0; k < noOfTrialVectors; k++)
    energyList[k] = energy_1el_list[k] + energy_2el_list[k] + nuclearEnergy;
}



void mult_by_CI_matrix(int nSOs,
		       int nEl,
		       const four_idx_SO_struct* g_SO, 
		       const two_idx_SO_struct* h_SO, 
		       int nSlaterDets, 
		       const SlaterDet_struct* SlaterDetList, 
		       int nSlaterDetPairs,
		       const SlaterDet_pair_struct* SlaterDet_pair_list,
		       const int* pairCountList,
		       const ergo_real* sourceVector,
		       ergo_real* resultVector,
		       ergo_real shift)
{
  memset(resultVector, 0, nSlaterDets*sizeof(ergo_real));
  for(int i = 0; i < nSlaterDetPairs; i++)
    {
      get_1el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets, 
				       SlaterDetList, 
				       &SlaterDet_pair_list[i],
				       h_SO,
				       sourceVector,
				       resultVector,
				       NULL);
      get_2el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets,
				       SlaterDetList,
				       &SlaterDet_pair_list[i],
				       g_SO,
				       sourceVector,
				       resultVector,
				       NULL);
    } // END FOR i
  // Apply shift
  for(int i = 0; i < nSlaterDets; i++)
    resultVector[i] -= shift * sourceVector[i];
}



void get_CI_matrix(int nSOs,
		   int nEl,
		   const four_idx_SO_struct* g_SO, 
		   const two_idx_SO_struct* h_SO, 
		   int nSlaterDets, 
		   const SlaterDet_struct* SlaterDetList, 
		   int nSlaterDetPairs,
		   const SlaterDet_pair_struct* SlaterDet_pair_list,
		   const int* pairCountList,
		   ergo_real* resultMatrix)
{
  memset(resultMatrix, 0, nSlaterDets*nSlaterDets*sizeof(ergo_real));
  get_1el_contribs(nSOs,
		   nEl,
		   nSlaterDets, 
		   SlaterDetList, 
		   nSlaterDetPairs,
		   SlaterDet_pair_list,
		   pairCountList,
		   h_SO,
		   resultMatrix);
  get_2el_contribs(nSOs,
		   nEl,
		   nSlaterDets,
		   SlaterDetList,
		   nSlaterDetPairs,
		   SlaterDet_pair_list,
		   pairCountList,
		   g_SO,
		   resultMatrix);
}



int get_combinations(SlaterDet_struct* SlaterDetList, 
		     int nEl, 
		     int nSOs)
{
  if(nEl == 0)
    return 1;

  int list[MAX_SOS];
  memset(list, 0, MAX_SOS * sizeof(int));

  int i;
  for(i = 0; i < nEl; i++)
    list[i] = i;

  int finishedFlag = 0;
  int count = 0;
  while(finishedFlag == 0)
    {
      if(SlaterDetList)
	{
	  // Create det for current configuration
	  for(i = 0; i < nEl; i++)
	    SlaterDetList[count].SO_list[i] = list[i];
	}
      count++;

      // Move to next configuration
      int breakFlag = 0;
      int level = nEl - 1;      
      while(breakFlag == 0)
	{
	  // Check if value at current level can be increased
	  if(list[level] < nSOs-(nEl-level))
	    {
	      // OK, increase value
	      list[level]++;
	      // Set all values at higher levels as low as possible
	      for(i = level+1; i < nEl; i++)
		list[i] = list[i-1] + 1;
	      breakFlag = 1;
	    }
	  else
	    {
	      // Go to next level
	      level--;
	    }
	  if(level < 0)
	    {
	      breakFlag = 1;
	      finishedFlag = 1;
	    }
	}
    }

  return count;  
}


int get_FCI_Slater_dets_alpha_beta(SlaterDet_struct* SlaterDetList, 
				   int nEl_a, 
				   int nEl_b, 
				   int nSOs)
{
  int nCombs_a = get_combinations(NULL, nEl_a, nSOs/2);  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "get_FCI_Slater_dets_alpha_beta, nCombs_a = %9i", nCombs_a);
  SlaterDet_struct* SlaterDetList_a = new SlaterDet_struct[nCombs_a];
  get_combinations(SlaterDetList_a, nEl_a, nSOs/2);

  int nCombs_b = get_combinations(NULL, nEl_b, nSOs/2);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "get_FCI_Slater_dets_alpha_beta, nCombs_b = %9i", nCombs_b);
  SlaterDet_struct* SlaterDetList_b = new SlaterDet_struct[nCombs_b];
  get_combinations(SlaterDetList_b, nEl_b, nSOs/2);

  int count = 0;
  int ia, ib;
  for(ia = 0; ia < nCombs_a; ia++)
    for(ib = 0; ib < nCombs_b; ib++)
      {
#if 0
	if(!((SlaterDetList_a[ia].SO_list[0] == 0 && nSOs / 2 + SlaterDetList_b[ib].SO_list[0] == 4) || (SlaterDetList_a[ia].SO_list[0] == 1 && nSOs / 2 + SlaterDetList_b[ib].SO_list[0] == 4)))
	  continue;
#endif
	if(SlaterDetList)
	  {
	    for(int i = 0; i < nEl_a; i++)
	      SlaterDetList[count].SO_list[i] = SlaterDetList_a[ia].SO_list[i];
	    for(int i = 0; i < nEl_b; i++)
	      SlaterDetList[count].SO_list[nEl_a+i] = nSOs / 2 + SlaterDetList_b[ib].SO_list[i];
	    SlaterDetList[count].startIndex = -1;
	    SlaterDetList[count].count = 0;
#if 0
	printf("new det: ");
	for(int i = 0; i < nEl_a+nEl_b; i++)
	  printf("%3i ", SlaterDetList[count].SO_list[i]);
	printf("\n");
#endif
	  }

	count++;
      }

  delete [] SlaterDetList_a;
  delete [] SlaterDetList_b;

  return count;
}



int get_FCI_Slater_dets_all(SlaterDet_struct* SlaterDetList, 
			    int nElTot, 
			    int nSOs)
{
  return get_combinations(SlaterDetList, nElTot, nSOs);
}



static ergo_real rand_m1_to_1()
{
  ergo_real randomNumber = (ergo_real)rand() / RAND_MAX;
  // Now randomNumber is between 0 and 1
  randomNumber *= 2;
  // Now randomNumber is between 0 and 2
  randomNumber -= 1;
  // Now randomNumber is between -1 and 1
  return randomNumber;
}



typedef struct
{
  int aDiffList[2];
  int bDiffList[2];
  int aDiffPosList[2];
  int bDiffPosList[2];
} pair_status_struct;





void get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(SlaterDet_struct** groupList, 
							    int nEl,
							    int level,
							    int groupIdx1,
							    int groupIdx2,
							    int ia,
							    int ib,
							    int aDiffCount,
							    int bDiffCount,
							    pair_status_struct* status,
							    int nSOs,
							    int nSlaterDets,
							    const SlaterDet_struct* SlaterDetList,
							    const four_idx_SO_struct* g_SO, 
							    const two_idx_SO_struct* h_SO, 
							    const ergo_real* sourceVector,
							    ergo_real* resultVector,
							    two_idx_SO_struct* result_dmat_1el,
							    four_idx_SO_struct* result_dmat_2el
							    )
{
  if(groupIdx1 > groupIdx2)
    return;

  if(ia != level-1 && ib != level-1 && level != 0)
    throw "Error in CI: (ia != level && ib != level)";

  // ia or ib or both point to an electron at this level.
  const SlaterDet_struct* group_1 = &groupList[level][groupIdx1];
  const SlaterDet_struct* group_2 = &groupList[level][groupIdx2];

  if(level == nEl)
    {
      // final level
      while(ia < level || ib < level)
	{
	  if(ia == nEl)
	    {
	      // Only b left, must be diff
	      if(bDiffCount == 2)
		return;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	      continue;
	    }
	  if(ib == nEl)
	    {
	      // Only a left, must be diff
	      if(aDiffCount == 2)
		return;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	      continue;
	    }
	  while(ia < level && ib < level)
	    {
	      if(group_1->SO_list[ia] == group_2->SO_list[ib])
		{
		  ia++;
		  ib++;
		}
	      else if(group_1->SO_list[ia] > group_2->SO_list[ib])
		{
		  if(bDiffCount == 2)
		    return;
		  status->bDiffList[bDiffCount] = group_2->SO_list[ib];
		  status->bDiffPosList[bDiffCount] = ib;
		  bDiffCount++;
		  ib++;
		}
	      else
		{	
		  if(aDiffCount == 2)
		    return;
		  status->aDiffList[aDiffCount] = group_1->SO_list[ia];
		  status->aDiffPosList[aDiffCount] = ia;
		  aDiffCount++;
		  ia++;
		}
	    } // END WHILE
	} // END WHILE

      // Do contribution to matrix-vector multiplication.
      SlaterDet_pair_struct SlaterDetPair;
      SlaterDetPair.a = groupIdx1;
      SlaterDetPair.b = groupIdx2;
      SlaterDetPair.nDiff = aDiffCount;
      for(int i = 0; i < aDiffCount; i++)
	{
	  SlaterDetPair.SOs_a[i] = status->aDiffList[i];
	  SlaterDetPair.SOs_b[i] = status->bDiffList[i];
	  SlaterDetPair.SOs_a_pos[i] = status->aDiffPosList[i];
	  SlaterDetPair.SOs_b_pos[i] = status->bDiffPosList[i];
	}
      get_1el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets, 
				       SlaterDetList, 
				       &SlaterDetPair,
				       h_SO,
				       sourceVector,
				       resultVector,
				       result_dmat_1el);
      get_2el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets,
				       SlaterDetList,
				       &SlaterDetPair,
				       g_SO,
				       sourceVector,
				       resultVector,
				       result_dmat_2el);
      
      return;
    } // END IF final level

  if(level > 0)
    {
      while(ia < level && ib < level)
	{
	  if(group_1->SO_list[ia] == group_2->SO_list[ib])
	    {
	      ia++;
	      ib++;
	    }
	  else if(group_1->SO_list[ia] > group_2->SO_list[ib])
	    {
	      if(bDiffCount == 2)
		return;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	    }
	  else
	    {
	      if(aDiffCount == 2)
		return;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	    }
	} // END WHILE
    } // END IF level > 0
 
  // No, we could not skip. Go to next level.
  for(int i = 0; i < group_1->count; i++)
    for(int j = 0; j < group_2->count; j++)
      {
	get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							       nEl,
							       level+1,
							       group_1->startIndex + i,
							       group_2->startIndex + j,
							       ia,
							       ib,
							       aDiffCount,
							       bDiffCount,
							       status,
							       nSOs,
							       nSlaterDets,
							       SlaterDetList,
							       g_SO, 
							       h_SO, 
							       sourceVector,
							       resultVector,
							       result_dmat_1el,
							       result_dmat_2el
							       );
      } // END FOR i j
  return;
}





void mult_by_CI_matrix_direct(int nSOs,
			      int nEl,
			      const four_idx_SO_struct* g_SO, 
			      const two_idx_SO_struct* h_SO, 
			      int nSlaterDets, 
			      const SlaterDet_struct* SlaterDetList, 
			      SlaterDet_struct** groupList, 
			      const ergo_real* sourceVector,
			      ergo_real* resultVector,
			      ergo_real shift)
{
  memset(resultVector, 0, nSlaterDets*sizeof(ergo_real));

  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							 nEl,
							 0,
							 0,
							 0,
							 0, 0, 0, 0,
							 &status,
							 nSOs,
							 nSlaterDets,
							 SlaterDetList,
							 g_SO, 
							 h_SO, 
							 sourceVector,
							 resultVector,
							 NULL,
							 NULL
							 );
  
  // Apply shift
  for(int i = 0; i < nSlaterDets; i++)
    resultVector[i] -= shift * sourceVector[i];
}



static ergo_real get_SlaterDet_energy(int nSOs,
				      int nEl,
				      const four_idx_SO_struct* g_SO, 
				      const two_idx_SO_struct* h_SO, 
				      const SlaterDet_struct* SlaterDet)
{
  // Create "Slater det pair" consisting of this determinant paired with itself.
  SlaterDet_pair_struct SlaterDetPair;
  SlaterDetPair.a = 0;
  SlaterDetPair.b = 0;
  SlaterDetPair.nDiff = 0;
  const int nSlaterDets = 1;
  ergo_real sourceVector[1];
  ergo_real resultVector[1];
  sourceVector[0] = 1.0;
  resultVector[0] = 0.0;
  get_1el_contribs_to_mult_or_dmat(nSOs,
				   nEl,
				   nSlaterDets, 
				   SlaterDet, 
				   &SlaterDetPair,
				   h_SO,
				   sourceVector,
				   resultVector,
				   NULL);
  get_2el_contribs_to_mult_or_dmat(nSOs,
				   nEl,
				   nSlaterDets,
				   SlaterDet,
				   &SlaterDetPair,
				   g_SO,
				   sourceVector,
				   resultVector,
				   NULL);
  return resultVector[0];
}



int get_relevant_SlaterDet_pairs_recursive_2(SlaterDet_struct** groupList, 
					     int nEl,
					     SlaterDet_pair_struct* resultList,
					     int level,
					     int groupIdx1,
					     int groupIdx2,
					     int ia,
					     int ib,
					     int aDiffCount,
					     int bDiffCount,
					     pair_status_struct* status)
{
  if(groupIdx1 > groupIdx2)
    return 0;

  if(ia != level-1 && ib != level-1 && level != 0)
    throw "Error in CI: (ia != level && ib != level)";

  // ia or ib or both point to an electron at this level.
  SlaterDet_struct* group_1 = &groupList[level][groupIdx1];
  SlaterDet_struct* group_2 = &groupList[level][groupIdx2];

  if(level == nEl)
    {
      // final level
      while(ia < level || ib < level)
	{
	  if(ia == nEl)
	    {
	      // Only b left, must be diff
	      if(bDiffCount == 2)
		return 0;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	      continue;
	    }
	  if(ib == nEl)
	    {
	      // Only a left, must be diff
	      if(aDiffCount == 2)
		return 0;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	      continue;
	    }
	  while(ia < level && ib < level)
	    {
	      if(group_1->SO_list[ia] == group_2->SO_list[ib])
		{
		  ia++;
		  ib++;
		}
	      else if(group_1->SO_list[ia] > group_2->SO_list[ib])
		{
		  if(bDiffCount == 2)
		    return 0;
		  status->bDiffList[bDiffCount] = group_2->SO_list[ib];
		  status->bDiffPosList[bDiffCount] = ib;
		  bDiffCount++;
		  ib++;
		}
	      else
		{	
		  if(aDiffCount == 2)
		    return 0;
		  status->aDiffList[aDiffCount] = group_1->SO_list[ia];
		  status->aDiffPosList[aDiffCount] = ia;
		  aDiffCount++;
		  ia++;
		}
	    } // END WHILE
	} // END WHILE
      
      if(resultList)
	{
	  resultList[0].a = groupIdx1;
	  resultList[0].b = groupIdx2;
	  resultList[0].nDiff = aDiffCount;
	  for(int i = 0; i < aDiffCount; i++)
	    {
	      resultList[0].SOs_a[i] = status->aDiffList[i];
	      resultList[0].SOs_b[i] = status->bDiffList[i];
	      resultList[0].SOs_a_pos[i] = status->aDiffPosList[i];
	      resultList[0].SOs_b_pos[i] = status->bDiffPosList[i];
	    }
	}
      return 1;
    } // END IF final level

  if(level > 0)
    {
      while(ia < level && ib < level)
	{
	  if(group_1->SO_list[ia] == group_2->SO_list[ib])
	    {
	      ia++;
	      ib++;
	    }
	  else if(group_1->SO_list[ia] > group_2->SO_list[ib])
	    {
	      if(bDiffCount == 2)
		return 0;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	    }
	  else
	    {
	      if(aDiffCount == 2)
		return 0;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	    }
	} // END WHILE
    } // END IF level > 0
 
  // No, we could not skip. Go to next level.
  int count = 0;
  for(int i = 0; i < group_1->count; i++)
    for(int j = 0; j < group_2->count; j++)
      {
	SlaterDet_pair_struct* resultListPtr = NULL;
	if(resultList)
	  resultListPtr = &resultList[count];
	int currCount = get_relevant_SlaterDet_pairs_recursive_2(groupList, 
								 nEl,
								 resultListPtr,
								 level+1,
								 group_1->startIndex + i,
								 group_2->startIndex + j,
								 ia,
								 ib,
								 aDiffCount,
								 bDiffCount,
								 status);
	count += currCount;
      } // END FOR i j
  return count;
}






int get_relevant_SlaterDet_pairs_recursive(int nSlaterDets, 
					   SlaterDet_struct* SlaterDetList, 
					   SlaterDet_struct** groupList, 
					   int nEl,
					   SlaterDet_pair_struct* resultList,
					   int level,
					   int groupIdx1,
					   int groupIdx2)
{
  //printf("get_relevant_SlaterDet_pairs_recursive, level = %i\n", level);

  // Check if this pair of groups can be skipped.
  SlaterDet_struct* group_1 = &groupList[level][groupIdx1];
  SlaterDet_struct* group_2 = &groupList[level][groupIdx2];
  if(level < nEl)
    {
      int nSame = 0;
      int ia = 0;
      int ib = 0;
      int aDiffCount = 0;
      int bDiffCount = 0;
      //int aDiffList[MAX_SOS];
      //int bDiffList[MAX_SOS];
      //int aDiffPosList[MAX_SOS];
      //int bDiffPosList[MAX_SOS];
      while(ia < level && ib < level)
	{
	  if(group_1->SO_list[ia] == group_2->SO_list[ib])
	    {
	      nSame++;
	      ia++;
	      ib++;
	    }
	  else
	    {
	      if(group_1->SO_list[ia] > group_2->SO_list[ib])
		{
		  //bDiffList[bDiffCount] = group_2->SO_list[ib];
		  //bDiffPosList[bDiffCount] = ib;
		  bDiffCount++;
		  if(bDiffCount > 2)
		    break;
		  ib++;
		}
	      else
		{
		  //aDiffList[aDiffCount] = group_1->SO_list[ia];
		  //aDiffPosList[aDiffCount] = ia;
		  aDiffCount++;
		  if(aDiffCount > 2)
		    break;
		  ia++;
		}
	    }
	}
      if(aDiffCount > 2 || bDiffCount > 2)
	{
	  //printf("aDiffCount = %i  bDiffCount = %i  level = %i\n", aDiffCount, bDiffCount, level);
	  return 0;
	}
 
      // No, we could not skip. Go to next level.
      int count = 0;
      for(int i = 0; i < group_1->count; i++)
	for(int j = 0; j < group_2->count; j++)
	  {
	    SlaterDet_pair_struct* resultListPtr = NULL;
	    if(resultList)
	      resultListPtr = &resultList[count];
	    int currCount = get_relevant_SlaterDet_pairs_recursive(nSlaterDets, 
								   SlaterDetList, 
								   groupList, 
								   nEl,
								   resultListPtr,
								   level+1,
								   group_1->startIndex + i,
								   group_2->startIndex + j);
	    count += currCount;
	  } // END FOR i j
      return count;
    } // END IF level < nEl
  
  // We are at lowest level.

  //printf("lowest level!\n");

  int a = groupIdx1;
  int b = groupIdx2;

  if(b < a)
    return 0;
  
  int nSame = 0;
  int ia = 0;
  int ib = 0;
  int aDiffCount = 0;
  int bDiffCount = 0;
  int aDiffList[MAX_SOS];
  int bDiffList[MAX_SOS];
  int aDiffPosList[MAX_SOS];
  int bDiffPosList[MAX_SOS];
  //printf("ia ib loop starting.\n");
  while(ia < nEl || ib < nEl)
    {
      if(ia == nEl)
	{
	  // Only b left, must be diff
	  bDiffList[bDiffCount] = SlaterDetList[b].SO_list[ib];
	  bDiffPosList[bDiffCount] = ib;
	  bDiffCount++;
	  ib++;
	  continue;
	}
      if(ib == nEl)
	{
	  // Only a left, must be diff
	  aDiffList[aDiffCount] = SlaterDetList[a].SO_list[ia];
	  aDiffPosList[aDiffCount] = ia;
	  aDiffCount++;
	  ia++;
	  continue;
	}
      if(SlaterDetList[a].SO_list[ia] == SlaterDetList[b].SO_list[ib])
	{
	  nSame++;
	  ia++;
	  ib++;
	}
      else
	{
	  if(SlaterDetList[a].SO_list[ia] > SlaterDetList[b].SO_list[ib])
	    {
	      bDiffList[bDiffCount] = SlaterDetList[b].SO_list[ib];
	      bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      if(bDiffCount > 2)
		break;
	      ib++;
	    }
	  else
	    {
	      aDiffList[aDiffCount] = SlaterDetList[a].SO_list[ia];
	      aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      if(aDiffCount > 2)
		break;
	      ia++;
	    }
	}
    }
  //printf("while ia ib loop done.\n");
#if 0
  if(aDiffCount != bDiffCount)
    {
      printf("ERROR: (aDiffCount != bDiffCount)\n");
      exit(0);
    }
  if(nSame + aDiffCount != nEl)
    {
      printf("ERROR: (nSame + aDiffCount != nEl)\n");
      exit(0);
    }
#endif
  if(nSame >= nEl - 2)
    {
      if(resultList)
	{
	  int pairCount = 0;
	  resultList[pairCount].a = a;
	  resultList[pairCount].b = b;
	  resultList[pairCount].nDiff = aDiffCount;
	  for(int i = 0; i < 2; i++)
	    {
	      resultList[pairCount].SOs_a[i] = -1;
	      resultList[pairCount].SOs_b[i] = -1;
	      resultList[pairCount].SOs_a_pos[i] = -1;
	      resultList[pairCount].SOs_b_pos[i] = -1;
	    }
	  for(int i = 0; i < aDiffCount; i++)
	    {
	      resultList[pairCount].SOs_a[i] = aDiffList[i];
	      resultList[pairCount].SOs_b[i] = bDiffList[i];
	      resultList[pairCount].SOs_a_pos[i] = aDiffPosList[i];
	      resultList[pairCount].SOs_b_pos[i] = bDiffPosList[i];
	    }
	}
      return 1;
    }
  else
    {
      return 0;
    }
}



int get_relevant_SlaterDet_pairs(int nSlaterDets, 
				 SlaterDet_struct* SlaterDetList, 
				 SlaterDet_struct** groupList, 
				 int nEl,
				 SlaterDet_pair_struct* resultList)
{
  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  return get_relevant_SlaterDet_pairs_recursive_2(groupList, 
						  nEl,
						  resultList,
						  0,
						  0,
						  0,
						  0, 0, 0, 0,
						  &status);
}




ergo_real get_eigs(int n, ergo_real* M, ergo_real* bestVector, ergo_real* eigValListResult)
{
  int lwork = 3*n*n;
  ergo_real* work = new ergo_real[lwork];
  ergo_real* eigvalList = new ergo_real[n];
  ergo_real* A = new ergo_real[n*n];
  memcpy(A, M, n*n*sizeof(ergo_real));
  int info = 0;
  mat::syev("V", "U", &n, A,
            &n, eigvalList, work, &lwork, 
            &info);
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "error in syev, info = %i", info);
      exit(0);
    }

#if 0
  for(int i = 0; i < n; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_CI, "get_eigs eigenValue %5i: %22.11f", i, (double)eigvalList[i]);
#endif

  if(bestVector)
    {
      for(int i = 0; i < n; i++)
	bestVector[i] = A[0*n+i];
    }

  if(eigValListResult)
    {
      for(int i = 0; i < n; i++)
	eigValListResult[i] = eigvalList[i];
    }

  ergo_real result = eigvalList[0];

  delete [] work;
  delete [] eigvalList;
  delete [] A;

  return result;
}






int get_Lowdin_orbitals(int n, const ergo_real* S, ergo_real* MOs)
{
  int lwork = 3*n*n;
  ergo_real* work = new ergo_real[lwork];
  ergo_real* eigvalList = new ergo_real[n];
  ergo_real* A = new ergo_real[n*n];
  memcpy(A, S, n*n*sizeof(ergo_real));
  int info = 0;

  mat::syev("V", "U", &n, A,
	    &n, eigvalList, work, &lwork, 
	    &info);
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "error in syev, info = %i", info);
      exit(0);
    }
  
  for(int i = 0; i < n; i++)
    {
      assert(eigvalList[i] > 0);
    }

  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < n; k++)
	  sum += A[k*n+i] * A[k*n+j] * (1.0 / std::sqrt(eigvalList[k]));
	MOs[i*n+j] = sum;
      }

  return 0;
}




typedef ergo_real* ergp_real_ptr;



ergo_real do_lanczos_method_direct(int n,
				   ergo_real* v,
				   int nSOs,
				   int nEl,
				   const four_idx_SO_struct* g_SO, 
				   const two_idx_SO_struct* h_SO, 
				   int nSlaterDets, 
				   const SlaterDet_struct* SlaterDetList, 
				   SlaterDet_struct** groupList,
				   int maxIterations_in,
				   ergo_real nucRepEnergy,
				   ergo_real shift)
{
  int maxIterations = maxIterations_in;
  if(maxIterations > n)
    maxIterations = n;
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "do_lanczos_method_direct, n = %9i, maxIterations = %3i, shift = %8.3f",
	    n, maxIterations, (double)shift);
  ergo_real** q = new ergp_real_ptr[n];
  q[0] = new ergo_real[n];
  for(int i = 0; i < n; i++)
    q[0][i] = 0;
  q[1] = new ergo_real[n];
  for(int i = 0; i < n; i++)
    q[1][i] = v[i];
  normalize_vector(n, q[1]);
  std::vector<ergo_real> z(n);
  std::vector<ergo_real> alpha(n);
  std::vector<ergo_real> beta(n);
  beta[0] = 0;
  ergo_real currEig = 0;
  ergo_real curr_E = 0;
  for(int j = 1; j < maxIterations; j++)
    {
      // Do matrix-vector multiplication
      mult_by_CI_matrix_direct(nSOs,
			       nEl,
			       g_SO, 
			       h_SO, 
			       nSlaterDets, 
			       SlaterDetList, 
			       groupList, 
			       q[j],
			       &z[0],
			       shift);
      // OK, matrix-vector multiplication done
      alpha[j] = 0;
      for(int i = 0; i < n; i++)
	alpha[j] += q[j][i] * z[i];
      for(int i = 0; i < n; i++)
	z[i] = z[i] - alpha[j] * q[j][i] - beta[j-1] * q[j-1][i];
      beta[j] = get_vector_norm(n, &z[0]);

      do_output(LOG_CAT_INFO, LOG_AREA_CI, "do_lanczos_method, j = %5i, alpha = %22.11f, beta = %22.11f", j, (double)alpha[j], (double)beta[j]);

      if(beta[j] < 1e-5)
	break;
      q[j+1] = new ergo_real[n];
      for(int i = 0; i < n; i++)
	q[j+1][i] = z[i] / beta[j];

      ergo_real* T = new ergo_real[j*j];
      for(int i = 0; i < j*j; i++)
	T[i] = 0;
      for(int i = 0; i < j; i++)
	T[i*j+i] = alpha[i+1];
      for(int i = 0; i < j-1; i++)
	{
	  T[i*j+(i+1)] = beta[i+1];
	  T[(i+1)*j+i] = beta[i+1];
	}

      ergo_real bestVector[maxIterations+1];
      currEig = get_eigs(j, T, bestVector, NULL);
      for(int k = 0; k < n; k++)
        {
          ergo_real sum = 0;
          for(int i = 1; i <= j; i++)
            sum += bestVector[i-1] * q[i][k];
          v[k] = sum;
        }

      curr_E = currEig + nucRepEnergy + shift;
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "current smallest eigenvalue: %22.11f  <->  E = %22.11f", (double)currEig, (double)curr_E);
    } // END FOR j

  return curr_E;
}




ergo_real do_lanczos_method(int n,
			    ergo_real* A,
			    ergo_real* v)
{
  ergo_real** q = new ergp_real_ptr[n];
  q[0] = new ergo_real[n];
  for(int i = 0; i < n; i++)
    q[0][i] = 0;
  q[1] = new ergo_real[n];
  for(int i = 0; i < n; i++)
    q[1][i] = v[i];
  normalize_vector(n, q[1]);
  ergo_real* z = new ergo_real[n];
  std::vector<ergo_real> alpha(n);
  std::vector<ergo_real> beta(n);
  beta[0] = 0;
  ergo_real currEig = 0;
  for(int j = 1; j < 22; j++)
    {
      for(int i = 0; i < n; i++)
	{
	  ergo_real sum = 0;
	  for(int k = 0; k < n; k++)
	    sum += A[i*n+k] * q[j][k];
	  z[i] = sum;
	} // END FOR i
      alpha[j] = 0;
      for(int i = 0; i < n; i++)
	alpha[j] += q[j][i] * z[i];
      for(int i = 0; i < n; i++)
	z[i] = z[i] - alpha[j] * q[j][i] - beta[j-1] * q[j-1][i];
      beta[j] = get_vector_norm(n, z);

      do_output(LOG_CAT_INFO, LOG_AREA_CI, "do_lanczos_method, j = %5i, alpha = %22.11f, beta = %22.11f", j, (double)alpha[j], (double)beta[j]);

      if(beta[j] < 1e-5)
	break;
      q[j+1] = new ergo_real[n];
      for(int i = 0; i < n; i++)
	q[j+1][i] = z[i] / beta[j];

      std::vector<ergo_real> T(j*j);
      for(int i = 0; i < j*j; i++)
	T[i] = 0;
      for(int i = 0; i < j; i++)
	T[i*j+i] = alpha[i+1];
      for(int i = 0; i < j-1; i++)
	{
	  T[i*j+(i+1)] = beta[i+1];
	  T[(i+1)*j+i] = beta[i+1];
	}
      currEig = get_eigs(j, &T[0], NULL, NULL);

    } // END FOR j

  return currEig;
}




int do_power_method(int n,
		    ergo_real* A,
		    ergo_real* v)
{
  ergo_real* v2 = new ergo_real[n];
  int iter = 0;
  while(1)
    {
      iter++;
      normalize_vector(n, v);
      for(int i = 0; i < n; i++)
	{
	  ergo_real sum = 0;
	  for(int j = 0; j < n; j++)
	    sum += A[i*n+j] * v[j];
	  v2[i] = sum;
	} // END FOR i
      ergo_real norm = get_vector_norm(n, v2);
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "power method, iter = %5i, norm = %22.11f", iter, (double)norm);
      for(int i = 0; i < n; i++)
	v[i] = v2[i];
    } // END WHILE
  return 0;
}


























#if 0

static ergo_real
do_dot_product(int n, const ergo_real* a, const ergo_real* b)
{
#if 1
  const int ONE = 1;
  return mat::dot(&n, a, &ONE, b, &ONE);
#else
  ergo_real sum = 0;
  for(int i = 0; i < n; i++)
    sum += a[i] * b[i];
  return sum;
#endif
}

static ergo_real
get_vector_length(int n, const ergo_real* a)
{
  return std::sqrt(do_dot_product(n, a, a));
}

static ergo_real
get_trace(int nSOs, two_idx_SO_struct A)
{
  ergo_real sum = 0;
  for(int i = 0; i < nSOs; i++)
    sum += A.x[i][i];
  return sum;
}

#endif





struct pqrs_struct
{
  int p;
  int q;
  int r;
  int s;
  void assign(int pp, int qq, int rr, int ss)
  {
    p = pp;
    q = qq;
    r = rr;
    s = ss;
  }
  int compare(int pp, int qq, int rr, int ss)
  {
    if(pp == p && qq == q && rr == r && ss == s)
      return 0;
    else
      return -1;
  }
};




#if 0
static int
do_experiment_1e2e(int nSlaterDets,
		   ergo_real* coeffs,
		   int nSOs,
		   int nEl,
		   const four_idx_SO_struct* g_SO, 
		   const two_idx_SO_struct* h_SO,
		   const SlaterDet_struct* SlaterDetList)
{
  
  printf("do_experiment_1e2e start, nSlaterDets = %5i.\n", nSlaterDets);
  
  // Randomize coeffs.
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] = rand_m1_to_1();
  normalize_vector(nSlaterDets, coeffs);

  // Compute 1el dens mat at base point.
  two_idx_SO_struct dmat1el;
  if(get_1e_density_matrix(nSOs,
			   nEl,
			   &dmat1el,
			   nSlaterDets, 
			   SlaterDetList, 
			   coeffs) != 0)
    {
      printf("error in get_1e_density_matrix\n");
      return -1;
    }
  printf("get_1e_density_matrix for base point OK.\n");

  // Compute 2el dens mat at base point.
  four_idx_SO_struct dmat2e;
  memset(&dmat2e, 0, sizeof(four_idx_SO_struct));
  if(get_2e_density_matrix(nSOs,
			   nEl,
			   &dmat2e,
			   nSlaterDets, 
			   SlaterDetList, 
			   coeffs) != 0)
    {
      printf("error in get_2e_density_matrix\n");
      return -1;
    }
  printf("get_2e_density_matrix for base point OK.\n");

  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < nSOs; k++)
	  sum += dmat2e.x[i][j][k][k];
	ergo_real guess = sum / (nEl-1);
	ergo_real diff = guess - dmat1el.x[i][j];
	ergo_real kvot = guess / dmat1el.x[i][j];
	printf("i j = %2i %2i, diff = %22.11f, kvot = %22.11f\n", i, j, (double)diff, (double)kvot);
      }

  printf("calling exit.\n");
  exit(0);

  return 0;
}


static int
do_experiment_1e(int nSlaterDets,
		 ergo_real* coeffs,
		 int nSOs,
		 int nEl,
		 const four_idx_SO_struct* g_SO, 
		 const two_idx_SO_struct* h_SO,
		 const SlaterDet_struct* SlaterDetList)
{
  
  printf("do_experiment_1e start, nSlaterDets = %5i.\n", nSlaterDets);
  
  // Randomize coeffs.
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] = rand_m1_to_1();
  normalize_vector(nSlaterDets, coeffs);

  // Compute 1el dens mat at base point.
  two_idx_SO_struct dmatbase;
  if(get_1e_density_matrix(nSOs,
			   nEl,
			   &dmatbase,
			   nSlaterDets, 
			   SlaterDetList, 
			   coeffs) != 0)
    {
      printf("error in get_1e_density_matrix\n");
      return -1;
    }
  printf("get_1e_density_matrix for base point OK.\n");

  two_idx_SO_struct dmatbase_2;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < nSOs; k++)
	  sum += dmatbase.x[i][k] * dmatbase.x[k][j];
	dmatbase_2.x[i][j] = sum;
      }
  
  two_idx_SO_struct dmatbase_4;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < nSOs; k++)
	  sum += dmatbase_2.x[i][k] * dmatbase_2.x[k][j];
	dmatbase_4.x[i][j] = sum;
      }
  
  printf("Tr(D ) = %22.11f\n", (double)get_trace(nSOs, dmatbase));
  printf("Tr(D2) = %22.11f\n", (double)get_trace(nSOs, dmatbase_2));
  printf("Tr(D4) = %22.11f\n", (double)get_trace(nSOs, dmatbase_4));

  

  
  const ergo_real h = 1e-6;

  const ergo_real lengthLimit = 0.001;

  const int maxiter = 444;

  ergo_real basisVectorList[maxiter][nSOs*nSOs];
  int nBasisVectors = 0;

  const int vectorN = nSOs*nSOs;
  printf("starting loop, nSOs = %5i, vectorN = %i\n", nSOs, vectorN);

  for(int iter = 0; iter < maxiter; iter++)
    {
      if(iter%10 == 0)
	printf("starting iter %5i\n", iter);

      // Generate random perturbation to coeffs
      ergo_real coeffsMod[nSlaterDets];
      for(int i = 0; i < nSlaterDets; i++)
	coeffsMod[i] = coeffs[i] + h * rand_m1_to_1();
      //normalize_vector(nSlaterDets, coeffsMod);
      
      // Create corresponding 1el dmat
      two_idx_SO_struct dmatmod;
      if(get_1e_density_matrix(nSOs,
			       nEl,
			       &dmatmod,
			       nSlaterDets, 
			       SlaterDetList, 
			       coeffsMod) != 0)
	{
	  printf("error in get_1e_density_matrix\n");
	  return -1;
	}
      //printf("get_1e_density_matrix for mod point OK, iter = %5i.\n", iter);

      ergo_real currVector[nSOs*nSOs];
      for(int i = 0; i < nSOs; i++)
	for(int j = 0; j < nSOs; j++)
	  currVector[i*nSOs+j] = dmatmod.x[i][j] - dmatbase.x[i][j];

      ergo_real originalLength  = get_vector_length(nSOs*nSOs, currVector);

      // Project away components of existing basis vectors.
      for(int i = 0; i < nBasisVectors; i++)
	{
	  ergo_real dotprod = do_dot_product(nSOs*nSOs, currVector, basisVectorList[i]);
	  //ergo_real len2 = get_vector_length(nSOs*nSOs, basisVectorList[i]);
	  ergo_real factor = dotprod;
	  //printf("i = %5i, len = %22.11f, len2 = %22.11f, factor = %22.11f\n", i, len, len2, factor);
	  for(int j = 0; j < nSOs*nSOs; j++)
	    currVector[j] -= factor * basisVectorList[i][j];
	}
      // Check length of remaining vector
      ergo_real remainingLength = get_vector_length(nSOs*nSOs, currVector);
      //printf("remainingLength = %22.11f\n", remainingLength);
      if(remainingLength/originalLength > lengthLimit)
	{
	  // Add new basis vector.
	  normalize_vector(nSOs*nSOs, currVector);
	  for(int i = 0; i < nSOs*nSOs; i++)
	    basisVectorList[nBasisVectors][i] = currVector[i];
	  nBasisVectors++;
	}
      else
	{
	  //printf("iter %i, too much projected away, skipping.\n", iter);
	}

    } // END FOR iter

  printf("Done, nBasisVectors = %5i\n", nBasisVectors);


  return 0;
}
#endif





#if 0
/* check4indeces: return 1 if this element should be stored, 0 otherwise. */
static int
check4indeces(int p, int q, int r, int s)
{
  if(p == r || q == s)
    {
      // element is always zero, should not be stored.
      return 0;
    }
  // Only one of pqrs rspq rqps psrq qpsr spqr qrsp srqp need to be stored.
  // p q r s   (1)
  // r s p q   (2)
  // r q p s   (3)
  // p s r q   (4)
  // q p s r   (5)
  // s p q r   (6)
  // q r s p   (7)
  // s r q p   (8)
  pqrs_struct list[8];
  list[0].assign(p, q, r, s);
  list[1].assign(r, s, p, q);
  list[2].assign(r, q, p, s);
  list[3].assign(p, s, r, q);
  list[4].assign(q, p, s, r);
  list[5].assign(s, p, q, r);
  list[6].assign(q, r, s, p);
  list[7].assign(s, r, q, p);
  // Sort list
  for(int i = 0; i < 7; i++)
    for(int j = i+1; j < 8; j++)
      {
	int doswitch = 0;
	if(list[i].p != list[j].p)
	  {
	    if(list[i].p > list[j].p)
	      doswitch = 1;
	  }
	else
	  {
	    if(list[i].q != list[j].q)
	      {
		if(list[i].q > list[j].q)
		  doswitch = 1;
	      }
	    else
	      {
		if(list[i].r != list[j].r)
		  {
		    if(list[i].r > list[j].r)
		      doswitch = 1;
		  }
		else
		  {
		    if(list[i].s != list[j].s)
		      {
			if(list[i].s > list[j].s)
			  doswitch = 1;
		      }		    
		  }
	      }
	  }
	if(doswitch == 1)
	  {
	    pqrs_struct temp = list[i];
	    list[i] = list[j];
	    list[j] = temp;
	  }
      }
  if(list[0].compare(p, q, r, s) == 0)
    return 1;
  else
    return 0;
}
#endif



#if 0
static int
do_experiment_2e(SlaterDet_struct** groupList,
		 int nSlaterDets,
		 ergo_real* coeffs,
		 int nSOs,
		 int nEl,
		 const four_idx_SO_struct* g_SO, 
		 const two_idx_SO_struct* h_SO,
		 const SlaterDet_struct* SlaterDetList)
{
  
  printf("do_experiment_2e start, nSlaterDets = %5i.\n", nSlaterDets);
  

#if 0
  // Randomize coeffs.
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] = rand_m1_to_1();
  normalize_vector(nSlaterDets, coeffs);
#endif


  // Compute 2el dens mat at base point.
  four_idx_SO_struct dmatbase2e;
  memset(&dmatbase2e, 0, sizeof(four_idx_SO_struct));
  two_idx_SO_struct dmatbase1e;
  memset(&dmatbase1e, 0, sizeof(two_idx_SO_struct));

#if 1
  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							 nEl,
							 0,
							 0,
							 0,
							 0,
							 0,
							 0,
							 0,
							 &status,
							 nSOs,
							 nSlaterDets,
							 SlaterDetList,
							 g_SO, 
							 h_SO, 
							 coeffs,
							 NULL,
							 &dmatbase1e,
							 &dmatbase2e
							 );
#else
  if(get_2e_density_matrix(nSOs,
			   nEl,
			   &dmatbase,
			   nSlaterDets, 
			   SlaterDetList, 
			   coeffs) != 0)
    {
      printf("error in get_2e_density_matrix\n");
      return -1;
    }
#endif

  printf("get_2e_density_matrix for base point OK.\n");

  const ergo_real h = 1e-6;

  const ergo_real lengthLimit = 0.005;

  const int maxiter = 14000;

  int vectorN = 0;

  int count = 0;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  {
	    if(i != j && i != k && i != l && j != k && j != l && k != l)
	      {
		// All indeces different. Store only if specific order.
		// First check symmetry.
		const ergo_real threshlimit = 1e-10;
		int errorflag = 0;
		if(std::fabs(dmatbase2e.x[i][j][k][l] - dmatbase2e.x[k][l][i][j]) > threshlimit)
		  errorflag++;
		if(std::fabs(dmatbase2e.x[i][j][k][l] + dmatbase2e.x[k][j][i][l]) > threshlimit)
		  errorflag++;
		if(std::fabs(dmatbase2e.x[i][j][k][l] + dmatbase2e.x[i][l][k][j]) > threshlimit)
		  errorflag++;
		if(std::fabs(dmatbase2e.x[i][j][k][l] - dmatbase2e.x[j][i][l][k]) > threshlimit)
		  errorflag++;
		if(std::fabs(dmatbase2e.x[i][j][k][l] + dmatbase2e.x[l][i][j][k]) > threshlimit)
		  errorflag++;
		if(std::fabs(dmatbase2e.x[i][j][k][l] + dmatbase2e.x[j][k][l][i]) > threshlimit)
		  errorflag++;
		if(std::fabs(dmatbase2e.x[i][j][k][l] - dmatbase2e.x[l][k][j][i]) > threshlimit)
		  errorflag++;
		if(i == k && std::fabs(dmatbase2e.x[i][j][k][l]) > threshlimit)
		  errorflag++;
		if(j == l && std::fabs(dmatbase2e.x[i][j][k][l]) > threshlimit)
		  errorflag++;
		if(errorflag != 0)
		  {
		    printf("ERROR!!!!!\n");
		    exit(0);
		  }
	      }
	    if(check4indeces(i, j, k, l) == 1)
	      {
		// Count this element.
		count++;
	      }
	  }
  vectorN = count;
  
  ergo_real basisVectorList[maxiter][vectorN];
  int nBasisVectors = 0;

  printf("starting loop, nSOs = %5i, vectorN = %i\n", nSOs, vectorN);

  for(int iter = 0; iter < maxiter; iter++)
    {
      if(iter%10 == 0)
	printf("starting iter %8i, nBasisVectors = %8i\n", iter, nBasisVectors);

      // Generate random perturbation to coeffs
      ergo_real coeffsMod[nSlaterDets];
      for(int i = 0; i < nSlaterDets; i++)
	coeffsMod[i] = coeffs[i] + h * rand_m1_to_1();
      //normalize_vector(nSlaterDets, coeffsMod);
      
      // Create corresponding 2el dmat
      four_idx_SO_struct dmatmod;
      memset(&dmatmod, 0, sizeof(four_idx_SO_struct));

#if 1
      pair_status_struct status;
      memset(&status, 0, sizeof(pair_status_struct));
      get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							     nEl,
							     0,
							     0,
							     0,
							     0,
							     0,
							     0,
							     0,
							     &status,
							     nSOs,
							     nSlaterDets,
							     SlaterDetList,
							     g_SO, 
							     h_SO, 
							     coeffsMod,
							     NULL,
							     &dmatbase1e,
							     &dmatmod
							     );
#else
      if(get_2e_density_matrix(nSOs,
			       nEl,
			       &dmatmod,
			       nSlaterDets, 
			       SlaterDetList, 
			       coeffsMod) != 0)
	{
	  printf("error in get_2e_density_matrix\n");
	  return -1;
	}
#endif

      ergo_real currVector[vectorN];
      int count = 0;
      for(int i = 0; i < nSOs; i++)
	for(int j = 0; j < nSOs; j++)
	  for(int k = 0; k < nSOs; k++)
	    for(int l = 0; l < nSOs; l++)
	      {
		if(i != j && i != k && i != l && j != k && j != l && k != l)
		  {
		    // All indeces different. Store only if specific order.
		    // First check symmetry.
		    const ergo_real threshlimit = 1e-10;
		    int errorflag = 0;
		    if(std::fabs(dmatmod.x[i][j][k][l] - dmatmod.x[k][l][i][j]) > threshlimit)
		      errorflag++;
		    if(std::fabs(dmatmod.x[i][j][k][l] + dmatmod.x[k][j][i][l]) > threshlimit)
		      errorflag++;
		    if(std::fabs(dmatmod.x[i][j][k][l] + dmatmod.x[i][l][k][j]) > threshlimit)
		      errorflag++;
		    if(errorflag != 0)
		      {
			printf("ERROR!!!!!\n");
			exit(0);
		      }
		  }	
		if(check4indeces(i, j, k, l) == 1)
		  {
		    // Store this element.
		    currVector[count] = dmatmod.x[i][j][k][l] - dmatbase2e.x[i][j][k][l];
		    count++;
		  }
	      }

      ergo_real originalLength  = get_vector_length(vectorN, currVector);
      //printf("originalLength = %22.11f\n", originalLength);


      // Project away components of existing basis vectors.
      for(int i = 0; i < nBasisVectors; i++)
	{
	  ergo_real dotprod = do_dot_product(vectorN, currVector, basisVectorList[i]);
	  //ergo_real len  = get_vector_length(vectorN, currVector);
	  //ergo_real len2 = get_vector_length(vectorN, basisVectorList[i]);
	  ergo_real factor = dotprod;
	  //printf("i = %5i, len = %22.11f, len2 = %22.11f, factor = %22.11f\n", i, len, len2, factor);
	  for(int j = 0; j < vectorN; j++)
	    currVector[j] -= factor * basisVectorList[i][j];
	}
      // Check length of remaining vector
      ergo_real remainingLength = get_vector_length(vectorN, currVector);
      //printf("remainingLength = %22.11f\n", remainingLength);
      if(remainingLength/originalLength > lengthLimit)
	{
	  // Add new basis vector.
	  normalize_vector(vectorN, currVector);
	  for(int i = 0; i < vectorN; i++)
	    basisVectorList[nBasisVectors][i] = currVector[i];
	  nBasisVectors++;
	  //printf("added new basis vector, nBasisVectors is now %5i\n", nBasisVectors);
	}
      else
	{
	  //printf("iter %i, too much projected away, skipping.\n", iter);
	}

    } // END FOR iter


  printf("Done, nBasisVectors = %5i\n", nBasisVectors);


  return 0;
}
#endif











#if 0

static int
get_2e_dmat_from_vector(four_idx_SO_struct* dmat2e, 
			int nSOs, 
			int vectorN, 
			const ergo_real* vector,
			const pqrs_struct* pqrsList)
{
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      for(int r = 0; r < nSOs; r++)
	for(int s = 0; s < nSOs; s++)
	  dmat2e->x[p][q][r][s] = 0;
  for(int i = 0; i < vectorN; i++)
    {
      int p = pqrsList[i].p;
      int q = pqrsList[i].q;
      int r = pqrsList[i].r;
      int s = pqrsList[i].s;
      ergo_real value = vector[i];
      dmat2e->x[p][q][r][s] = value;
      dmat2e->x[r][s][p][q] = value;
      dmat2e->x[r][q][p][s] = value * -1;
      dmat2e->x[p][s][r][q] = value * -1;
      dmat2e->x[q][p][s][r] = value;
      dmat2e->x[s][p][q][r] = value * -1;
      dmat2e->x[q][r][s][p] = value * -1;
      dmat2e->x[s][r][q][p] = value;
    } // END FOR i
  return 0;
}

#endif












#if 0

static int
do_experiment_2e_b(SlaterDet_struct** groupList,
		   int nSlaterDets,
		   ergo_real* coeffs,
		   int nSOs,
		   int nEl,
		   const four_idx_SO_struct* g_SO, 
		   const two_idx_SO_struct* h_SO,
		   const SlaterDet_struct* SlaterDetList)
{
  
  printf("do_experiment_2e_b start, nSlaterDets = %5i.\n", nSlaterDets);
  
  // Randomize coeffs.
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] = rand_m1_to_1();
  normalize_vector(nSlaterDets, coeffs);

  const ergo_real h = 1e-6;

  const ergo_real lengthLimit = 0.005;

  const int maxiter = 18;

  int vectorN = 0;

  int count = 0;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  {
	    if(check4indeces(i, j, k, l) == 1)
	      {
		// Count this element.
		count++;
	      }
	  }
  vectorN = count;
  printf("vectorN = %9i\n", vectorN);

  pqrs_struct pqrsList[vectorN];
  count = 0;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  {
	    if(check4indeces(i, j, k, l) == 1)
	      {
		// Count this element.
		pqrsList[count].p = i;
		pqrsList[count].q = j;
		pqrsList[count].r = k;
		pqrsList[count].s = l;
		count++;
	      }
	  }
  
  
  ergo_real list[maxiter][vectorN];
  int nMats = 0;

  printf("starting loop, nSOs = %5i, vectorN = %i\n", nSOs, vectorN);

  for(int iter = 0; iter < maxiter; iter++)
    {
      printf("iter = %5i\n", iter);

      // Randomize coeffs.
      for(int i = 0; i < nSlaterDets; i++)
	coeffs[i] = rand_m1_to_1();
      normalize_vector(nSlaterDets, coeffs);
      
      // Compute 2el dens mat at base point.
      four_idx_SO_struct dmatbase2e;
      memset(&dmatbase2e, 0, sizeof(four_idx_SO_struct));
      two_idx_SO_struct dmatbase1e;
      memset(&dmatbase1e, 0, sizeof(two_idx_SO_struct));

      pair_status_struct status;
      memset(&status, 0, sizeof(pair_status_struct));
      get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							     nEl,
							     0,
							     0,
							     0,
							     0,
							     0,
							     0,
							     0,
							     &status,
							     nSOs,
							     nSlaterDets,
							     SlaterDetList,
							     g_SO, 
							     h_SO, 
							     coeffs,
							     NULL,
							     &dmatbase1e,
							     &dmatbase2e
							     );

      printf("get_2e_density_matrix for base point OK.\n");

      int count = 0;
      for(int i = 0; i < nSOs; i++)
	for(int j = 0; j < nSOs; j++)
	  for(int k = 0; k < nSOs; k++)
	    for(int l = 0; l < nSOs; l++)
	      {
		if(check4indeces(i, j, k, l) == 1)
		  {
		    // Store this element.
		    list[iter][count] = dmatbase2e.x[i][j][k][l];
		    count++;
		  }
	      }
      if(count != vectorN)
	{
	  printf("ERROR: (count != vectorN)\n");
	  exit(0);
	}

      // test restoring 2el dmat
      four_idx_SO_struct dmat2;
      get_2e_dmat_from_vector(&dmat2, nSOs, vectorN, list[iter], pqrsList);
      for(int i = 0; i < nSOs; i++)
	for(int j = 0; j < nSOs; j++)
	  for(int k = 0; k < nSOs; k++)
	    for(int l = 0; l < nSOs; l++)
	      {
		ergo_real absdiff = std::fabs(dmatbase2e.x[i][j][k][l] - dmat2.x[i][j][k][l]);
		if(absdiff > 1e-11)
		  {
		    printf("ERRORRRRRRRRRRRR in restoring 2el dmat!!\n");
		    exit(0);
		  }
	      }

    } // END FOR iter
  printf("loop done, %i dmats created.\n", maxiter);

  const ergo_real tolerance = 1e-9;

  int keffCount = 0;
  for(int i = 0; i < vectorN; i++)
    for(int j = i+1; j < vectorN; j++)
      {
	ergo_real maxabsdiff = 0;
	for(int k = 0; k < maxiter; k++)
	  {
	    ergo_real absdiff = std::fabs(std::fabs(list[k][i]) - std::fabs(list[k][j]));
	    if(absdiff > maxabsdiff)
	      maxabsdiff = absdiff;
	  }
	if(maxabsdiff < tolerance)
	  {
	    printf("WARNING: elements %5i %5i are of nearly same magnitude in all cases.\n", i, j);
	    printf("%2i %2i %2i %2i\n", pqrsList[i].p, pqrsList[i].q, pqrsList[i].r, pqrsList[i].s);
	    printf("%2i %2i %2i %2i\n", pqrsList[j].p, pqrsList[j].q, pqrsList[j].r, pqrsList[j].s);
	    keffCount++;
	  }
      }
  printf("keffCount = %9i\n", keffCount);
  printf("DONE. Calling exit.\n");
  exit(0);

  return 0;
}

#endif






#if 0

static int
get_1e_dmat_from_2e_dmat(two_idx_SO_struct* dmat1e, 
			 const four_idx_SO_struct* dmat2e, 
			 int nSOs,
			 int nEl)
{
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < nSOs; k++)
	  sum += dmat2e->x[i][j][k][k];
	ergo_real dmat1e_element = sum / (nEl-1);
	dmat1e->x[i][j] = dmat1e_element;
      }
  return 0;
}

#endif





#if 0

static ergo_real 
get_trace_of_1e_dmat(int nSOs, 
		     const two_idx_SO_struct* dmat1e)
{
  ergo_real sum = 0;
  for(int i = 0; i < nSOs; i++)
    sum += dmat1e->x[i][i];
  return sum;
}

static void 
scale_1e_dmat(int nSOs, 
	      two_idx_SO_struct* dmat1e,
	      ergo_real scaleFactor)
{
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      dmat1e->x[i][j] *= scaleFactor;
}

static void 
scale_2e_dmat(int nSOs, 
	      four_idx_SO_struct* dmat2e,
	      ergo_real scaleFactor)
{
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  dmat2e->x[i][j][k][l] *= scaleFactor;
}

#endif



#if 0

static int
get_2idx_matrix_from_4idx_matrix(int nSOs, int n, ergo_real* A, const four_idx_SO_struct* dmat2e)
{
  memset(A, 0, n*n*sizeof(ergo_real));
  int* Atest = new int[n*n];
  memset(Atest, 0, n*n*sizeof(int));
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      for(int r = 0; r < nSOs; r++)
	for(int s = 0; s < nSOs; s++)
	  {
	    if(p > q && r > s)
	      {
		int i = q*nSOs + p - (q+1)*(q+2)/2;
		int j = s*nSOs + r - (s+1)*(s+2)/2;
		A[i*n+j] = dmat2e->x[p][r][q][s];
		Atest[i*n+j]++;
	      }
	  }
  for(int i = 0; i < n*n; i++)
    {
      if(Atest[i] != 1)
	{
	  printf("ERROR: (Atest[i] != 1)\n");
	  exit(0);
	}
    }
  // Check that A is symmetric
  ergo_real maxabsdiff = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real absdiff = std::fabs(A[i*n+j] - A[j*n+i]);
	if(absdiff > maxabsdiff)
	  maxabsdiff = absdiff;
      }
  if(maxabsdiff > 1e-11)
    {
      printf("ERROR: A not symmetric!\n");
      exit(0);
    }
  //printf("symmerty check of A: maxabsdiff = %22.11g\n");

  delete [] Atest;

  return 0;
}

#endif



#if 0

static void
get_Q_and_G_from_dmats(int nSOs, 
		       four_idx_SO_struct* Q, 
		       four_idx_SO_struct* G, 
		       const two_idx_SO_struct* dmat1e,
		       const four_idx_SO_struct* dmat2e)
{
  two_idx_SO_struct id1e;
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      {
	if(p == q)
	  id1e.x[p][q] = 1;
	else
	  id1e.x[p][q] = 0;
      }
  four_idx_SO_struct id2e;
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      for(int r = 0; r < nSOs; r++)
	for(int s = 0; s < nSOs; s++)
	  {
	    ergo_real id2epqrs = 0;
	    if(p == q && r == s)
	      id2epqrs += 0.5;
	    if(q == r && p == s)
	      id2epqrs += 0.5;
	    id2e.x[p][q][r][s] = id2epqrs;
	  }
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      for(int r = 0; r < nSOs; r++)
	for(int s = 0; s < nSOs; s++)
	  {
	    // Get element for Q
	    ergo_real Qpqrs = 0;
	    // First term: 2*identity
	    Qpqrs += 2 * id2e.x[p][q][r][s];
	    // Add wedge product, 4 terms
	    Qpqrs -= dmat1e->x[p][q] * id1e.x[r][s];
	    Qpqrs += dmat1e->x[r][q] * id1e.x[p][s];
	    Qpqrs += dmat1e->x[p][s] * id1e.x[r][q];
	    Qpqrs -= dmat1e->x[r][s] * id1e.x[p][q];
	    // Final term: dmat2e
	    Qpqrs += dmat2e->x[p][q][r][s];
	    Q->x[p][q][r][s] = Qpqrs;

	    // Get element for G
	    ergo_real Gpqrs = 0;
	    // First term: identity * dmat1e
	    Gpqrs += id1e.x[r][s] * dmat1e->x[p][q];
	    // Second term: dmat2e with switched indeces
	    Gpqrs -= dmat2e->x[p][q][s][r];
	    G->x[p][q][r][s] = Gpqrs;
	  }
}


static ergo_real
get_penalty(int nSOs, const four_idx_SO_struct* dmat2e)
{
  // Get 2-idx matrix A from dmat2e
  int n = nSOs*(nSOs-1) / 2;
  ergo_real* A = new ergo_real[n*n];
  get_2idx_matrix_from_4idx_matrix(nSOs, n, A, dmat2e);
  const ergo_real penaltyParam = 1e2;
  // Get all eigenvalues of A
  ergo_real* eigValList = new ergo_real[n];
  get_eigs(n, A, NULL, eigValList);
  ergo_real eigValMin =  1e22;
  ergo_real eigValMax = -1e22;
  int negCount = 0;
  ergo_real penalty = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real curr = eigValList[i];
      if(curr > eigValMax)
	eigValMax = curr;
      if(curr < eigValMin)
	eigValMin = curr;
      if(curr < 0)
	{
	  negCount++;
	  penalty += penaltyParam * std::fabs(curr) * std::fabs(curr);
	}
    }
  //printf("eigValMin = %22.14f\n", eigValMin);
  //printf("eigValMax = %22.14f\n", eigValMax);
  //printf("n = %8i, negCount = %8i\n", n, negCount);
  delete [] eigValList;
  delete [] A;
  return penalty;
}

#endif





#if 0

static ergo_real get_energy_from_2e_vector(int nSOs,
					   int nEl,
					   const four_idx_SO_struct* g_SO, 
					   const two_idx_SO_struct* h_SO,
					   int vectorN,
					   const ergo_real* vector,
					   const pqrs_struct* pqrsList,
					   ergo_real nuclearEnergy,
					   int do_eigs,
					   int do_penalty)
{
  // Create 2e dmat from vector.
  four_idx_SO_struct dmat2e;
  get_2e_dmat_from_vector(&dmat2e, nSOs, vectorN, vector, pqrsList);

  // Create 1e dmat from 2e dmat
  two_idx_SO_struct dmat1e;
  get_1e_dmat_from_2e_dmat(&dmat1e, &dmat2e, nSOs, nEl);

  // Compute trace of 1e dmat
  ergo_real trace = get_trace_of_1e_dmat(nSOs, &dmat1e);
  ergo_real scaleFactor = (ergo_real)nEl / trace;
   
  // Normalize 1el and 2el dmats
  scale_1e_dmat(nSOs, &dmat1e, scaleFactor);
  scale_2e_dmat(nSOs, &dmat2e, scaleFactor);

  // Compute 1e energy
  ergo_real energy_1el = 0;
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      energy_1el += dmat1e.x[p][q] * h_SO->x[p][q];

  // Compute 2e energy
  ergo_real energy_2el = 0;
  for(int p = 0; p < nSOs; p++)
    for(int q = 0; q < nSOs; q++)
      for(int r = 0; r < nSOs; r++)
	for(int s = 0; s < nSOs; s++)
	  energy_2el += 0.5 * dmat2e.x[p][q][r][s] * g_SO->x[p][q][r][s];

  ergo_real E = energy_1el + energy_2el + nuclearEnergy;


#if 0
  // Get 2-idx matrix A from dmat2e
  int n = nSOs*(nSOs-1) / 2;
  ergo_real* A = new ergo_real[n*n];
  get_2idx_matrix_from_4idx_matrix(nSOs, n, A, &dmat2e);
  const ergo_real penaltyParam = 1e3;
  if(do_eigs)
    {
      // Get all eigenvalues of A
      ergo_real* eigValList = new ergo_real[n];
      get_eigs(n, A, NULL, eigValList);
      ergo_real eigValMin =  1e22;
      ergo_real eigValMax = -1e22;
      int negCount = 0;
      for(int i = 0; i < n; i++)
	{
	  ergo_real curr = eigValList[i];
	  if(curr > eigValMax)
	    eigValMax = curr;
	  if(curr < eigValMin)
	    eigValMin = curr;
	  if(curr < 0)
	    {
	      negCount++;
	      if(do_penalty)
		E += penaltyParam * std::fabs(curr) * std::fabs(curr);
	    }
	}
      //printf("eigValMin = %22.14f\n", eigValMin);
      //printf("eigValMax = %22.14f\n", eigValMax);
      //printf("n = %8i, negCount = %8i\n", n, negCount);
      delete [] eigValList;
    }
  delete [] A;
#endif


  // Compute Q and G
  four_idx_SO_struct Q;
  four_idx_SO_struct G;
  get_Q_and_G_from_dmats(nSOs, 
			 &Q, 
			 &G, 
			 &dmat1e,
			 &dmat2e);

  E += get_penalty(nSOs, &dmat2e);
  E += get_penalty(nSOs, &Q);
  E += get_penalty(nSOs, &G);
   
  return E;
}

#endif







#if 0
static int
do_experiment_2e_optimize(SlaterDet_struct** groupList,
			  int nSlaterDets,
			  ergo_real* coeffs,
			  int nSOs,
			  int nEl,
			  const four_idx_SO_struct* g_SO, 
			  const two_idx_SO_struct* h_SO,
			  const SlaterDet_struct* SlaterDetList,
			  ergo_real nuclearEnergy)
{
  
  printf("do_experiment_2e_optimize start, nSlaterDets = %5i.\n", nSlaterDets);
  
#if 0
  // HF guess coeffs.
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] = 0;
  coeffs[0] = 1;
  normalize_vector(nSlaterDets, coeffs);
#endif

#if 0
  // Randomize coeffs.
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] = rand_m1_to_1();
  normalize_vector(nSlaterDets, coeffs);
#endif

#if 0
  // Disturb coeffs
  for(int i = 0; i < nSlaterDets; i++)
    coeffs[i] += 0.05 * rand_m1_to_1();
  normalize_vector(nSlaterDets, coeffs);
#endif


  // Compute 2el dens mat at base point.
  four_idx_SO_struct dmatbase2e;
  memset(&dmatbase2e, 0, sizeof(four_idx_SO_struct));
  two_idx_SO_struct dmatbase1e;
  memset(&dmatbase1e, 0, sizeof(two_idx_SO_struct));

  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							 nEl,
							 0,
							 0,
							 0,
							 0,
							 0,
							 0,
							 0,
							 &status,
							 nSOs,
							 nSlaterDets,
							 SlaterDetList,
							 g_SO, 
							 h_SO, 
							 coeffs,
							 NULL,
							 &dmatbase1e,
							 &dmatbase2e
							 );

  printf("get 1e and 2e density_matrix for base point OK.\n");






  int count = 0;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  {
	    if(check4indeces(i, j, k, l) == 1)
	      {
		// Count this element.
		count++;
	      }
	  }
  int vectorN = count;
  printf("vectorN = %9i\n", vectorN);

  pqrs_struct pqrsList[vectorN];
  count = 0;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  {
	    if(check4indeces(i, j, k, l) == 1)
	      {
		// Count this element.
		pqrsList[count].p = i;
		pqrsList[count].q = j;
		pqrsList[count].r = k;
		pqrsList[count].s = l;
		count++;
	      }
	  }


  ergo_real currVector[vectorN];


  count = 0;
  for(int i = 0; i < nSOs; i++)
    for(int j = 0; j < nSOs; j++)
      for(int k = 0; k < nSOs; k++)
	for(int l = 0; l < nSOs; l++)
	  {
	    if(check4indeces(i, j, k, l) == 1)
	      {
		// Store this element.
		currVector[count] = dmatbase2e.x[i][j][k][l];
		count++;
	      }
	  }
  if(count != vectorN)
    {
      printf("ERROR: (count != vectorN)\n");
      exit(0);
    }


  // Now start optimization.
  
  printf("init done, now starting optimization...\n");

  // currVector is the current guess.


  const ergo_real h       = 0.000005;
  ergo_real stepLen = 0.0005;
  ergo_real Eprev = 1e22;
  int step = 0;
  
  while(1)
    {
      step++;

      normalize_vector(vectorN, currVector);

      // Get current energy.
      ergo_real E = get_energy_from_2e_vector(nSOs,
					      nEl,
					      g_SO,
					      h_SO,
					      vectorN,
					      currVector,
					      pqrsList,
					      nuclearEnergy,
					      1,
					      1);
      //printf("E = %22.11f\n", E);
      
      if(E > Eprev)
	{
	  stepLen *= 0.5;
	  printf("Energy increased; shortening stepLen to %22.11g\n", (double)stepLen);
	}
      else
	stepLen *= 1.05;

      if(stepLen < 1e-22)
	{
	  printf("stepLen very small, stopping now.\n");
	  break;
	}

      Eprev = E;

      // Get gradient.
      ergo_real Emin = E;
      int iEmin = 0;
      ergo_real gradient[vectorN];
      for(int i = 0; i < vectorN; i++)
	{
	  if((i+1) % 10000 == 0)
	    printf("getting gradient, i = %9i\n", i);
	  currVector[i] += h;	  
	  ergo_real E2 = get_energy_from_2e_vector(nSOs,
						   nEl,
						   g_SO,
						   h_SO,
						   vectorN,
						   currVector,
						   pqrsList,
						   nuclearEnergy,
						   1,
						   1);
	  if(E2 < Emin)
	    {
	      Emin = E2;
	      iEmin = i;
	    }
	  gradient[i] = (E2 - E) / h;
	  currVector[i] -= h;
	  
	}
      ergo_real normsum = 0;
      for(int i = 0; i < vectorN; i++)
	normsum += gradient[i]*gradient[i];
      ergo_real gradnorm = std::sqrt(normsum);
      printf("step %7i, E = %22.11f, gradnorm = %22.11f.\n", step, (double)E, (double)gradnorm);

      int bestIndex = 0;
      for(int i = 0; i < vectorN; i++)
	{
	  if(std::fabs(gradient[i]) > std::fabs(gradient[bestIndex]))
	    bestIndex = i;
	}

      normalize_vector(vectorN, gradient);
            
      for(int i = 0; i < vectorN; i++)
	currVector[i] -= stepLen * gradient[i];

    } // END WHILE


  return 0;
}
#endif


































int do_CI(
	  const BasisInfoStruct & basisInfo, 
	  const IntegralInfo & integralInfo,
	  const CI::Options& options,
	  const ergo_real* S,
	  const ergo_real* h_AO,
	  const ergo_real* F_a,
	  const ergo_real* F_b,
	  int n_el_a,
	  int n_el_b,
	  ergo_real nuclearEnergy,
	  ergo_real HF_energy
	  )
{
  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "do_CI, n = %3i, HF_energy = %22.11f", n, (double)HF_energy);

  if(n > MAX_AOS)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "Error in do_CI: (n > MAX_AOS).");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_CI, "beginning of do_CI");

  // Use basisInfo to get g_AO
  four_idx_AO_struct* g_AO = new four_idx_AO_struct;
  int p, q, r, s;
  for(p = 0; p < n; p++)
    for(q = 0; q < n; q++)
      for(r = 0; r < n; r++)
	for(s = 0; s < n; s++)
	  {
	    g_AO->x[p][q][r][s] = do_2e_integral(p, q, r, s, 
						 basisInfo, 
						 integralInfo);
	  }
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "g_AO done");
  
  ergo_real* MOs_a = new ergo_real[n*n];
  ergo_real* MOs_b = new ergo_real[n*n];
  ergo_real* eigv_a = new ergo_real[n];
  ergo_real* eigv_b = new ergo_real[n];


  if(options.use_random_orbitals == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using random orbitals, different for alpha and beta.");
      ergo_real* F_random = new ergo_real[n*n];
      for(p = 0; p < n*n; p++)
	F_random[p] = rand_m1_to_1();
      get_F_orbs(n, F_random, S, MOs_a, eigv_a);
      for(p = 0; p < n*n; p++)
	F_random[p] = rand_m1_to_1();
      get_F_orbs(n, F_random, S, MOs_b, eigv_b);
    }
  else if(options.use_lowdin_orbitals == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using Lowdin orbitals, same for alpha and beta.");
      get_Lowdin_orbitals(n, S, MOs_a);
      get_Lowdin_orbitals(n, S, MOs_b);
    }
  else
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using HF orbitals obtained from alpha and beta Fock matrices.");
      get_F_orbs(n, F_a, S, MOs_a, eigv_a);  
      get_F_orbs(n, F_b, S, MOs_b, eigv_b);
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "MOs_a MOs_b done.");

  int nSOs = 2 * n;
  
  SO_struct* SOs = new SO_struct[nSOs];
  int count = 0;
  int i;
  // ALPHA
  for(i = 0; i < n; i++)
    {
      SOs[count].spin = SPIN_A;
      int j;
      for(j = 0; j < n; j++)
	SOs[count].coeffs[j] = MOs_a[i*n+j];
      count++;
    }
  // BETA
  for(i = 0; i < n; i++)
    {
      SOs[count].spin = SPIN_B;
      int j;
      for(j = 0; j < n; j++)
	SOs[count].coeffs[j] = MOs_b[i*n+j];
      count++;
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "SOs done.");

  four_idx_SO_struct* g_SO = new four_idx_SO_struct;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      for(r = 0; r < nSOs; r++)
	for(s = 0; s < nSOs; s++)
	  {
	    ergo_real value = 0;
	    if(SOs[p].spin == SOs[q].spin && SOs[r].spin == SOs[s].spin)
	      {
		int a, b, c, d;
		for(a = 0; a < n; a++)
		  for(b = 0; b < n; b++)
		    for(c = 0; c < n; c++)
		      for(d = 0; d < n; d++)
			value += SOs[p].coeffs[a] * SOs[q].coeffs[b] * SOs[r].coeffs[c] * SOs[s].coeffs[d] * g_AO->x[a][b][c][d];
	      }
	    g_SO->x[p][q][r][s] = value;
	  }
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "g_SO done");

  two_idx_SO_struct* h_SO = new two_idx_SO_struct;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      {
	ergo_real value = 0;
	if(SOs[p].spin == SOs[q].spin)
	  {
	    int a, b;
	    for(a = 0; a < n; a++)
	      for(b = 0; b < n; b++)
		value += SOs[p].coeffs[a] * SOs[q].coeffs[b] * h_AO[a*n+b];
	  }
	h_SO->x[p][q] = value;
      }
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "h_SO done.");

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "n_el_a = %2i, n_el_b = %2i", n_el_a, n_el_b);

  int nElTot = n_el_a + n_el_b;

  if(nElTot > MAX_ELECTRONS)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "Error in do_CI: (nElTot > MAX_ELECTRONS).");
      return -1;
    }

  int nSlaterDets1 = get_FCI_Slater_dets_alpha_beta(NULL, n_el_a, n_el_b, nSOs);
  SlaterDet_struct* SlaterDetList1 = new SlaterDet_struct[nSlaterDets1];
  get_FCI_Slater_dets_alpha_beta(SlaterDetList1, n_el_a, n_el_b, nSOs);

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "Slater determinants done, nSlaterDets = %5i", nSlaterDets1);

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "Computing energy for each Slater determinant, nSlaterDets = %9i.", nSlaterDets1);
  std::vector<ergo_real> SlaterDetEnergyList(nSlaterDets1);
  for(int i = 0; i < nSlaterDets1; i++)
    {
      SlaterDetEnergyList[i] = get_SlaterDet_energy(nSOs,
						    nElTot,
						    g_SO, 
						    h_SO, 
						    &SlaterDetList1[i]);
    }
#if 0
  for(int i = 0; i < nSlaterDets1; i++)
    printf("SlaterDetEnergyList[%9i] = %22.11f  <->  E = %22.11f\n", i, SlaterDetEnergyList[i], SlaterDetEnergyList[i] + nuclearEnergy);
#endif

  SlaterDet_struct* SlaterDetList = new SlaterDet_struct[nSlaterDets1];
  int nSlaterDets;
  if(options.use_energy_diff_limit == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Selecting Slater dets using energy_diff_limit = %9.5f.", (double)options.energy_diff_limit);
      int detCount = 0;
      for(int i = 0; i < nSlaterDets1; i++)
	{
	  if(SlaterDetEnergyList[i] + nuclearEnergy - HF_energy < options.energy_diff_limit)
	    {
	      SlaterDetList[detCount] = SlaterDetList1[i];
	      detCount++;
	    }
	}
      nSlaterDets = detCount;
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Selected Slater dets low enough in energy, nSlaterDets = %9i.", nSlaterDets);
    }
  else
    {
      for(int i = 0; i < nSlaterDets1; i++)
	SlaterDetList[i] = SlaterDetList1[i];
      nSlaterDets = nSlaterDets1;
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using all Slater dets, nSlaterDets = %9i.", nSlaterDets);
    }

  if(nSlaterDets < 1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "Error: nSlaterDets = %i", nSlaterDets);
      return -1;
    }
  

  // Define groups of Slater determinants, one level at a time
  SlaterDet_struct* groupList[MAX_ELECTRONS+1];
  int groupCountList[MAX_ELECTRONS+1];
  groupList[nElTot] = SlaterDetList;
  groupCountList[nElTot] = nSlaterDets;
  SlaterDet_struct* list = new SlaterDet_struct[nSlaterDets];
  for(int i = 0; i < nElTot; i++)
    {
      // Define group with nElTot-1-i electrons.
      int nElCurrLevel = nElTot-1-i;
      int prevGroupIndex = nElTot-i;
      // Go through prev level to find groups at curr level.
      SlaterDet_struct* prevGroupList = groupList[prevGroupIndex];
      int prevCount = groupCountList[prevGroupIndex];
      int j = 0;
      int nGroupsCurrLevel = 0;
      while(j < prevCount)
	{
	  // Now j points to beginning of a new group
	  SlaterDet_struct newGroup;
	  memset(&newGroup, 0, sizeof(SlaterDet_struct));
	  for(int k = 0; k < nElCurrLevel; k++)
	    newGroup.SO_list[k] = prevGroupList[j].SO_list[k];
	  newGroup.startIndex = j;
	  int count = 0;
	  while(j < prevCount)
	    {
	      int ok = 1;
	      for(int k = 0; k < nElCurrLevel; k++)
		{
		  if(prevGroupList[j].SO_list[k] != newGroup.SO_list[k])
		    ok = 0;
		}
	      if(ok)
		{
		  j++;
		  count++;
		}
	      else
		break;
	    }
	  // Now j points to the beginning of next group or to the end of the list.
	  newGroup.count = count;
	  list[nGroupsCurrLevel] = newGroup;
	  nGroupsCurrLevel++;
	} // END WHILE
      int groupIndex = nElTot-1-i;
      groupList[groupIndex] = new SlaterDet_struct[nGroupsCurrLevel];
      memcpy(groupList[groupIndex], list, nGroupsCurrLevel*sizeof(SlaterDet_struct));
      groupCountList[groupIndex] = nGroupsCurrLevel;
    }
  delete list;
  list = NULL;
  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "Slater determinant groups done:");
  for(int i = 0; i < nElTot+1; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_CI, "groupCountList[%2i] = %9i", i, groupCountList[i]);

  

#if 0
  // Get number of relevant pairs of SlaterDets
  int nRelevantSlaterDetPairs = get_relevant_SlaterDet_pairs(nSlaterDets, SlaterDetList, groupList, nElTot, NULL);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "nRelevantSlaterDetPairs = %15i", nRelevantSlaterDetPairs);
  SlaterDet_pair_struct* SlaterDet_pair_list = new SlaterDet_pair_struct[nRelevantSlaterDetPairs];
  get_relevant_SlaterDet_pairs(nSlaterDets, SlaterDetList, groupList, nElTot, SlaterDet_pair_list);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "get_relevant_SlaterDet_pairs done (twice)");
  
  int* pairCountList = new int[nSlaterDets];
  for(int i = 0; i < nSlaterDets; i++)
    pairCountList[i] = 0;
  for(int i = 0; i < nRelevantSlaterDetPairs; i++)
    pairCountList[SlaterDet_pair_list[i].a]++;
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "pairCountList done.");
#endif


  ergo_real* coeffList = new ergo_real[nSlaterDets];
  if(options.use_random_starting_guess == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using random CI starting guess.");
      for(i = 0; i < nSlaterDets; i++)
	coeffList[i] = rand_m1_to_1();
    }
  else
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using vector [ 1 0 0 ... 0 ] as CI starting guess.");
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "If HF orbitals are used, this should correspond to the HF state.");
      for(i = 0; i < nSlaterDets; i++)
	coeffList[i] = 0;
      coeffList[0] = 1;
    }
  normalize_vector(nSlaterDets, coeffList);

  ergo_real CI_energy = do_lanczos_method_direct(nSlaterDets,
						 coeffList,
						 nSOs,
						 nElTot,
						 g_SO, 
						 h_SO, 
						 nSlaterDets, 
						 SlaterDetList, 
						 groupList,
						 options.max_no_of_iterations,
						 nuclearEnergy,
						 options.shift
						 );


  do_output(LOG_CAT_INFO, LOG_AREA_CI, "FINAL CI ENERGY = %22.11f", (double)CI_energy);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "FINAL CI CORRELATION ENERGY = %22.11f", (double)(CI_energy - HF_energy));

  return 0;
}


