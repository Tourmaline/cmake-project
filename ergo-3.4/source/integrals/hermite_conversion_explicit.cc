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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <cmath>
#include "hermite_conversion_symb.h"
#include "output.h"



int 
get_hermite_conversion_matrix(const monomial_info_struct* monomial_info,
			      int nmax,
			      int inverseFlag,
			      ergo_real exponent,
			      ergo_real* result)
{
  int noOfMonomials = monomial_info->no_of_monomials_list[nmax];
  symb_matrix_element result_symb[noOfMonomials*noOfMonomials];
  get_hermite_conversion_matrix_symb(monomial_info,
				     nmax,
				     inverseFlag,
				     result_symb);
  for(int i = 0; i < noOfMonomials*noOfMonomials; i++)
    {
      ergo_real factor = std::pow(exponent, result_symb[i].ia);
      result[i] = result_symb[i].coeff * factor;
    }
  return 0;
}





#if 0

typedef struct
{
  int ix;
  ergo_real coeff;
} poly_1d_term_struct;

#define MAX_NO_OF_1D_TERMS 888

typedef struct
{
  int noOfTerms;
  poly_1d_term_struct termList[MAX_NO_OF_1D_TERMS];
} poly_1d_struct;


typedef struct
{
  int monomialInts[3];
  ergo_real coeff;
} poly_3d_term_struct;

#define MAX_NO_OF_3D_TERMS 888

typedef struct
{
  int noOfTerms;
  poly_3d_term_struct termList[MAX_NO_OF_3D_TERMS];
} poly_3d_struct;


static int
get_1d_hermite_poly(poly_1d_struct* result, int n, ergo_real a)
{
  switch(n)
    {
    case 0:
      result->noOfTerms = 1;
      result->termList[0].ix = 0;
      result->termList[0].coeff = 1;
      break;
    case 1:
      result->noOfTerms = 1;
      result->termList[0].ix = 1;
      result->termList[0].coeff = 2 * a;
      break;
    default:
      {
	// Create polys for n-1 and n-2
	poly_1d_struct poly_n_m_1;
	poly_1d_struct poly_n_m_2;
	get_1d_hermite_poly(&poly_n_m_1, n - 1, a);
	get_1d_hermite_poly(&poly_n_m_2, n - 2, a);
	if(poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms >= MAX_NO_OF_1D_TERMS)
	  {
	    printf("error: (poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms >= MAX_NO_OF_1D_TERMS)\n");
	    exit(0);
	  }
	// Now the result is 2*a*x*poly_n_m_1 - (n-1)*2*a*poly_n_m_2
	int i;
	for(i = 0; i < poly_n_m_1.noOfTerms; i++)
	  {
	    result->termList[i] = poly_n_m_1.termList[i];
	    result->termList[i].ix++;
	    result->termList[i].coeff *= 2 * a;
	  }
	int nn = poly_n_m_1.noOfTerms;
	for(i = 0; i < poly_n_m_2.noOfTerms; i++)
	  {
	    result->termList[nn+i] = poly_n_m_2.termList[i];
	    result->termList[nn+i].coeff *= -2 * (n-1) * a;
	  }
	result->noOfTerms = poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms;
      }
    }
  return 0;
}

static int
get_1d_hermite_poly_inv(poly_1d_struct* result, int n, ergo_real a)
{
  switch(n)
    {
    case 0:
      result->noOfTerms = 1;
      result->termList[0].ix = 0;
      result->termList[0].coeff = 1;
      break;
    case 1:
      result->noOfTerms = 1;
      result->termList[0].ix = 1;
      result->termList[0].coeff = 0.5 * a;
      break;
    default:
      {
	// Create polys for n-1 and n-2
	poly_1d_struct poly_n_m_1;
	poly_1d_struct poly_n_m_2;
	get_1d_hermite_poly_inv(&poly_n_m_1, n - 1, a);
	get_1d_hermite_poly_inv(&poly_n_m_2, n - 2, a);
	if(poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms >= MAX_NO_OF_1D_TERMS)
	  {
	    printf("error: (poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms >= MAX_NO_OF_1D_TERMS)\n");
	    exit(0);
	  }
	// Now the result is 0.5*a*x*poly_n_m_1 + (n-1)*a*0.5*poly_n_m_2
	int i;
	for(i = 0; i < poly_n_m_1.noOfTerms; i++)
	  {
	    result->termList[i] = poly_n_m_1.termList[i];
	    result->termList[i].ix++;
	    result->termList[i].coeff *= 0.5 * a;
	  }
	int nn = poly_n_m_1.noOfTerms;
	for(i = 0; i < poly_n_m_2.noOfTerms; i++)
	  {
	    result->termList[nn+i] = poly_n_m_2.termList[i];
	    result->termList[nn+i].coeff *= (n-1) * 0.5 * a;
	  }
	result->noOfTerms = poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms;
      }
    }
  return 0;
}


static int
create_3d_poly_from_1d_poly(poly_3d_struct* poly_3d, poly_1d_struct* poly_1d, int coordIndex)
{
  memset(poly_3d, 0, sizeof(poly_3d_struct));
  int i;
  for(i = 0; i < poly_1d->noOfTerms; i++)
    {
      poly_3d->termList[i].coeff = poly_1d->termList[i].coeff;
      poly_3d->termList[i].monomialInts[coordIndex] = poly_1d->termList[i].ix;
    }
  poly_3d->noOfTerms = poly_1d->noOfTerms;
  return 0;
}


static int
compute_product_of_3d_polys(poly_3d_struct* result, 
			    poly_3d_struct* poly_1, 
			    poly_3d_struct* poly_2)
{
  int termCount = 0;
  int termidx_1, termidx_2;
  for(termidx_1 = 0; termidx_1 < poly_1->noOfTerms; termidx_1++)
    for(termidx_2 = 0; termidx_2 < poly_2->noOfTerms; termidx_2++)
      {
	poly_3d_term_struct* term_1 = &poly_1->termList[termidx_1];
	poly_3d_term_struct* term_2 = &poly_2->termList[termidx_2];
	
	// Create product term
	poly_3d_term_struct newTerm;
	newTerm.coeff = term_1->coeff * term_2->coeff;
	int k;
	for(k = 0; k < 3; k++)
	  newTerm.monomialInts[k] = term_1->monomialInts[k] + term_2->monomialInts[k];

	result->termList[termCount] = newTerm;
	termCount++;
	if(termCount >= MAX_NO_OF_3D_TERMS)
	  {
	    printf("ERROR: (termCount >= MAX_NO_OF_3D_TERMS)\n");
	    exit(0);
	  }
      } // END FOR termidx_1 termidx_2
  result->noOfTerms = termCount;
  return 0;
}


static int 
get_hermite_conversion_matrix(const IntegralInfo* integralInfo,
			      int nmax, 
			      int inverseFlag,
			      ergo_real exponent,
			      ergo_real* result)
{
  int noOfMonomials = integralInfo->monomial_info.no_of_monomials_list[nmax];
  memset(result, 0, noOfMonomials*noOfMonomials*sizeof(ergo_real));

#if 1
  symb_matrix_element result_symb[noOfMonomials*noOfMonomials];
  get_hermite_conversion_matrix_symb(&integralInfo->monomial_info,
				     nmax,
				     inverseFlag,
				     result_symb);
  for(int i = 0; i < noOfMonomials*noOfMonomials; i++)
    {
      ergo_real factor = std::pow(exponent, result_symb[i].ia);
      result[i] = result_symb[i].coeff * factor;
    }
  return 0;
#endif

  int monomialIndex;
  for(monomialIndex = 0; monomialIndex < noOfMonomials; monomialIndex++)
    {

      // get monomialInts
      int ix = integralInfo->monomial_info.monomial_list[monomialIndex].ix;
      int iy = integralInfo->monomial_info.monomial_list[monomialIndex].iy;
      int iz = integralInfo->monomial_info.monomial_list[monomialIndex].iz;

      //printf("monomialIndex = %5i, ix, iy, iz = %2i %2i %2i\n", monomialIndex, ix, iy, iz);

      // Get x y z 1-d Hermite polynomials
      poly_1d_struct hermitePoly_1d_x;
      poly_1d_struct hermitePoly_1d_y;
      poly_1d_struct hermitePoly_1d_z;
      if(inverseFlag == 1)
	{
	  get_1d_hermite_poly_inv(&hermitePoly_1d_x, ix, 1.0/exponent);
	  get_1d_hermite_poly_inv(&hermitePoly_1d_y, iy, 1.0/exponent);
	  get_1d_hermite_poly_inv(&hermitePoly_1d_z, iz, 1.0/exponent);
	}
      else
	{
	  get_1d_hermite_poly(&hermitePoly_1d_x, ix, exponent);
	  get_1d_hermite_poly(&hermitePoly_1d_y, iy, exponent);
	  get_1d_hermite_poly(&hermitePoly_1d_z, iz, exponent);
	}

      // Store x y z Hermite polys as 3-d polys
      poly_3d_struct hermitePoly_3d_x;
      poly_3d_struct hermitePoly_3d_y;
      poly_3d_struct hermitePoly_3d_z;
      create_3d_poly_from_1d_poly(&hermitePoly_3d_x, &hermitePoly_1d_x, 0);
      create_3d_poly_from_1d_poly(&hermitePoly_3d_y, &hermitePoly_1d_y, 1);
      create_3d_poly_from_1d_poly(&hermitePoly_3d_z, &hermitePoly_1d_z, 2);

      // Compute product
      poly_3d_struct hermitePoly_3d_xy;
      poly_3d_struct hermitePoly_3d_xyz;
      compute_product_of_3d_polys(&hermitePoly_3d_xy, &hermitePoly_3d_x, &hermitePoly_3d_y);
      compute_product_of_3d_polys(&hermitePoly_3d_xyz, &hermitePoly_3d_xy, &hermitePoly_3d_z);
      
      // Go through result product poly, for each term get monomialIndex and add 
      // coeff to final result at position given by monomialIndex.
      int i;
      for(i = 0; i < hermitePoly_3d_xyz.noOfTerms; i++)
	{
	  poly_3d_term_struct* currTerm = &hermitePoly_3d_xyz.termList[i];
	  // Get monomialIndex2
	  int ix = currTerm->monomialInts[0];
	  int iy = currTerm->monomialInts[1];
	  int iz = currTerm->monomialInts[2];
	  int monomialIndex2 = integralInfo->monomial_info.monomial_index_list[ix][iy][iz];
	  //printf("currTerm->coeff = %22.11f\n", currTerm->coeff);
	  result[monomialIndex * noOfMonomials + monomialIndex2] += currTerm->coeff;
	} // END FOR i
      
    } // END FOR monomialIndex

  return 0;
}


#endif
