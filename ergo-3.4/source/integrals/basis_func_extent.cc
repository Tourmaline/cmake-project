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

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "basis_func_extent.h"
#include "output.h"
#include "integrals_general.h"
#include "pi.h"
#include "integrals_2el_single.h"
#include "exponent_list.h"


static ergo_real
get_M(const IntegralInfo & integralInfo,
      const BasisInfoStruct & basisInfo)
{
  const JK::ExchWeights CAM_params_not_used;
  int n = basisInfo.noOfBasisFuncs;
  ergo_real M = 0;
  for(int i = 0; i < n; i++)
    {
      BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
      // go through all primitives for this basis function.
      int nPrims = basisFunc->noOfSimplePrimitives;
      int start  = basisFunc->simplePrimitiveIndex;
      int j;
      for(j = 0; j < nPrims; j++)
	{
	  DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
	  ergo_real currValue = do_2e_integral_using_symb_info(CAM_params_not_used, prim, prim, integralInfo);
	  if(currValue > M)
	    M = currValue;
	} // END FOR j
    } // END FOR i
  return M;
}


static int
compute_extent_for_all_basis_funcs_core(const BasisInfoStruct & basisInfo, 
					ergo_real* basisFuncExtentList, 
					ergo_real threshold,
					ExponentList exponentList,
					ergo_real M,
					ergo_real maxAbsDensMatElement)
{
  ergo_real twotopow1o4 = std::pow(2, 0.25);
  ergo_real pitopow5o4 = std::pow(pi, 1.25);
  
  int n = basisInfo.noOfBasisFuncs;
  for(int i = 0; i < n; i++)
    {
      BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
      ergo_real largestExtentSoFar = 0;
      // go through all primitives for this basis function.
      int nPrims = basisFunc->noOfSimplePrimitives;
      int start  = basisFunc->simplePrimitiveIndex;
      for(int j = 0; j < nPrims; j++)
	{
	  DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
	  ergo_real currExponent = prim->exponent;
	  ergo_real currAbsCoeff = std::fabs(prim->coeff);
	  ergo_real a = currExponent;
	  ergo_real c_a = currAbsCoeff;
	  // now go through all available exponents
	  for(int ii = 0; ii < exponentList.noOfExponents; ii++)
	    {
	      ergo_real b = exponentList.list[ii].exponent;
	      ergo_real c_b = exponentList.list[ii].maxAbsCoeff;
	      if(c_b > 0)
		{
		  ergo_real A = M * twotopow1o4 * pitopow5o4 * c_a * c_b * std::pow(a+b, (ergo_real)-1.25);
		  ergo_real R = std::sqrt( ((a + b) / (a * b)) * std::log(maxAbsDensMatElement * A / threshold));
		  if(R > largestExtentSoFar)
		    largestExtentSoFar = R;
		}
	    } // END FOR ii
	} // END FOR j
      basisFuncExtentList[i] = largestExtentSoFar;
    } // END FOR i
  return 0;
}


int
compute_extent_for_all_basis_funcs_2el(const IntegralInfo & integralInfo,
				       const BasisInfoStruct & basisInfo, 
				       ergo_real* basisFuncExtentList, 
				       ergo_real threshold,
				       ergo_real maxAbsDensMatElement)
{
  // Compute M = max sqrt[ (cc|cc) ]
  ergo_real M = get_M(integralInfo, basisInfo);

  // Create discretized list of exponents with maxAbsCoeff for each unique exponent.
  ExponentList exponentList;
  if(exponentList.get_list_of_available_exponents(basisInfo) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_available_exponents");
      return -1;
    }

  // Compute extent of each basis func by taking worst case of all available exponents.
  return compute_extent_for_all_basis_funcs_core(basisInfo, 
						 basisFuncExtentList, 
						 threshold,
						 exponentList,
						 M,
						 maxAbsDensMatElement);
}

