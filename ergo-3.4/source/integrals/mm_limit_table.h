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

#ifndef MM_LIMIT_TABLE_HEADER
#define MM_LIMIT_TABLE_HEADER

#include "realtype.h"
#include "multipole.h"



void 
mm_limits_init(ergo_real maxDistance);


ergo_real 
mm_limits_get_max_abs_mm_contrib(int degree1,
				 const ergo_real* maxMomentVectorNormList1,
				 int degree2,
				 const ergo_real* maxMomentVectorNormList2,
				 ergo_real distance);


int 
mm_limits_get_minimum_multipole_degree_needed(ergo_real distance,
					      const multipole_struct_large* boxMultipole, 
					      int maxDegreeForDistrs, 
					      const ergo_real* maxMomentVectorNormForDistrsList, 
					      ergo_real threshold);



#endif
