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

#ifndef BOYSFUNCTION_HEADER
#define BOYSFUNCTION_HEADER

#include "realtype.h"
#include "polydegree.h"

/* We need Boys functions up to a degree 4 times the highest degree of basis functions. */
#define BOYS_N_MAX (BASIS_FUNC_POLY_MAX_DEGREE*4+1)

#define BOYS_TAB_DEGREE 12
#define BOYS_X_MAX 160.0
#define BOYS_NO_OF_INTERVALS 200

typedef struct
{
  ergo_real midx;
  ergo_real A[BOYS_TAB_DEGREE];
} BoysFuncIntervalStruct;

class BoysFunctionManager {
 private:
  BoysFuncIntervalStruct Boys_list[BOYS_N_MAX][BOYS_NO_OF_INTERVALS];
  int Boys_init_flag;
  ergo_real BoysFunction_pretabulated(int n, ergo_real x) const;
 public:
  BoysFunctionManager() : Boys_init_flag(0) { }
  void init();
  ergo_real BoysFunction(int n, ergo_real x) const;
  // Stuff needed for Chunks&Tasks usage
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};

#endif
