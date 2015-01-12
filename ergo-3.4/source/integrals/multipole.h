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

#ifndef MULTIPOLE_HEADER
#define MULTIPOLE_HEADER


#include "realtype.h"
#include "integral_info.h"
#include "basisinfo.h"

#define MAX_MULTIPOLE_DEGREE 15
#define MAX_NO_OF_MOMENTS_PER_MULTIPOLE ((MAX_MULTIPOLE_DEGREE+1)*(MAX_MULTIPOLE_DEGREE+1))

#define MAX_MULTIPOLE_DEGREE_BASIC BASIS_FUNC_POLY_MAX_DEGREE
#define MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC ((MAX_MULTIPOLE_DEGREE_BASIC+1)*(MAX_MULTIPOLE_DEGREE_BASIC+1))


typedef struct
{
  ergo_real centerCoords[3];
  int degree;
  int noOfMoments;
  ergo_real momentList[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
  ergo_real maxAbsMomentList[MAX_MULTIPOLE_DEGREE+1];
  ergo_real euclideanNormList[MAX_MULTIPOLE_DEGREE+1];
} multipole_struct_large;

typedef struct
{
  ergo_real centerCoords[3];
  int degree;
  int noOfMoments;
  ergo_real momentList[MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC];
} multipole_struct_small;


int init_multipole_code();

int
compute_multipole_moments(const IntegralInfo& integralInfo,
			  const DistributionSpecStruct* distr,
			  multipole_struct_small* result);

class MMTranslator {
  static const int MMDP1 = MAX_MULTIPOLE_DEGREE+1;
  ergo_real *buffer_W_cc;
  ergo_real *buffer_W_cs;
  ergo_real *buffer_W_sc;
  ergo_real *buffer_W_ss;
 public:
  MMTranslator();
  ~MMTranslator();
  int getTranslationMatrix(ergo_real dx,
                           ergo_real dy,
                           ergo_real dz,
                           int l_1,
                           int l_2,
                           ergo_real* result_W) const;
};

class MMInteractor {
  static const int MMDP1 = MAX_MULTIPOLE_DEGREE+1;
  ergo_real *buffer_T_cc;
  ergo_real *buffer_T_cs;
  ergo_real *buffer_T_sc;
  ergo_real *buffer_T_ss;
 public:
  MMInteractor();
  ~MMInteractor();
  int getInteractionMatrix(ergo_real dx,
                           ergo_real dy,
                           ergo_real dz,
                           int l_1,
                           int l_2,
                           ergo_real* result_T);
};


int
setup_multipole_maxAbsMomentList(multipole_struct_large* multipole);

#endif
