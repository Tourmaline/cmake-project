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

#ifndef INTEGRALS_2EL_HEADER
#define INTEGRALS_2EL_HEADER

#include "basisinfo.h"


namespace JK {

struct Params
{
  int use_naive_fockmatrix_construction;
  ergo_real threshold_J;
  ergo_real threshold_K;
  ergo_real multipole_threshold_factor;
  int use_densfit_for_J;
  int use_differential_density;
  int use_fmm;
  ergo_real fmm_box_size;
  int fmm_no_of_branches;
  ergo_real fmm_branch_splitter_extent_1;
  ergo_real fmm_branch_splitter_extent_2;
  ergo_real fmm_branch_splitter_extent_3;
  ergo_real fmm_branch_splitter_extent_4;
  ergo_real fmm_branch_splitter_extent_5;
  ergo_real exchange_box_size;
  int noOfThreads_J;
  int noOfThreads_K;

  Params() : use_naive_fockmatrix_construction(0),
       threshold_J(1e-12),
       threshold_K(1e-12),
       multipole_threshold_factor(1),
       use_densfit_for_J(0),
       use_differential_density(0),
       use_fmm(1),
       fmm_box_size(5.0),
       fmm_no_of_branches(0),
       fmm_branch_splitter_extent_1(0),
       fmm_branch_splitter_extent_2(0),
       fmm_branch_splitter_extent_3(0),
       fmm_branch_splitter_extent_4(0),
       fmm_branch_splitter_extent_5(0),
       exchange_box_size(5.0),
       noOfThreads_J(1),
       noOfThreads_K(1)
  {}
    
};

 

}


#if 0

typedef struct
{
  int shell_idx_1;
  int shell_idx_2;
  int noOfExponentPairs;
  int exponentPairIndex;
} shell_pair_struct;


typedef struct
{
  int idx_1;
  int idx_2;
  int productOrgIdx;
  ergo_real exponent_1;
  ergo_real exponent_2;
  ergo_real coeff_1;
  ergo_real coeff_2;
  ergo_real coeff_12, exponent_12;
  ergo_real sizeOfProduct;
  ergo_real centerCoords_1[3];
  ergo_real centerCoords_2[3];
  ergo_real centerCoords_12[3];
  int n1max;
  int n2max;
} exponent_pair_struct;

#endif


#endif
