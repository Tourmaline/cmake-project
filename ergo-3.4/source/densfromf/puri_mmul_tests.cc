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

#include <sstream>
#include <cmath>
#include "puri_mmul_tests.h"
#include "output.h"
#include "utilities.h"

enum approachType {standard, modified};

static void run_one_puri_mmul_test( symmMatrix const & F_ort, 
				    intervalType const & allEigsInterval,
				    const std::vector<int> & poly_choices,
				    const std::vector<ergo_real> & thresh_values,
				    symmMatrix const & D_ort_ref,
				    int const globalCounter,
				    approachType approach) {
  Util::TimeMeter timeMeterTot;
  int nSteps = poly_choices.size();
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "run_puri_mmul_tests, nSteps = %d.", nSteps);
  for(int i = 0; i < nSteps; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "poly_choices[%2d] = %d, thresh_values[%2d] = %8.4g", 
	      i, poly_choices[i], i, thresh_values[i]);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "run_puri_mmul_tests, rescalingand shifting to [0,1] ...");
  symmMatrix X(F_ort);
  ergo_real lmin = allEigsInterval.low();
  ergo_real lmax = allEigsInterval.upp();
  X.add_identity(-lmax);      /* Scale to [0, 1] interval and negate */
  X *= ((ergo_real)1.0 / (lmin - lmax));
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "run_puri_mmul_tests, now doing iterations...");
  double accumulatedTimeMmul = 0;
  symmMatrix X2(X);
  size_t nnzMax = 0;
  for(int i = 0; i < nSteps; i++) {
    ergo_real currThresh = thresh_values[i];
    // Now compute X2 and truncate it. This can be done using different approaches.
    if(approach == standard) {
      Util::TimeMeter timeMeterMmul;
      X2 = (ergo_real)1.0 * X * X;
      accumulatedTimeMmul += Util::TimeMeter::get_wall_seconds() - timeMeterMmul.get_start_time_wall_seconds();
      size_t nnz = X2.nnz();
      if(nnz > nnzMax)
	nnzMax = nnz;
      X2.thresh(currThresh, mat::euclNorm);
    }
    else if(approach == modified) {
      // Split X into two parts, A and B.
      symmMatrix A(X);
      symmMatrix B(X);
      symmMatrix A2(X);
      symmMatrix B2(X);
      A = X;
      A.eucl_thresh(sqrt(currThresh/2));
      B = X - A;
      {
	Util::TimeMeter timeMeterMmul;
	A2 = (ergo_real)1.0 * A * A;
	accumulatedTimeMmul += Util::TimeMeter::get_wall_seconds() - timeMeterMmul.get_start_time_wall_seconds();
      }
      B2 = 0;//(ergo_real)1.0 * B * B;
      normalMatrix AB;
      {
	Util::TimeMeter timeMeterMmul;
	AB = (ergo_real)1.0 * A * B;
	accumulatedTimeMmul += Util::TimeMeter::get_wall_seconds() - timeMeterMmul.get_start_time_wall_seconds();
      }
      normalMatrix BA;
      BA = transpose(AB);
      normalMatrix ABplusBA;
      ABplusBA = AB + BA;
      symmMatrix ABplusBAsymm(ABplusBA);
      X2 = A2 + B2;
      X2 += ABplusBAsymm;
      size_t nnz = X2.nnz();
      if(nnz > nnzMax)
	nnzMax = nnz;
      X2.thresh(currThresh/2, mat::euclNorm);
    }
    else {
      throw "Error: unimplemented approach in run_one_puri_mmul_test";
    }
    ergo_real frobDiffX2 = symmMatrix::frob_diff(X, X2);
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "frobDiffX2 = %8.4g", frobDiffX2);
    symmMatrix Xnew;
    if(poly_choices[i] == 0) { // X*X
      Xnew = X2;
    }
    else { // 2*X - X*X
      Xnew = (ergo_real)2.0*X;
      Xnew += (ergo_real)-1.0*X2;
    }
    X = Xnew;
  }
  ergo_real frobDiff = symmMatrix::frob_diff(X, D_ort_ref);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "frobDiff = %8.4g", frobDiff);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "run_one_puri_mmul_test: accumulatedTimeMmul = %12.6f wall seconds", accumulatedTimeMmul);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "run_one_puri_mmul_test: nnzMax = %12d", nnzMax);
  timeMeterTot.print(LOG_AREA_DENSFROMF, "run_one_puri_mmul_test");
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "End of run_one_puri_mmul_test.");
}

void run_puri_mmul_tests( symmMatrix const & F_ort, 
			  intervalType const & allEigsInterval,
			  const std::vector<int> & poly_choices,
			  const std::vector<ergo_real> & thresh_values,
			  symmMatrix const & D_ort_ref,
			  int const globalCounter) {
  run_one_puri_mmul_test( F_ort, 
			  allEigsInterval,
			  poly_choices,
			  thresh_values,
			  D_ort_ref,
			  globalCounter,
			  standard);
  run_one_puri_mmul_test( F_ort, 
			  allEigsInterval,
			  poly_choices,
			  thresh_values,
			  D_ort_ref,
			  globalCounter,
			  modified);
}

