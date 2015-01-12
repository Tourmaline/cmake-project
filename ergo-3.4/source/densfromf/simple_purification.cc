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
#include "simple_purification.h"
#include "output.h"
#include "utilities.h"

typedef mat::Interval<ergo_real> intervalType;

void simple_purification( symmMatrix const & F_ort, 
			  int const noOfOccupiedOrbs,
			  intervalType const & allEigsInterval,
			  ergo_real const threshold,
			  ergo_real const conv_crit_frob_norm_diff,
			  int const maxiter,
			  symmMatrix & X,
			  int & iter,
			  std::vector<double> & wall_seconds,
			  std::vector<size_t> & nnz_X,
			  std::vector<size_t> & nnz_X2) {
  Util::TimeMeter timeMeterFirst;
  X = F_ort;
  symmMatrix X2(F_ort);
  ergo_real lmin = allEigsInterval.low();
  ergo_real lmax = allEigsInterval.upp();
  X.add_identity(-lmax);      /* Scale to [0, 1] interval and negate */
  X *= ((ergo_real)1.0 / (lmin - lmax));
  X.simple_blockwise_frob_thresh(threshold);
  X2 = (ergo_real)1.0 * X * X;
  
  std::vector<double> e(maxiter);

  nnz_X[0]  = X.nnz();
  nnz_X2[0] = X2.nnz();
  
  ergo_real trace_X  = X.trace();
  ergo_real trace_X2 = X2.trace();

  e[0] = symmMatrix::frob_diff(X,X2);

  wall_seconds[0] = timeMeterFirst.get_wall_seconds() - timeMeterFirst.get_start_time_wall_seconds();
  
  /*
    It is difficult to set the convergence criterion.  If we set it
    based on the desired accuracy there is a risk of never converging.

    If we quit "too early" the error may be dominated by eigenvalue errors.

    It is difficult to relate the convergence criterion to the threshold value
    
    Possible criterions include:
    Tr(X-X2)   < eps
    Tr(X-nocc) < eps
    ||X_i-X_{i-1}|| < eps
  */
  iter = 1; // iter == number of multiplications

  while(1) {
    // Check if it is time to break
    if(conv_crit_frob_norm_diff > 0) {
      if(symmMatrix::frob_diff(X,X2) < conv_crit_frob_norm_diff || iter >= maxiter )
	break;
    }
    else {
      if( ! ( iter < 1 || (( e[iter-1] > 0.1 || e[iter-1] < e[iter-2]) && iter < maxiter) ) )
	break;
    }
    Util::TimeMeter timeMeterIter;
    iter++;
    if ( std::fabs(noOfOccupiedOrbs - (2*trace_X-trace_X2)) < 
	 std::fabs(noOfOccupiedOrbs - trace_X2            )   ) {
      /*  Polynomial 2 * x - x * x  */
      X2 *= (ergo_real)-1.0;
      X2 += (ergo_real)2.0 * X;
      // Now X2 contains 2*X-X*X
    }
    X2.transfer(X); /* Transfer result to X by swapping pointers.
		     * This clears X2. */
    X.simple_blockwise_frob_thresh(threshold);
    X2 = (ergo_real)1.0 * X * X;
    trace_X  = X.trace();
    trace_X2 = X2.trace();    
    
    wall_seconds[iter-1] = timeMeterIter.get_wall_seconds() - timeMeterIter.get_start_time_wall_seconds();
    nnz_X [iter-1] = X.nnz();
    nnz_X2[iter-1] = X2.nnz();
    e[iter-1] = symmMatrix::frob_diff(X,X2);

  } // end while
}


static void get_E_norm_and_sin_theta(symmMatrix const & X, 
				     int const noOfOccupiedOrbs, 
				     intervalType const & allEigsInterval, 
				     symmMatrix const & D_exact, 
				     ergo_real threshold, 
				     double & E_norm,
				     double & sin_theta) {
  symmMatrix X_truncated(X);
  X_truncated.simple_blockwise_frob_thresh(threshold);
  // Get E_norm
  E_norm = (double)symmMatrix::eucl_diff(X, X_truncated, 1e-10);
  // Get sinTheta
  symmMatrix projector;
  int nmul;
  int maxiter = 200;
  std::vector<double> wallSecs(maxiter);
  std::vector<size_t> nnz_X(maxiter), nnz_X2(maxiter);
  simple_purification( X_truncated, 
		       noOfOccupiedOrbs,
		       allEigsInterval,
		       1e-10,
		       1e-10,
		       200,
		       projector,
		       nmul,
		       wallSecs,
		       nnz_X,
		       nnz_X2 );
  // Now we have the projector corresponding to X_truncated.
  sin_theta = (double)symmMatrix::eucl_diff(projector, D_exact, 1e-10);
  // Done!
}


static void output_m_vector(std::vector<double> const & v, std::string const & s, int n, std::ofstream & ff) {
  ff << s << " = [";
  for (int i=0; i<n; i++) 
    ff << v[i] << "  ";
  ff << "];\n";
}


void run_comparison_to_simple_purification
( symmMatrix const & F_ort, 
  int const noOfOccupiedOrbs,
  intervalType const & allEigsInterval,
  symmMatrix const & D_ort,
  int const globalCounter ) {
  char ffname[888];
  sprintf(ffname, "simple_puri_comparison_%i.m", globalCounter);
  std::ofstream ff(ffname);

  int maxiter = 200;
  // First compute "exact" density matrix
  symmMatrix D_exact;
  {
    ergo_real threshold_simple_blocked_frob = 0;
    ergo_real conv_crit_frob_norm_diff      = 1e-10;
    int nmul;
    std::vector<double> wall_seconds_all_iters(maxiter);
    std::vector<size_t> nnz_X(maxiter);
    std::vector<size_t> nnz_X2(maxiter);
    simple_purification( F_ort, 
			 noOfOccupiedOrbs,
			 allEigsInterval,
			 threshold_simple_blocked_frob,
			 conv_crit_frob_norm_diff,
			 maxiter,
			 D_exact,
			 nmul,
			 wall_seconds_all_iters,
			 nnz_X,
			 nnz_X2 );
    ergo_real D_eucl_diff = symmMatrix::eucl_diff(D_ort,D_exact,1e-10);
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, 
	      "Comparison to exact matrix, eucl_norm_error:  %22.11f.", 
	      (double)D_eucl_diff);
    ff << "D_from_homolumopuri_eucl_diff = " << D_eucl_diff  << ";\n";
  }
 
  symmMatrix D_simple;
  int const n_thresh = 6;
  ergo_real threshold_simple_blocked_frob[n_thresh] = {1e-2,1e-3,1e-4,1e-5,1e-6,1e-7};
  double wall_seconds_tot[n_thresh];
  double eucl_norm_error[n_thresh];
  double eucl_norm_idempotency_error[n_thresh];
  int nmul[n_thresh];
  symmMatrix D_simple_squared;
  std::vector< std::vector<double> > wall_seconds_matrix(n_thresh);
  std::vector< std::vector<size_t> > nnz_X_matrix(n_thresh);
  std::vector< std::vector<size_t> > nnz_X2_matrix(n_thresh);
  std::vector<double> E_norm_F_ort(n_thresh);
  std::vector<double> sin_theta_F_ort(n_thresh);
  std::vector<double> E_norm_D_ort(n_thresh);
  std::vector<double> sin_theta_D_ort(n_thresh);

  for (int i=0; i<n_thresh; i++) {
    //    ergo_real conv_crit_frob_norm_diff = 10*threshold_simple_blocked_frob[i];
    
    Util::TimeMeter timeMeterSimplePurification;
    wall_seconds_matrix[i].resize(maxiter);
    nnz_X_matrix[i].resize(maxiter);
    nnz_X2_matrix[i].resize(maxiter);
    simple_purification( F_ort, 
			 noOfOccupiedOrbs,
			 allEigsInterval,
			 threshold_simple_blocked_frob[i],
			 0, // do not use explicit conv crit value
			 maxiter,
			 D_simple,
			 nmul[i],
			 wall_seconds_matrix[i],
			 nnz_X_matrix[i],
			 nnz_X2_matrix[i] );
    wall_seconds_tot[i]    = timeMeterSimplePurification.get_wall_seconds() - timeMeterSimplePurification.get_start_time_wall_seconds();
    eucl_norm_error[i] = (double)symmMatrix::eucl_diff(D_exact,D_simple,1e-10);
    D_simple_squared = (ergo_real)1.0 * D_simple * D_simple;
    eucl_norm_idempotency_error[i] = (double)symmMatrix::eucl_diff(D_simple_squared,D_simple,1e-10);

    // Also get norm of error matrix and sinTheta for F_ort and D_ort.
    ergo_real currThresh = threshold_simple_blocked_frob[i];
    // First for F_ort.
    get_E_norm_and_sin_theta(F_ort, 
			     noOfOccupiedOrbs, 
			     allEigsInterval, 
			     D_exact, 
			     currThresh, 
			     E_norm_F_ort[i],
			     sin_theta_F_ort[i]);
    // Now for D_ort.
    symmMatrix D_exact_minus(D_exact);
    D_exact_minus *= (ergo_real)-1;
    intervalType allEigsIntervalForD_exact_minus(-1, 0);
    get_E_norm_and_sin_theta(D_exact_minus, 
			     noOfOccupiedOrbs, 
			     allEigsIntervalForD_exact_minus, 
			     D_exact, 
			     currThresh, 
			     E_norm_D_ort[i],
			     sin_theta_D_ort[i]);
  }

  ff << "threshold_simple_blocked_frob = [";
  for (int i=0; i<n_thresh; i++) 
    ff << threshold_simple_blocked_frob[i] << "  ";
  ff << "];\n";
  ff << "eucl_norm_error = [";
  for (int i=0; i<n_thresh; i++) 
    ff << eucl_norm_error[i] << "  ";
  ff << "];\n";
  ff << "eucl_norm_idempotency_error = [";
  for (int i=0; i<n_thresh; i++) 
    ff << eucl_norm_idempotency_error[i] << "  ";
  ff << "];\n";
  ff << "wall_seconds = [";
  for (int i=0; i<n_thresh; i++) 
    ff << wall_seconds_tot[i] << "  ";
  ff << "];\n";
  ff << "n_multiplications = [";
  for (int i=0; i<n_thresh; i++) 
    ff << nmul[i] << "  ";
  ff << "];\n";

  output_m_vector(E_norm_F_ort, "E_norm_F_ort", n_thresh, ff);
  output_m_vector(sin_theta_F_ort, "sin_theta_F_ort", n_thresh, ff);
  output_m_vector(E_norm_D_ort, "E_norm_D_ort", n_thresh, ff);
  output_m_vector(sin_theta_D_ort, "sin_theta_D_ort", n_thresh, ff);

  ff << "wall_seconds_matrix = [";
  for (int i=0; i<n_thresh; i++) {
    for(int j = 0; j < maxiter; j++) {
      double wallsecs = 0;
      if( j < nmul[i] )
	wallsecs = wall_seconds_matrix[i][j];
      ff << wallsecs << "  ";
    }
    ff << "\n";
  }
  ff << "];\n";
  ff << "nnz_X = [";
  for (int i=0; i<n_thresh; i++) {
    for(int j = 0; j < maxiter; j++) {
      size_t nnz = 0;
      if( j < nmul[i] )
	nnz = nnz_X_matrix[i][j];
      ff << nnz << "  ";
    }
    ff << "\n";
  }
  ff << "];\n";
  ff << "nnz_X2 = [";
  for (int i=0; i<n_thresh; i++) {
    for(int j = 0; j < maxiter; j++) {
      size_t nnz = 0;
      if( j < nmul[i] )
	nnz = nnz_X2_matrix[i][j];
      ff << nnz << "  ";
    }
    ff << "\n";
  }
  ff << "];\n";    
  
  ff.close();
}
