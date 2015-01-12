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

/** @file Purification_scaled.h Purification_scaled class
 *
 * @author Emanuel Rubensson
 * @date January 2011
 *
 */


#ifndef PUR_PURIFICATION_SCALED
#define PUR_PURIFICATION_SCALED
#include <iomanip>
#include <sys/time.h>
#include <stdexcept>
#include "matrix_proxy.h"
#include "matInclude.h"
#include "Step.h"
namespace pur {
  template<typename Tmatrix>
    class Purification_scaled {
  public:
    typedef typename Tmatrix::real real;

    Purification_scaled(Tmatrix & F_and_D,  /**< (input/output) 
					     *	  F on input, 
					     *	  D on output */ 
			mat::Interval<real> const & eigFInt,
			/**< (input) Interval containing all eigenvalues. */
			mat::Interval<real> const & hoF, 
			/**< (input) Interval containing the HOMO eigenvalue. */
			mat::Interval<real> const & luF, 
			/**< (input) Interval containing the LUMO eigenvalue. */
			real const toleratedEigenvalError, 
			/**< Tolerated deviation of eigenvalues from 0 and 1 
			     in the final result. */
			real const toleratedSubspaceError,
			/**< Tolerated error in the occupied subspace  
			     in the final result (measured by sin(theta)). */
			int const max_steps, /**< Max nr of iterations. */
			mat::normType normForTruncation, 
			/**< Norm to use for truncation. */
			bool const use_scaling
			/**< Use scaling technique to speed up calculation. */
			); 


    void purify();

    void get_homo_lumo_intervals( mat::Interval<real> & hoF,
				  mat::Interval<real> & luF );
    /**< Computed eigenvalues of F. Can only be called after
     *   convergence. 
     */

    void mInfo(std::ostream & file) const;
    
    void mTime(std::ostream & file) const;

    class Time;
  private:
    
    static int get_poly( real const homo, real const lumo ) {
      //      return fabs(homo*homo + lumo*lumo - 1) > 
      //	fabs(2*homo - homo*homo + 2*lumo - lumo*lumo - 1);
      return homo + lumo < 1; 
    }

    static int estimated_steps_left( mat::Interval<real> eig_homo, 
				     mat::Interval<real> eig_lumo,
				     real const toleratedEigenvalError,
				     int const max_steps,
				     bool const use_scaling ); 
    /**< Estimates number of remaining iterations.  Possible movements
     *   of eigenvalues due to truncation is not taken into account.
     *   The computed number is actually the number of truncations
     *   left and since the function is always called just before a
     *   truncation, it is always larger than or equal to 1.
     */
    
    static real get_threshold( real const gap, 
			       int const steps_left, 
			       real const acc_error,
			       real const toleratedSubspaceError );
    /**< Compute acceptable error due to truncation in an iteration
     *   based on the current band gap, the number of steps left, and
     *   the total tolerated error in the occupied subspace (
     *   tolSubspaceError ) as measured by the largest canonical angle
     *   between exact and perturbed subspaces. Information about the
     *   total error accumulated so far is also used.
     */

    void improve_homo_lumo_based_on_normXmX2( mat::Interval<real>& homo,
					      mat::Interval<real>& lumo,
					      bool & homo_was_computed,
					      bool & lumo_was_computed );
    /**< Computes interval containing spectral norm of X-X^2 matrix.
     *   Computed accurately if HOMO or LUMO eigenvalue has not
     *   previously been computed and is the eigenvalue closest to
     *   0.5, and if the norm is small (which it usually is only in
     *   later iterations when there is larger separation between
     *   eigenvalues).  The HOMO and LUMO eigenvalues may be improved.
     */
    real accumulated_error(); 
    /**< Total accumulated error due to removal of small matrix
     *   elements up to current iteration.  This value is (cheaply)
     *   recalculated in each iteration since information about
     *   HOMO/LUMO eigenvalues may have changed.
     */

    void propagate_homo_information();
    /**< Go through step vector and improve homo information in each
     *   step possible.
     */
    void propagate_lumo_information();
    /**< Go through step vector and improve lumo information in each
     *   step possible.
     */

    bool converged();

    Tmatrix & X;
    Tmatrix   X2;
    int n;
    std::vector<Step<real> > puri_steps; 
    mat::Interval<real> const & eigFInt;
    real const tolEigenvalError; 
    /**< Tolerated error in eigenvalues at convergence.  */
    real const tolSubspaceError; 
    /**< Tolerated subspace error. */
    unsigned int current_step; 
    mat::normType normTruncation;    
    bool use_scaling;
    bool homo_is_computed;
    bool lumo_is_computed;
  };

  template<typename Tmatrix>
    Purification_scaled<Tmatrix>::
    Purification_scaled( Tmatrix & F_and_D,
			 mat::Interval<real> const & eigFInt,
			 mat::Interval<real> const & hoF, 
			 mat::Interval<real> const & luF, 
			 real const toleratedEigenvalError, 
			 real const toleratedSubspaceError,
			 int const max_steps, 
			 mat::normType normForTruncation,
			 bool const use_scaling) 
    : X(F_and_D),  
    puri_steps(max_steps+1), 
    eigFInt(eigFInt),
    tolEigenvalError(toleratedEigenvalError), 
    tolSubspaceError(toleratedSubspaceError),
    current_step(0),
    normTruncation(normForTruncation),
    use_scaling(use_scaling),
    homo_is_computed(0),
    lumo_is_computed(0) {
      Time t_tot; t_tot.tic();
      n = X.get_nrows();
      using namespace mat;
      real lmin = eigFInt.low();
      real lmax = eigFInt.upp();
      X.add_identity(-lmax);      /* Scale to [0, 1] interval and negate */
      X *= ((real)1.0 / (lmin - lmax));
      /* Compute transformed homo and lumo eigenvalues. */
      mat::Interval<real> homo = hoF;
      mat::Interval<real> lumo = luF;
      // homo and lumo must be in the [lmin, lmax] interval
      homo.intersect( mat::Interval<real>(lmin,lmax) );
      lumo.intersect( mat::Interval<real>(lmin,lmax) );
      if ( homo.empty() )
	throw "homo interval empty in Purification_scaled constructor due to incorrect input.";
      if ( lumo.empty() )
	throw "lumo interval empty in Purification_scaled constructor due to incorrect input.";
      homo = (homo - lmax) / (lmin - lmax);
      lumo = (lumo - lmax) / (lmin - lmax);
      
      int steps_left     = estimated_steps_left( homo, lumo,
						 tolEigenvalError,
						 2*puri_steps.size(),
						 use_scaling );
      real acc_error     = accumulated_error(); 
      real chosen_thresh = get_threshold( homo.low() - lumo.upp(), 
					  steps_left, acc_error,
					  tolSubspaceError );
      Time t_thr; t_thr.tic();
      real actual_thresh = X.thresh( chosen_thresh, normTruncation );
      double wall_sec_thresh = t_thr.toc();
      // Eigenvalues may have moved due to truncation. Therefore, we
      // increase the intervals containing the HOMO and LUMO
      // eigenvalues according to Weyl's theorem.
      homo.increase( actual_thresh );
      lumo.increase( actual_thresh );
      // We do not believe that homo or lumo eigenvalues are outside
      // [0, 1] here
      mat::Interval<real> zero_one_int(0.0,1.0);
      homo.intersect(zero_one_int);
      lumo.intersect(zero_one_int);
      real const ONE = 1.0;
      Time t_sqr; t_sqr.tic();
      X2 = ONE * X * X;
      double wall_sec_square = t_sqr.toc();
      double wall_sec_total = t_tot.toc();
      puri_steps[current_step] = Step<real>( -1, 0, 
					     homo, lumo, 
					     X.trace(), X2.trace(),
					     chosen_thresh,
					     actual_thresh,
					     X.nnz(), X2.nnz(),
					     wall_sec_thresh,
					     wall_sec_square,
					     0,
					     wall_sec_total);      
    } 
  
  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::purify() {
    using namespace mat;
    while ( !converged() && current_step < puri_steps.size()-1 ) {
      Time t_tot; t_tot.tic();
      mat::Interval<real> homo = puri_steps[current_step].eig_homo;
      mat::Interval<real> lumo = puri_steps[current_step].eig_lumo;

      int poly = get_poly( homo.low(), lumo.upp() );
      // Scaling
      real alpha = 1; // alpha == 1 means no scaling 
      if (use_scaling) {
	if (poly) {
	  alpha = 2 / ( 1+homo.upp() );
	}
	else {
	  alpha = 2 / ( 2-lumo.low() );
	}      
      }
      if (poly) {
	/*  Polynomial 2 * alpha * x - alpha^2 * x^2  */
	X2 *= (real)-(alpha*alpha);
        X2 += ((real)2.0 * alpha) * X;
      }
      else {
	/*  Polynomial ( alpha * x + (1-alpha) )^2  */
	X2 *= (real)(alpha*alpha);
        X2 += ((real)2.0 * alpha * (1-alpha)) * X;
	X2.add_identity( (1-alpha)*(1-alpha) );
      }
      // Apply polynomial to homo and lumo intervals
      homo.puriStep(poly, alpha);
      lumo.puriStep(poly, alpha);
      // Transfer content of X2 to X clearing previous content of X if any
      // In current implementation this is needed regardless of which
      // polynomial is used
      X2.transfer(X); 

      int steps_left     = estimated_steps_left( homo, lumo,
						 tolEigenvalError,
						 2*puri_steps.size(),
						 use_scaling );
      real acc_error     = accumulated_error();
      assert(acc_error >= 0);
      real chosen_thresh = get_threshold( homo.low() - lumo.upp(), 
					  steps_left, acc_error,
					  tolSubspaceError );
      Time t_thr; t_thr.tic();
      real actual_thresh = X.thresh( chosen_thresh, normTruncation );
      double wall_sec_thresh = t_thr.toc();
      // Eigenvalues may have moved due to truncation. Therefore, we
      // increase the intervals containing the HOMO and LUMO
      // eigenvalues according to Weyl's theorem.
      homo.increase( actual_thresh );
      lumo.increase( actual_thresh );
      
      real const ONE = 1.0;
      Time t_sqr; t_sqr.tic();
      X2 = ONE * X * X;
      double wall_sec_square = t_sqr.toc();
      bool homo_was_computed = false;
      bool lumo_was_computed = false;
      Time t_norm; t_norm.tic();
      improve_homo_lumo_based_on_normXmX2( homo, lumo, 
					   homo_was_computed,
					   lumo_was_computed );
      double wall_sec_XmX2norm = t_norm.toc();
      ++current_step;
      double wall_sec_total = t_tot.toc();
      puri_steps[current_step] = Step<real>( poly, alpha, 
					     homo, lumo, 
					     X.trace(), X2.trace(),
					     chosen_thresh,
					     actual_thresh,
					     X.nnz(), X2.nnz(),
					     wall_sec_thresh,
					     wall_sec_square,
					     wall_sec_XmX2norm,
					     wall_sec_total);      
      if ( homo_was_computed ) {
	propagate_homo_information();
	homo_is_computed = true;
      }
      if ( lumo_was_computed ) {
	propagate_lumo_information();
	lumo_is_computed = true;
      }
    } // end while
  } // end purify()
  
  template<typename Tmatrix>
    int Purification_scaled<Tmatrix>::
    estimated_steps_left( mat::Interval<real> homo, 
			  mat::Interval<real> lumo,
			  real const toleratedEigenvalError,
			  int const max_steps,
			  bool const use_scaling) {

    int steps = 1; 
    // There is always at least one truncation left to do and to be on
    // the safe side we expect that truncation will result in an extra
    // step which means that the estimated number of steps is always
    // larger than or equal to 2.
    while ((1 - homo.low() > toleratedEigenvalError && 
	    lumo.upp() > toleratedEigenvalError) && steps < max_steps) {
      steps++;
      int poly = get_poly( homo.low(), lumo.upp() );
      real alpha = 1;
      if (use_scaling) {
	if (poly) {
	  alpha = 2 / ( 1+homo.upp() );
	}
	else {
	  alpha = 2 / ( 2-lumo.low() );
	}      
      }
      homo.puriStep(poly, alpha);
      lumo.puriStep(poly, alpha);
    } // end while
    return steps;
  } 
  
  template<typename Tmatrix>
    typename Tmatrix::real Purification_scaled<Tmatrix>::
    get_threshold( real const gap, 
		   int const steps_left, 
		   real const acc_error,
		   real const toleratedSubspaceError ) {
    real acceptable_error = toleratedSubspaceError - acc_error;
    /* 
      Set number of steps left to be at least smallest_factor since
      then there will always be some room left for another
      truncation. That is, we can make some error in the last
      iterations in case the number of steps was underestimated. The
      number smallest_factor can be replaced but should always be
      larger than 1. This number can be used to control the
      convergence of homo and lumo eigenvalues in the final
      iterations.
    */
    real smallest_factor = 1.01;
    real steps_left_pessimistic = steps_left > smallest_factor ? steps_left : smallest_factor;
    assert(acceptable_error >= 0);
    return ( gap * acceptable_error/steps_left_pessimistic ) / 
             ( 1 + acceptable_error/steps_left_pessimistic );
  }

#if 1
  // The following new version differs from the previous one in that
  // the norm ||X-X^2|| is computed accurately not only if we are sure to
  // compute a wanted eigenvalue but also if there is a chance that we
  // do so. Sometimes, if the computed eigenvalue is only in one of
  // the previously known homo and lumo eigenvalue intervals, one can
  // after the computation know which eigenvalue was
  // computed. However, sometimes it is not possible to know even after
  // computation which eigenvalue was computed, and then the
  // computation was done in vain. In the previous purification
  // version of 2008 (or earlier?) the information (the norm value)
  // was stored since further information may become available in
  // later iterations. That approach is not used here.
  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::
    improve_homo_lumo_based_on_normXmX2( mat::Interval<real>& homo,
					 mat::Interval<real>& lumo,
					 bool & homo_was_computed,
					 bool & lumo_was_computed ) {
    homo_was_computed = false;
    lumo_was_computed = false;
    bool computeAccurately = false;
    bool homo_may_be_computed = false;
    bool lumo_may_be_computed = false;
    if ( homo.low() > 0.5 && homo.low() + lumo.low() < 1 ) {
      homo_may_be_computed = true;
      // Necessary but not sufficient conditions to compute homo
      // eigenvalue satisfied
      if ( !homo_is_computed ) 
	computeAccurately = true;
    }
    if ( lumo.upp() < 0.5 && homo.upp() + lumo.upp() > 1 ) {
      lumo_may_be_computed = true;
      // Necessary but not sufficient conditions to compute lumo
      // eigenvalue satisfied
      if ( !lumo_is_computed ) 
	computeAccurately = true;
    }
    //    real const XmX2ENIsSmallValue = 0.207106781186547;
    // A value here of 0.2475 means that the interval [0.45, 0.55]
    // must be empty from eigenvalues. 
    // (0.45-0.45*0.45 == 0.55-0.55*0.55 == 0.2475)
    real const XmX2ENIsSmallValue = 0.2475; 
    real diffAcc = puri_steps[current_step].chosen_thresh / 100;
    mat::Interval<real> normXmX2;
    if ( computeAccurately ) 
      normXmX2 = Tmatrix::diffIfSmall(X, X2,  
				      mat::euclNorm, diffAcc, 
				      XmX2ENIsSmallValue);
    else
      normXmX2 = Tmatrix::diffIfSmall(X, X2,  
				      mat::frobNorm, diffAcc, 
				      XmX2ENIsSmallValue);
#if 0
    // Code for computing upper bound using mixed norm. Not clear if
    // this would have any effect.
    real norm_upper_bound = Tmatrix::mixed_diff(X, X2, diffAcc);
    norm_upper_bound = norm_upper_bound < normXmX2.upp() ? norm_upper_bound : normXmX2.upp();
    normXmX2 = mat::Interval<real>(normXmX2.low(), norm_upper_bound);
#endif
    // The norm normXmX2 must be smaller than or equal to 0.25
    {
      real lowerBound = normXmX2.low();
      real upperBound = normXmX2.upp() < (real)0.25 ? normXmX2.upp() : (real)0.25; 
      normXmX2 = mat::Interval<real>( lowerBound , upperBound );
    }
    
    // Try to improve homo or lumo
    mat::Interval<real> tmp = mat::sqrtInt((normXmX2 * (real)(-4.0)) + (real)1.0);
    // Let's see if we can see if we have computed the homo or lumo
    // eigenvalue (or perhaps none of them)
    //
    // if we have computed homo and homo > 1/2  =>
    mat::Interval<real> homoTmp = (tmp + (real)1.0) / (real)2.0;
    // if we have computed lumo and lumo < 1/2  =>
    mat::Interval<real> lumoTmp = 
      ((tmp * (real)(-1.0))  + (real)1.0) / (real)2.0;
    if (homo_may_be_computed && !lumo.overlap(lumoTmp) ) {
      homo.intersect_always_non_empty(homoTmp);
      if ( computeAccurately && normXmX2.length() <= 2.1*diffAcc ) 
	homo_was_computed = true;            
    }
    if (lumo_may_be_computed  && !homo.overlap(homoTmp)) {
      lumo.intersect_always_non_empty(lumoTmp);
      if ( computeAccurately && normXmX2.length() <= 2.1*diffAcc ) 
	lumo_was_computed = true;      
    }
  } // end new version of improve_homo_lumo_based_on_normXmX2

#else
  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::
    improve_homo_lumo_based_on_normXmX2( mat::Interval<real>& homo,
					 mat::Interval<real>& lumo,
					 bool & homo_was_computed,
					 bool & lumo_was_computed ) {
    homo_was_computed = false;
    lumo_was_computed = false;
    bool computeAccurately = false;
    bool homoIsComputable = false;
    bool lumoIsComputable = false;
    if ( homo.upp() + lumo.upp() < 1 ) {
      // Possibly homo can be computed
      // Check if homo can be and needs to be computed
      if (homo.low() > 0.5) {
	homoIsComputable = true;
	if ( !homo_is_computed ) 
	  computeAccurately = true;
      }
    }
    if ( homo.low() + lumo.low() > 1 ){
      // Possibly lumo can be computed
      // Check if lumo can be and needs to be computed
      if (lumo.upp() < 0.5) {
	lumoIsComputable = true;
	if ( !lumo_is_computed ) 
	  computeAccurately = true; 
      }     
    }

    //    real const XmX2ENIsSmallValue = 0.207106781186547;
    // A value here of 0.2475 means that the interval [0.45, 0.55]
    // must be empty from eigenvalues. 
    // (0.45-0.45*0.45 == 0.55-0.55*0.55 == 0.2475)
    real const XmX2ENIsSmallValue = 0.2475; 
    real diffAcc = puri_steps[current_step].chosen_thresh / 100;
    mat::Interval<real> normXmX2;
    if ( computeAccurately ) 
      normXmX2 = Tmatrix::diffIfSmall(X, X2,  
				      mat::euclNorm, diffAcc, 
				      XmX2ENIsSmallValue);
    else
      normXmX2 = Tmatrix::diffIfSmall(X, X2,  
				      mat::frobNorm, diffAcc, 
				      XmX2ENIsSmallValue);
    // The norm normXmX2 must be smaller than or equal to 0.25
    {
      real lowerBound = normXmX2.low();
      real upperBound = normXmX2.upp() < (real)0.25 ? normXmX2.upp() : (real)0.25; 
      normXmX2 = mat::Interval<real>( lowerBound , upperBound );
    }
    // Try to improve homo or lumo
    mat::Interval<real> tmp = mat::sqrtInt((normXmX2 * (real)(-4.0)) + (real)1.0);
    if (homoIsComputable) {
      // homo > 1/2  =>
      mat::Interval<real> homoTmp = (tmp + (real)1.0) / (real)2.0;
      homo.intersect_always_non_empty(homoTmp);
      if ( computeAccurately && normXmX2.length() <= 2.1*diffAcc ) 
	homo_was_computed = true;      
    }
    if (lumoIsComputable) {
      // lumo < 1/2  =>
      mat::Interval<real> lumoTmp = 
	((tmp * (real)(-1.0))  + (real)1.0) / (real)2.0;
      lumo.intersect_always_non_empty(lumoTmp);
      if ( computeAccurately && normXmX2.length() <= 2.1*diffAcc )  
	lumo_was_computed = true;      
    }
  } // end get_normXmX2
#endif

  template<typename Tmatrix>
    typename Tmatrix::real Purification_scaled<Tmatrix>::accumulated_error() {
    real total_error = 0;
    for (unsigned int step = 0; step <= current_step; step++) {
      real normE = puri_steps[step].actual_thresh;
      real gap   = puri_steps[step].eig_homo.low() - puri_steps[step].eig_lumo.upp();
      total_error += normE / (gap - normE);
    }
    return total_error;
  }

  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::propagate_homo_information() {
    for (unsigned int step = current_step; step >= 1; step-- ) 
      puri_steps[step].propagate_homo_to_previous( puri_steps[step-1] );
  }
  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::propagate_lumo_information() {
    for (unsigned int step = current_step; step >= 1; step-- ) 
      puri_steps[step].propagate_lumo_to_previous( puri_steps[step-1] );
  }
  
  template<typename Tmatrix>
    bool Purification_scaled<Tmatrix>::converged() {
    mat::Interval<real> homo = puri_steps[current_step].eig_homo;
    mat::Interval<real> lumo = puri_steps[current_step].eig_lumo;
    return ( 1 - homo.low() <= tolEigenvalError && 
	     lumo.upp() <= tolEigenvalError );
  }

  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::get_homo_lumo_intervals
    (mat::Interval<real> & hoF,
     mat::Interval<real> & luF) {
    if ( !converged() ) 
      throw "Attempt to get homo and lumo eigenvalue intervals when purification has not converged.";
    mat::Interval<real> homo = puri_steps[0].eig_homo;
    mat::Interval<real> lumo = puri_steps[0].eig_lumo;
    homo.increase( puri_steps[0].actual_thresh );
    lumo.increase( puri_steps[0].actual_thresh );
    real lmin = eigFInt.low();
    real lmax = eigFInt.upp();
    hoF = homo * (lmin - lmax) + lmax;
    luF = lumo * (lmin - lmax) + lmax;
  }

  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::
    mInfo(std::ostream & s) const {
    s<<"% PURIFICATION INFO IN MATLAB/OCTAVE FILE FORMAT"<<std::endl;
    s<<"% Norm for truncation: "<< mat::getNormTypeString(normTruncation) 
     << std::endl;
    s<<"n = "<<n<<";"<<std::endl;
    s<<"tolSubspaceError = "<<tolSubspaceError<<";"<<std::endl;
    s<<"tolEigenvalError = "<<tolEigenvalError<<";"<<std::endl;
    s<<"nIter = "<<current_step+1<<";"<<std::endl;

    s<<std::setprecision(15);
    s<<"% traceX and traceX2 \n";
    s<<"traces = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s<< puri_steps[ind].traceX  << "  "
       << puri_steps[ind].traceX2 << std::endl;
    }
    s<<"];\n";
    s<<"% poly \n";
    s<<"poly = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s<< puri_steps[ind].poly << std::endl;
    }
    s<<"];\n";
    s<<"% chosenThresh actualThresh \n";
    s<<"thresh = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s<< puri_steps[ind].chosen_thresh << "  "
       << puri_steps[ind].actual_thresh << std::endl;
    }
    s<<"];\n";
    s<<"% nnzX nnzX2 \n";
    s<<"nnzVals = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s<< puri_steps[ind].nnzX << "  "
       << puri_steps[ind].nnzX2 << std::endl;
    }
    s<<"];\n";
    s<<"% lumo_low lumo_upp homo_low homo_upp  \n";
    s<<"homo_lumo_eigs = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s<< puri_steps[ind].eig_lumo_orig.low() << "  " 
       << puri_steps[ind].eig_lumo_orig.upp() << "  " 
       << puri_steps[ind].eig_homo_orig.low() << "  " 
       << puri_steps[ind].eig_homo_orig.upp() << std::endl;
    }
    s<<"];\n";

    s<<"figure; \n"
     <<"plot(1:nIter,(homo_lumo_eigs(:,3)+homo_lumo_eigs(:,4))/2,'-b',"
     <<     "1:nIter,(homo_lumo_eigs(:,1)+homo_lumo_eigs(:,2))/2,'-r')\n"
     <<"legend('HOMO','LUMO'),xlabel('Iteration')\n"
     <<"hold on\n"
     <<"for ind = 1:nIter \n"
     <<"  plot([ind ind],[homo_lumo_eigs(ind,3) homo_lumo_eigs(ind,4)], '-b')\n"
     <<"  plot([ind ind],[homo_lumo_eigs(ind,1) homo_lumo_eigs(ind,2)], '-r')\n"
     <<"end\n"
     <<"axis([0 nIter 0 1])\n";

    s<<"figure; \n"
     <<"subplot(211)\n"
     <<"plot(1:nIter,100*nnzVals(:,1)/(n*n),'o-r',...\n"
     <<"1:nIter, 100*nnzVals(:,2)/(n*n),'x-b')\n"
     <<"legend('nnz(X)','nnz(X^2)'),xlabel('Iteration') \n"
     <<"ylabel('Percentage')\n"
     <<"axis([0 nIter 0 100])\n";
    s<<"subplot(212)\n"
     <<"semilogy(1:nIter,thresh(:,1),'x-r',1:nIter,thresh(:,2),'o-b')\n"
     <<"xlabel('Iteration'), ylabel('Threshold')\n"
     <<"legend('Chosen threshold', 'Actual threshold')\n"
     <<"axis([0 nIter min(thresh(:,2))/10 max(thresh(:,1))*10])\n";
    
  }

  template<typename Tmatrix>
    void Purification_scaled<Tmatrix>::
    mTime(std::ostream & s) const {
    s<<"% PURIFICATION TIMINGS IN MATLAB/OCTAVE FILE FORMAT"<<std::endl;
    s<<"% Norm for truncation: "<< mat::getNormTypeString(normTruncation) 
     << std::endl;
    s<<"n = "<<n<<";"<<std::endl;
    s<<"tolSubspaceError = "<<tolSubspaceError<<";"<<std::endl;
    s<<"tolEigenvalError = "<<tolEigenvalError<<";"<<std::endl;
    s<<"wall_sec_thresh = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s << puri_steps[ind].wall_sec_thresh 
	<< "  " << std::endl;
    }
    s<<"];\n";
    s<<"wall_sec_square = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s << puri_steps[ind].wall_sec_square
	<< "  " << std::endl;
    }
    s<<"];\n";
    s<<"wall_sec_XmX2norm = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s << puri_steps[ind].wall_sec_XmX2norm
	<< "  " << std::endl;
    }
    s<<"];\n";
    s<<"wall_sec_total = [";
    for (unsigned int ind = 0; ind <= current_step; ++ind) {
      s << puri_steps[ind].wall_sec_total
	<< "  " << std::endl;
    }
    s<<"];\n";
    // Plotting
    s<<"puriTime = [wall_sec_square wall_sec_thresh wall_sec_XmX2norm];\n";
    s<<"figure; bar(puriTime(:,1:3),'stacked')"<<std::endl<<
      "legend('Matrix Square', 'Truncation', '||X-X^2||'),"<<
      " xlabel('Iteration'), ylabel('Time (seconds)')"<<std::endl;
    s<<"figure; plot(wall_sec_total,'-x')"<<std::endl;       
  }

  template<typename Tmatrix>
    class Purification_scaled<Tmatrix>::Time {
  public:
    Time() {}
    void tic() {ticTime = get_wall_seconds();}
    double toc() {return get_wall_seconds() - ticTime;}
  private:
    double ticTime;
    static double get_wall_seconds() {
      struct timeval tv;
      if(gettimeofday(&tv, NULL) != 0)
	throw std::runtime_error("Error in get_wall_seconds(), in gettimeofday().");
      double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
      return seconds;
    }    
  };


} /* end namespace pur */
#endif
