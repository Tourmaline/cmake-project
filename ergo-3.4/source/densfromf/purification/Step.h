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

/** @file Step.h Step class
 *
 * @author Emanuel Rubensson
 * @date January 2011
 *
 */


#ifndef PUR_STEP
#define PUR_STEP
#include <limits>
namespace pur {
  
  template<typename Treal>
    struct Step {
      typedef Treal real;
      
      int poly;   /**< The step just taken, 0 for x*x and 1 for 2*x-x*x */
      real alpha; /**< Scaling parameter just before the step just taken. */
      mat::Interval<real> eig_homo; /**< Interval containing the HOMO 
				       eigenvalue of X. */
      mat::Interval<real> eig_lumo; /**< Interval containing the LUMO 
				       eigenvalue of X. */
      mat::Interval<real> eig_homo_orig; /**< Original interval. */
      mat::Interval<real> eig_lumo_orig; /**< Original interval. */
      real traceX;
      real traceX2;
      real chosen_thresh;
      real actual_thresh; /**< If the matrix has been truncated, the
			   *   error inflicted measured by the
			   *   spectral norm. Otherwise it should be
			   *   zero.
			   */
      size_t nnzX;
      size_t nnzX2;
      double wall_sec_thresh;
      double wall_sec_square;
      double wall_sec_XmX2norm;
      double wall_sec_total;
      Step() {}
      Step(int const poly, real const alpha, 
	   mat::Interval<real> eig_homo,
	   mat::Interval<real> eig_lumo,
	   real traceX,
	   real traceX2,
	   real chosen_thresh,
	   real actual_thresh,
	   size_t nnzX,
	   size_t nnzX2,
	   double wall_sec_thresh,
	   double wall_sec_square,
	   double wall_sec_XmX2norm,
	   double wall_sec_total) 
      : poly(poly), alpha(alpha), 
	eig_homo(eig_homo), eig_lumo(eig_lumo), 
	eig_homo_orig(eig_homo), eig_lumo_orig(eig_lumo), 
	traceX(traceX), traceX2(traceX2),
	chosen_thresh(chosen_thresh),
	actual_thresh(actual_thresh),
        nnzX(nnzX), nnzX2(nnzX2), 
	wall_sec_thresh(wall_sec_thresh),
	wall_sec_square(wall_sec_square),
	wall_sec_XmX2norm(wall_sec_XmX2norm),
	wall_sec_total(wall_sec_total) {}
      
      void propagate_homo_to_previous(Step<real> & previous) const;
      void propagate_lumo_to_previous(Step<real> & previous) const;
    }; // end struct Step  

  template<typename Treal>
    void Step<Treal>::
    propagate_homo_to_previous( Step<real> & previous ) const {
    mat::Interval<real> zeroOneInt(0.0,1.0);
    mat::Interval<real> homo = eig_homo;
    homo.increase( actual_thresh );
    homo.intersect_always_non_empty(zeroOneInt);
    homo.invPuriStep(poly, alpha);
    previous.eig_homo.intersect_always_non_empty(homo);
  }
  template<typename Treal>
    void Step<Treal>::
    propagate_lumo_to_previous( Step<real> & previous ) const {
    mat::Interval<real> zeroOneInt(0.0,1.0);
    mat::Interval<real> lumo = eig_lumo;
    lumo.increase( actual_thresh );
    lumo.intersect_always_non_empty(zeroOneInt);
    lumo.invPuriStep(poly, alpha);    
    previous.eig_lumo.intersect_always_non_empty(lumo);
  }

  
} // end namespace pur
#endif
