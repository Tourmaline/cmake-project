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

/** @file PuriStepInfoDebug.h PuriStepInfo base class for tests
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2007
 *
 */
#ifndef MAT_PURISTEPINFODEBUG
#define MAT_PURISTEPINFODEBUG
#include "matInclude.h"
#include "Interval.h"
#define __ID__ "$Id$"
namespace mat {
  template<typename Treal, typename TdebugPolicy>
    class PuriStepInfoDebug : public TdebugPolicy {
  public:
    
    inline void checkIntervals(Interval<Treal> const & eigInterval, 
			       Interval<Treal> const & homo, 
			       Interval<Treal> const & lumo, 
			       Interval<Treal> const & XmX2EuclNorm,
			       const char* descriptionString) const {}
    template<typename Tmatrix>
      inline void computeExactValues(Tmatrix const & X, Tmatrix const & X2,
				     int const n, int const nocc) {}
  protected:
    PuriStepInfoDebug(){}
  private:    
  };
  
  

/* Specialization for high debug level */
  template<typename Treal>
    class PuriStepInfoDebug<Treal, DebugLevelHigh> : public DebugLevelHigh {
  public:
    void checkIntervals(Interval<Treal> const & eigInterval, 
			Interval<Treal> const & homo, 
			Interval<Treal> const & lumo, 
			Interval<Treal> const & XmX2EuclNorm,
			const char* descriptionString) const;    
    template<typename Tmatrix>
      void computeExactValues(Tmatrix const & X, Tmatrix const & X2,
			      int const n, int const nocc);  
  protected:
    PuriStepInfoDebug()
      :homoExact(), lumoExact(), 
      lmaxExact(), lminExact(),
      XmX2EuclNormExact()
	{}
      Interval<Treal> homoExact; 
      Interval<Treal> lumoExact; 
      Interval<Treal> lmaxExact;
      Interval<Treal> lminExact;
      Interval<Treal> XmX2EuclNormExact;
  private:
  };

  template<typename Treal>
    void PuriStepInfoDebug<Treal, DebugLevelHigh>::
    checkIntervals(Interval<Treal> const & eigInterval, 
		   Interval<Treal> const & homo, 
		   Interval<Treal> const & lumo, 
		   Interval<Treal> const & XmX2EuclNorm,
		   const char* descriptionString) const {
    /*  Failure here probably means that checkIntervals was called before
     *  the exact values were calculated.
     */
    if(lminExact.empty())
      std::cout << "PuriStepInfoDebug::checkIntervals failed, descriptionString = '" 
		<< descriptionString << "'" << std::endl;
    ASSERTALWAYS(!lminExact.empty()); 
    Interval<Treal> zeroOneInt(0.0,1.0); 
    if (Interval<Treal>::intersect(eigInterval, lminExact).empty())
      std::cout<<std::endl<<" eigInterval "<<
	std::setprecision(25)<<eigInterval<<std::endl
	       <<" lminExact "<<std::setprecision(25)
	       <<lminExact<<std::endl;
    ASSERTDEBUG(!Interval<Treal>::intersect(eigInterval, lminExact).empty());
    if (Interval<Treal>::intersect(eigInterval, lmaxExact).empty())
      std::cout<<std::endl<<" eigInterval "<<
	std::setprecision(25)<<eigInterval<<std::endl
	       <<" lmaxExact "<<lmaxExact<<std::endl;
    ASSERTDEBUG(!Interval<Treal>::intersect(eigInterval, lmaxExact).empty());
    /* The homo and lumo intervals only provides information if 
     * the exact values lie in [0, 1] otherwise no check is done.
     */
    if (!Interval<Treal>::intersect(zeroOneInt, homoExact).empty()) {
      if (Interval<Treal>::intersect(homo, homoExact).empty())
	std::cout<<std::endl<<" homo "<<
	  std::setprecision(25)<<homo<<std::endl
		 <<" homoExact "<<homoExact<<std::endl;
      if(Interval<Treal>::intersect(homo, homoExact).empty())
        std::cout << "PuriStepInfoDebug::checkIntervals failed due to intersect(homo, homoExact) , descriptionString = '"
                  << descriptionString << "'" << std::endl;
      ASSERTDEBUG(!Interval<Treal>::intersect(homo, homoExact).empty());
    }
    if (!Interval<Treal>::intersect(zeroOneInt, lumoExact).empty()) {
      if (Interval<Treal>::intersect(lumo, lumoExact).empty())
	std::cout<<std::endl<<" lumo "<<
	  std::setprecision(25)<<lumo<<std::endl
		 <<" lumoExact "<<lumoExact<<std::endl;
      ASSERTDEBUG(!Interval<Treal>::intersect(lumo, lumoExact).empty());
    }
    if (Interval<Treal>::
	intersect(XmX2EuclNorm, XmX2EuclNormExact).empty()) {
      std::cout<<std::endl<<" XmX2EuclNorm "<<
	std::setprecision(25)<<XmX2EuclNorm<<std::endl
	       <<" XmX2EuclNormExact "<<XmX2EuclNormExact<<std::endl;
      std::cout << "PuriStepInfoDebug::checkIntervals failed, descriptionString = '" 
		<< descriptionString << "'" << std::endl;
    }
    ASSERTDEBUG(!Interval<Treal>::
		intersect(XmX2EuclNorm, XmX2EuclNormExact).empty());   
  }

  template<typename Treal>
    template<typename Tmatrix>
    void PuriStepInfoDebug<Treal, DebugLevelHigh>::
    computeExactValues(Tmatrix const & X, Tmatrix const & X2,
		       int const n, int const nocc) {
    /* Set up work space */
    std::vector<Treal> full(n*n);
    //    Treal* full = new Treal[n*n];
    Treal* eigvals = new Treal[n];
    int lwork = 3*n-1;
    Treal* work = new Treal[lwork];
    int info = 0;
    Treal euclNormMatrix;
    Treal precSyev;
    /* Compute lumoExact, homoExact, lminExact, and lmaxExact.  */
    X.fullMatrix(full);
    syev("N", "U", &n, &full[0], &n, eigvals, work, &lwork, &info);
    ASSERTALWAYS(!info);
    euclNormMatrix = template_blas_fabs(eigvals[0]) > template_blas_fabs(eigvals[n-1]) ? 
      template_blas_fabs(eigvals[0]) : template_blas_fabs(eigvals[n-1]);
    precSyev = euclNormMatrix * getRelPrecision<Treal>();
    lumoExact = Interval<Treal>(eigvals[n-nocc-1] - precSyev, 
				eigvals[n-nocc-1] + precSyev);
    homoExact = Interval<Treal>(eigvals[n-nocc] - precSyev, 
				eigvals[n-nocc] + precSyev);
    lminExact = Interval<Treal>(eigvals[0] - precSyev, 
				eigvals[0] + precSyev);
    lmaxExact = Interval<Treal>(eigvals[n-1] - precSyev, 
				eigvals[n-1] + precSyev);
    /* Compute largest magnitude eigenvalue of X-X^2 */
    Tmatrix XmX2(X);
    XmX2 += -1.0 * X2;
    XmX2.fullMatrix(full);
    syev("N", "U", &n, &full[0], &n, eigvals, work, &lwork, &info);
    ASSERTALWAYS(!info);
    euclNormMatrix = template_blas_fabs(eigvals[0]) > template_blas_fabs(eigvals[n-1]) ? 
      template_blas_fabs(eigvals[0]) : template_blas_fabs(eigvals[n-1]);
    precSyev = euclNormMatrix * getRelPrecision<Treal>();

#if 0
    std::cout<<"EIGS XmX2: "<<std::endl;
    for(int ind = 0; ind < n; ++ind)
      std::cout<<std::setprecision(20)<<eigvals[ind]<<std::endl;
    std::cout<<"end EIGS XmX2: "<<std::endl<<std::endl;
#endif
    
    Treal lowBound = euclNormMatrix - precSyev > 0 ? 
      euclNormMatrix - precSyev : 0;
    XmX2EuclNormExact = Interval<Treal>(lowBound, euclNormMatrix + precSyev);
    /* Free memory */
    delete[] eigvals;
    delete[] work;
  }


} /* end namespace mat */
#undef __ID__
#endif
