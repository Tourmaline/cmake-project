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

/** @file PuriStepInfo.h PuriStepInfo class
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2007
 *
 */
#ifndef MAT_PURISTEPINFO
#define MAT_PURISTEPINFO
#include <math.h>
#include <iomanip>
#include "PuriStepInfoDebug.h"
#define __ID__ "$Id$"
namespace mat {
  /* IDEA: Adjust threshold (looser) in later iteration if threshold in early
   * iteration was tighter than expected.
   */

  /** Contains information about the matrix before a purification step and 
   *  about the step.
   *  All info is for the truncated matrix in the current step.
   *  We use inheritance for test class so that empty base class 
   *  optimization can be used in case of an empty test class.
   */

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    class PuriStepInfo : public PuriStepInfoDebug<Treal, TdebugPolicy> {
  public:
    explicit PuriStepInfo(int nn = -1, int noc = -1, Treal eigvalConvCrit = 0.0)
      : n(nn),nocc(noc), traceX(-1.0), traceX2(-1.0), poly(-1), 
      chosenThresh(0.0), actualThresh(0.0), 
      estimatedStepsLeft(-1),
      eigInterval(-1.0,2.0), delta(-1.0), correctOccupation(0),
      XmX2EuclNorm(0.0, 10.0), eigVecPtr(0), 
      lumoWasComputed(0), homoWasComputed(0),
      n0(0.0), n1(0.0), 
      homo(0.0, 1.0), lumo(0.0, 1.0), eigConvCrit(eigvalConvCrit), 
      nnzX(0), nnzX2(0), eigAccLoss(0),
      timeThresh(0),  timeSquare(0),  timeXmX2Norm(0), timeTotal(0),
      timeXX2Write(0), timeXX2Read(0) {}
    
    ~PuriStepInfo() {
      delete eigVecPtr;
    } 
      bool converged() const {
	bool homoLumo_converged = (1 - homo.low() < eigConvCrit && 
				   lumo.upp() < eigConvCrit);
	bool XmX2norm_converged = getXmX2EuclNorm().upp() < eigConvCrit;
	// FIXME: Possibly, propagating XmX2norm between the
	// iterations can give more accurate values since this norm is not
	// computed accurately in each iteration.
	return homoLumo_converged || XmX2norm_converged;
      }
      
      inline void setChosenThresh(Treal const thr) {chosenThresh = thr;}
      inline Treal getChosenThresh() const { return chosenThresh;}
      inline void setActualThresh(Treal const thr) {actualThresh = thr;}
      inline Treal getActualThresh() const { return actualThresh;}
      inline void setEstimatedStepsLeft(int const stepsleft) {
	estimatedStepsLeft = stepsleft;
      }
      inline int getEstimatedStepsLeft() const { return estimatedStepsLeft;}
      inline void setTraceX(Treal const trX) {traceX = trX;}
      inline void setTraceX2(Treal const trX2) {traceX2 = trX2;}
      inline Treal getTraceX() const {return traceX;}
      inline Treal getTraceX2() const {return traceX2;}
      void setPoly();
      inline int getPoly() const {return poly;}
      /** Sets XmX2EuclNorm bounds  */
      inline void setXmX2EuclNorm(Interval<Treal> const XmX2eucl) {
	XmX2EuclNorm = XmX2eucl;
      }
      
      inline Interval<Treal> getXmX2EuclNorm() const {return XmX2EuclNorm;}
      /** Improves homo and lumo bounds if the new ones are better. 
       *  Uses XmX2EuclNorm if possible.
       */

      void setEigVecPtr(Tvector * eigVecPtr_) {
	delete eigVecPtr;
	eigVecPtr = eigVecPtr_;
      }
      Tvector const * getEigVecPtr() const {return eigVecPtr;}

      void improveHomoLumo(Interval<Treal> const homoInt, 
			   Interval<Treal> const lumoInt);
      inline Interval<Treal> const & getHomo() const {
	return homo;
      }
      inline Interval<Treal> const & getLumo() const {
	return lumo;
      }
      inline Interval<Treal> const & getEigInterval() const {
	return eigInterval;
      }
      void exchangeInfoWithNext(PuriStepInfo<Treal, Tvector, TdebugPolicy> & next);
      /** Set correct occ.
       */
      inline void setCorrectOccupation() {correctOccupation = 1;}
      inline int getCorrectOccupation() const {return correctOccupation;}
      Treal subspaceError() const;
      /** Improve eigenvalue bounds and delta if possible.
       *  Returns delta.
       */
      void improveEigInterval(Interval<Treal> const eInt);
      
      inline void setNnzX(size_t const nzX) {nnzX = nzX;}
      inline size_t getNnzX() const {return nnzX;}
      inline void setNnzX2(size_t const nzX2) {nnzX2 = nzX2;}
      inline size_t getNnzX2() const {return nnzX2;}
      void computeEigAccLoss();
      inline Treal getEigAccLoss() const {return eigAccLoss;}
      
      inline Treal getN0() const {return n0;}
      inline Treal getN1() const {return n1;}
      inline Treal getDelta() const {return delta;}
      

      inline void checkIntervals(const char* descriptionString) const {
	PuriStepInfoDebug<Treal, TdebugPolicy>::
	  checkIntervals(eigInterval, homo, lumo, XmX2EuclNorm, descriptionString);    
      }
      template<typename Tmatrix>
	inline void computeExactValues(Tmatrix const & X, Tmatrix const & X2) {
	PuriStepInfoDebug<Treal, TdebugPolicy>::
	  computeExactValues(X, X2, n, nocc);
      }
  
      /* Functions to set and get timings and mem usage info for the different steps: */
      void setMemUsageBeforeTrunc() {
	MemUsage::getMemUsage(memUsageBeforeTrunc);
      }
      void setMemUsageInXmX2Diff(MemUsage::Values & memUsage) {
	memUsageInXmX2Diff = memUsage;
      }
      void setTimeThresh(float t) {timeThresh = t;}
      void setTimeSquare(float t) {timeSquare = t;}
      void setTimeXmX2Norm(float t) {timeXmX2Norm = t;}
      void setTimeTotal(float t) {timeTotal = t;}
      void setTimeXX2Write(float t) {timeXX2Write = t;}
      void setTimeXX2Read(float t) {timeXX2Read = t;}

      MemUsage::Values getMemUsageBeforeTrunc() { return memUsageBeforeTrunc; }
      MemUsage::Values getMemUsageInXmX2Diff()  { return memUsageInXmX2Diff;  }
      float getTimeThresh() {return timeThresh;}
      float getTimeSquare() {return timeSquare;}
      float getTimeXmX2Norm() {return timeXmX2Norm;}
      float getTimeTotal() {return timeTotal;}
      float getTimeXX2Write() {return timeXX2Write;}
      float getTimeXX2Read() {return timeXX2Read;}
      
      inline bool homoIsAccuratelyKnown
	(Treal accuracyLimit /** HOMO estimation is considered to be accurate
			      *  if the error is smaller than this value.
			      */
	 ) const {
	return homo.length() < accuracyLimit;
      }
      inline bool lumoIsAccuratelyKnown
	(Treal accuracyLimit /** LUMO estimation is considered to be accurate
			      *  if the error is smaller than this value.
			      */
	 ) const {
	return lumo.length() < accuracyLimit;
      }

      inline bool getLumoWasComputed() {return lumoWasComputed;}
      inline bool getHomoWasComputed() {return homoWasComputed;}
      
  protected:
      /** Compute n0 and n1 
       *  Called by improveEigInterval
       */
      void computen0n1();
      
      int n; /**< System size */
      int nocc; /**< Number of occupied orbitals. */  
      /* FIXME?: Upper and lower bound on traces? */
      Treal traceX; /**< Trace of the matrix X */ 
      Treal traceX2; /**< Trace of the squared matrix X^2 */
      int poly; /**< Choice of polynomial 0 for x^2 and 1 for 2 * x - x^2 */
      Treal chosenThresh; /**< Chosen threshold value applied before step. */
      Treal actualThresh; /**< Actual threshold value applied before step. */
      int estimatedStepsLeft; /**< Estimated steps left in purification. 
			       *   Used to chose threshold.
			       *   -1 indicates no estimation possible.
 			       */
      Interval<Treal> eigInterval; /**< Interval containing the 
				    *   eigenvalue spectrum. */
      Treal delta; /**< Largest possible deviation from the [0 1] interval. */
      int correctOccupation; 
      /**< Takes the values:
       *    1 meaning all eigenvalues supposed to go to 1 are larger than 0.5 
       *      and  all eigenvalues supposed to go to 0 are smaller than 0.5
       *      for sure.
       *    0 otherwise
       */
      Interval<Treal> XmX2EuclNorm; /**< Interval containing the euclidean norm 
				     *   ||X-X^2||_2 before step. 
				     */

      Tvector * eigVecPtr; /**< Eigenvector possibly containing the homo or 
			    *   lumo eigenvector. 
			    */
      bool lumoWasComputed; /**< Flag indicating if lumo was computed. (Eigenvector exists.) */
      bool homoWasComputed; /**< Flag indicating if homo was computed. */
      Treal n0; /**< Lower bound on the number of eigenvalues in 
		 *   [lambdaMin, 0.5]. 
		 */
      Treal n1; /**< Lower bound on the number of eigenvalues in 
		 *   [0.5, 1 + delta].
		 */
    
      Interval<Treal> homo; /**< Interval containing the homo eigenvalue. */
      Interval<Treal> lumo; /**< Interval containing the lumo eigenvalue. */
      Treal eigConvCrit; /**< Tolerated deviation from 0 and 1 of converged 
			  *   eigenvalues. 
			  */
      size_t nnzX;  /**< Number of nonzeros in the matrix X. */
      size_t nnzX2; /**< Number of nonzeros in the matrix X2. */
    
      /** A probable upper bound of the accuracy that is lost in the 
       *  eigenvalues of X * X because of limited relative precision in the 
       *  storage of X. 
       */    
      Treal eigAccLoss;
      
      /* Variables for the time and mem usage of various operations */
      MemUsage::Values memUsageBeforeTrunc;
      MemUsage::Values memUsageInXmX2Diff;
      float timeThresh;
      float timeSquare;
      float timeXmX2Norm;
      float timeTotal;
      float timeXX2Write;
      float timeXX2Read;
  private:
  };

#if 1
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    std::ostream& operator<<(std::ostream& s, 
			     PuriStepInfo<Treal, Tvector, TdebugPolicy> const & psi) {
    s<<" Trace X: "<<psi.getTraceX()
     <<" Trace X2: "<<psi.getTraceX2()
     <<" poly: "<<psi.getPoly()
     <<" ||E||_2: "<<psi.getActualThresh()
     <<" Eiginterval: "<<psi.getEigInterval()
     <<" correctOccupation: "<<psi.getCorrectOccupation()
     <<std::endl
     <<" XmX2EuclNorm: "<<psi.getXmX2EuclNorm()
     <<" n0: "<<psi.getN0()
     <<" n1: "<<psi.getN1()
     <<" homo: "<<psi.getHomo()
     <<" lumo: "<<psi.getLumo()
     <<" nnzX: "<<psi.getNnzX()
     <<" nnzX2: "<<psi.getNnzX2();
    return s;
  }

#endif

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriStepInfo<Treal, Tvector, TdebugPolicy>::setPoly() {
    if (Interval<Treal>::intersect(homo,lumo).empty()) {
      ASSERTALWAYS(homo.low() > lumo.upp());
      if (homo.low() + lumo.upp() > 1)
	poly = 0; /* x*x */
      else
	poly = 1; /* 2*x - x*x */
    }
    else {
      if (template_blas_fabs(traceX2 - nocc) < template_blas_fabs(2 * traceX - traceX2 - nocc))
	poly = 0; /* x*x */
      else
	poly = 1; /* 2*x - x*x */
    }
  }
  
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriStepInfo<Treal, Tvector, TdebugPolicy>::
    improveHomoLumo(Interval<Treal> const homoInt, 
		    Interval<Treal> const lumoInt) {
    checkIntervals("PuriStepInfo::improveHomoLumo A0");
    homo.intersect(homoInt);
    checkIntervals("PuriStepInfo::improveHomoLumo A1");
    lumo.intersect(lumoInt);
    checkIntervals("PuriStepInfo::improveHomoLumo A2");
    ASSERTALWAYS(!homo.empty());
    ASSERTALWAYS(!lumo.empty());
    if (homo.low() > 0.5 && lumo.upp() < 0.5)
      this->setCorrectOccupation();
    
    if (correctOccupation && 1.0 - XmX2EuclNorm.upp() * 4.0 > 0) {
      Interval<Treal> tmp = 
	sqrtInt((XmX2EuclNorm * (Treal)(-4.0)) + (Treal)1.0);
      ASSERTALWAYS(tmp.length() > 0);

      ASSERTALWAYS(!homo.empty());
      ASSERTALWAYS(!lumo.empty());

      homo.intersect(Interval<Treal>((1.0 + tmp.low()) / 2.0, homo.upp()));

      ASSERTALWAYS(!homo.empty());
      ASSERTALWAYS(!lumo.empty());

      lumo.intersect(Interval<Treal>(lumo.low(), (1.0 - tmp.low()) / 2.0));
      checkIntervals("PuriStepInfo::improveHomoLumo B");
      ASSERTALWAYS(!homo.empty());
      ASSERTALWAYS(!lumo.empty());
      if (!eigInterval.cover((1 + template_blas_sqrt(1.0 + 4.0 * XmX2EuclNorm.low())) / 2) &&
	  !eigInterval.cover((1 - template_blas_sqrt(1.0 + 4.0 * XmX2EuclNorm.low())) / 2)){
	/* Either homo lies in homoTmp or lumo lies in lumoTmp. */
	Interval<Treal> homoTmp = 
	  (tmp + (Treal)1.0) / (Treal)2.0;
	Interval<Treal> lumoTmp = 
	  ((tmp * (Treal)(-1.0))  + (Treal)1.0) / (Treal)2.0;
	ASSERTALWAYS(!(Interval<Treal>::intersect(homo, homoTmp).empty() && 
		       Interval<Treal>::intersect(lumo, lumoTmp).empty()));
	if (Interval<Treal>::intersect(homo, homoTmp).empty()) {
	  // ok, lumo was computed in this iteration
	  if (eigVecPtr)
	    lumoWasComputed = 1;
	  lumo.intersect(lumoTmp);
	}
	if (Interval<Treal>::intersect(lumo, lumoTmp).empty()) {
	  // ok, homo was computed in this iteration
	  if (eigVecPtr)
	    homoWasComputed = 1;
	  homo.intersect(homoTmp);
	}
	checkIntervals("PuriStepInfo::improveHomoLumo C");
      }
    }
  }
  
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriStepInfo<Treal, Tvector, TdebugPolicy>::
    exchangeInfoWithNext(PuriStepInfo<Treal, Tvector, TdebugPolicy> & next) {
    Interval<Treal> zeroOneInt(0.0,1.0);

    next.checkIntervals("PuriStepInfo::exchangeInfoWithNext A");

    /* Improve homo/lumo and eig-bounds for next */
    Interval<Treal> homoForNext(homo);
    Interval<Treal> lumoForNext(lumo);
    Interval<Treal> eigIntForNext(eigInterval);
    homoForNext.puriStep(poly);
    lumoForNext.puriStep(poly);
    eigIntForNext.puriStep(poly);
    ASSERTALWAYS(!homoForNext.empty());
    ASSERTALWAYS(!lumoForNext.empty());
    ASSERTALWAYS(!eigIntForNext.empty());
    /* Increase intervals because relative precision in matrix-matrix 
     * multiplication can result in loss of accuracy in eigenvalues
     */
    homoForNext.increase(eigAccLoss);
    lumoForNext.increase(eigAccLoss);
    eigIntForNext.increase(eigAccLoss);
    /* Increase intervals because of thresholding */
    homoForNext.increase(next.actualThresh);
    lumoForNext.increase(next.actualThresh);
    homoForNext.intersect(zeroOneInt);
    lumoForNext.intersect(zeroOneInt);
    eigIntForNext.increase(next.actualThresh);
    next.improveEigInterval(eigIntForNext);
    next.improveHomoLumo(homoForNext, lumoForNext);
    
    /* Improve homo/lumo for this */
    /* FIXME: Consider improving also eigInterval from next
     * This could possibly only be done in one end of the interval since
     * for example information about negative eigenvalues is lost in case
     * of an x*x step. 
     */
    Interval<Treal> homoTmp(next.homo);
    Interval<Treal> lumoTmp(next.lumo);
    ASSERTALWAYS(!homoTmp.empty());
    ASSERTALWAYS(!lumoTmp.empty());
    /* Increase intervals because of thresholding. */
    homoTmp.increase(next.actualThresh);
    lumoTmp.increase(next.actualThresh);
    homoTmp.intersect(zeroOneInt);
    lumoTmp.intersect(zeroOneInt);
    homoTmp.invPuriStep(poly);
    lumoTmp.invPuriStep(poly);
    /* Increase intervals because relative precision in matrix-matrix 
     * multiplication can result in loss of accuracy in eigenvalues
     */
    homoTmp.increase(eigAccLoss);
    lumoTmp.increase(eigAccLoss);

    this->improveHomoLumo(homoTmp, lumoTmp);
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    Treal PuriStepInfo<Treal, Tvector, TdebugPolicy>::subspaceError() const {
    Interval<Treal> gap = Interval<Treal>(lumo.upp(),homo.low()); 
    if (actualThresh >= gap.length())
      return 1.0; /* 1.0 means no accuracy. */
    else {
      Treal error = actualThresh / (gap.length() - actualThresh); 
      return error < 1.0 ? error : (Treal)1.0;
    }
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriStepInfo<Treal, Tvector, TdebugPolicy>::
    improveEigInterval(Interval<Treal> const eInt) {
    eigInterval.intersect(eInt);
    Treal delta1 = eigInterval.upp() - 1;
    Treal delta0 = -eigInterval.low();
    delta = delta1 > delta0 ? delta1 : delta0;
    this->computen0n1();
  }
  

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriStepInfo<Treal, Tvector, TdebugPolicy>::computen0n1() {
    Treal beta = 0.5; 
    /* Increase beta if possible */
    if (XmX2EuclNorm.upp() < 1 / (Treal)4)
      beta = (1 + template_blas_sqrt(1 - 4 * XmX2EuclNorm.upp())) / 2;
    n1 = (traceX2 - delta * (1 - beta) * n - 
	  (1 - delta - beta) * traceX) / 
      ((1 + 2 * delta) * (delta + beta));
    n0 = (traceX2 + beta * (1 + delta) * n - 
	  (1 + delta + beta) * traceX) / 
      ((1 + 2 * delta) * (delta + beta));
    if (n1 > nocc -1 && 
	n0 > n - nocc - 1)
      correctOccupation = 1;
  }

  /** Computes a probable upper bound of the accuracy that is lost in 
   *  the eigenvalues of X * X because of limited relative precision in 
   *  the storage of X. 
   */
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriStepInfo<Treal, Tvector, TdebugPolicy>::computeEigAccLoss() {
    Treal nnzPerRowX = nnzX / (Treal)n;
    Treal maxAbsErrPerElX2 = getRelPrecision<Treal>() * nnzPerRowX;
    /*  mah = max(abs(h_ij)) \approx relPrec * (nnz(X) / n)
     *  e is the exact eigenvalue of X^2, e' is the eigenvalue of the 
     *  computed X^2.
     *  | e - e' | <= || H ||_2 <=  || H ||_F = sqrt( sum h_ij^2 ) <=
     *  sqrt( mah^2 * nnz(X2))
     */
    eigAccLoss = maxAbsErrPerElX2 * template_blas_sqrt((Treal)nnzX2);
    ASSERTALWAYS(eigAccLoss >= 0);
  }
  


} /* end namespace mat */
#undef __ID__
#endif
