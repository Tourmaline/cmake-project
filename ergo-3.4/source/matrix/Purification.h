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

/** @file Purification.h Purification class
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2007
 *
 */
#ifndef MAT_PURIFICATION
#define MAT_PURIFICATION
#include <math.h>
#include <iomanip>
#define __ID__ "$Id$"
namespace mat {
  template<typename Treal, typename Tmatrix, typename TdebugPolicy>
    class Purification :public TdebugPolicy {
  public: 
    typedef typename Tmatrix::VectorType VectorType;
    Purification(Tmatrix& M, /**< Fock/Kohn-Sham matrix (input)
			      *   Density matrix        (output) 
			      */ 
		 normType const normXmX2, /**< Norm for calculation of norm of X-X2. */
		 PuriInfo<Treal, VectorType, TdebugPolicy>& info
		 /**< Purification information object. Works as input with
		  *   system size (n), number of occupied orbitals (nocc),
		  *   and number of allowed multiplications.
		  *   Outputs all kind of information about the purification
		  *   process, see PuriInfo.h.
		  */
		 ); 
    void step();    
    void purify();
  protected:
    Tmatrix & X;
    Tmatrix X2;
    normType const normTruncation;
    normType const normXmX2;
    PuriInfo<Treal, VectorType, TdebugPolicy>& info;
    int niter;
    
    void stepComputeInfo(PuriStepInfo<Treal, VectorType, TdebugPolicy> & currentStep);
    
  private:
  };
  
  template<typename Treal, typename Tmatrix, typename TdebugPolicy>
    Purification<Treal, Tmatrix, TdebugPolicy>::
    Purification(Tmatrix& M, 
		 normType const normXmX2_inp, 
		 PuriInfo<Treal, VectorType, TdebugPolicy>& infop)
    :X(M), X2(M), normTruncation( infop.getNormForTruncation() ), normXmX2(normXmX2_inp), info(infop),niter(0) {
    Time tStepTotal;
    tStepTotal.tic();

    Treal lmin(0);
    Treal lmax(0);
    PuriStepInfo<Treal, VectorType, TdebugPolicy> & currentStep = info.getNext();
    currentStep.improveEigInterval(Interval<Treal>(0.0,1.0));
    Interval<Treal> eigFInt = info.getEigFInterval();
    lmin = eigFInt.low();
    lmax = eigFInt.upp();
    X.add_identity(-lmax);      /* Scale to [0, 1] interval and negate */
    X *= ((Treal)1.0 / (lmin - lmax));
    /* Compute transformed homo and lumo eigenvalues. */
    Interval<Treal> homo = info.getHomoF();
    Interval<Treal> lumo = info.getLumoF();
    homo = (homo - lmax) / (lmin - lmax);
    lumo = (lumo - lmax) / (lmin - lmax);

    Treal chosenThresh = info.getOptimalThresh();
    currentStep.setChosenThresh(chosenThresh);
    /* Consider convergence of eigs in getOptimalThresh */
    currentStep.setMemUsageBeforeTrunc();
    Time t;
    t.tic();
    Treal actualThresh = X.thresh(chosenThresh, normTruncation);
    currentStep.setTimeThresh(t.toc());
    currentStep.setActualThresh(actualThresh);
    if (homo.empty() || lumo.empty()) {
      throw Failure("Purification<Treal, Tmatrix, TdebugPolicy>::"
		    "Purification(Tmatrix&, normType const, " 
		    "PuriInfo<Treal, VectorType, TdebugPolicy>&): "
		    "HOMO or LUMO empty in purification constructor. ");
    }
    
    homo.increase(actualThresh);
    lumo.increase(actualThresh);
    stepComputeInfo(currentStep);
    currentStep.improveHomoLumo(homo,lumo);
    currentStep.setTimeTotal(tStepTotal.toc());
  } /**< Constructor */
    
    template<typename Treal, typename Tmatrix, typename TdebugPolicy>
      void Purification<Treal, Tmatrix, TdebugPolicy>::step() {
      Time tStepTotal;
      tStepTotal.tic();
      PuriStepInfo<Treal, VectorType, TdebugPolicy> & currentStep = info.getNext();
      /* It may happen that X2 has many more nonzeros than X, for
	 example 5 times as many.  Therefore it makes sense to try
	 having only one "big" matrix in memory at a time. However,
	 file operations have proved to be quite expensive and should
	 be avoided if possible. Hence we want to achieve having only
	 one big matrix in memory without unnecessary file
	 operations. We are currently hoping that it will be ok to add
	 a "small" matrix to a "big" one, that the memory usage after
	 that operation will be like the memory usage for one big
	 matrix + one small matrix. Therefore we are adding X to X2 (X
	 is truncated, a "small" matrix) instead of the opposite when
	 the 2*X-X*X polynomial is evaluated.
      */
      if (info(niter).getPoly()) {
	/*  Polynomial 2 * x - x * x  */
	X2 *= (Treal)-1.0;
	X2 += (Treal)2.0 * X;
	// Now X2 contains 2*X-X*X
      }
      /* In case of polynomial x * x we only need to transfer the
	 content of X2 to X.  */
      // Transfer content of X2 to X clearing previous content of X if any
      // In current implementation this is needed regardless of which
      // polynomial is used
      X2.transfer(X); 

     /* Obtain homo/lumo information from previous before thresh choice. 
       *  Default value of 0.0 for thresh is used. getOptimalThresh can
       *  try different thresh values.
       */
      Treal chosenThresh = info.getOptimalThresh();
      currentStep.setChosenThresh(chosenThresh);
      /* Consider convergence of eigs in getOptimalThresh? */
      currentStep.setMemUsageBeforeTrunc();
      Time t;
      t.tic();
      currentStep.setActualThresh(X.thresh(chosenThresh, normTruncation));
      currentStep.setTimeThresh(t.toc());
      stepComputeInfo(currentStep);
      if (niter > 5 &&
	  currentStep.getHomo().low() > 0.9 &&
	  currentStep.getLumo().upp() < 0.1 &&
	  info(niter-5).getHomo().low() > 0.9 &&
	  info(niter-5).getLumo().upp() < 0.1 &&
	  ((1 - currentStep.getHomo().low()) > 
	   (1 - info(niter-5).getHomo().low()) - 
	   currentStep.getEigAccLoss()  ||
	   currentStep.getLumo().upp() >
	   info(niter-5).getLumo().upp() - 
	   currentStep.getEigAccLoss())) {	
	throw Failure("Purification<Treal, Tmatrix, TdebugPolicy>"
		      "::step(): HOMO-LUMO gap has not increased in the "
		      "last five iterations. Desired accuracy too high for"
		      "current floating point precision?");
      }
      
      ++niter;
      currentStep.setTimeTotal(tStepTotal.toc());
    }
    
    template<typename Treal, typename Tmatrix, typename TdebugPolicy>
      void Purification<Treal, Tmatrix, TdebugPolicy>::purify() {
      
      while (niter < info.getMaxSteps() - 1 && !info.converged() ) {
	step();
      }
      if ( info.converged() ) {
	// Only if converged because forcing correct occupation can have 
	// strange effects otherwise.
	info.forceCorrectOccupation();	
	info.improveCorrectOccupation();
	info.improveInfo();
	info.improveHomoLumoF();
      }
    }

    template<typename Treal, typename Tmatrix, typename TdebugPolicy>
      void Purification<Treal, Tmatrix, TdebugPolicy>::
      stepComputeInfo(PuriStepInfo<Treal, VectorType, TdebugPolicy> & currentStep) {
      Treal const XmX2ENIsSmallValue = 0.207106781186547;
      Treal const ONE = 1.0;
      
      Time t;
      t.tic();
      X2 = ONE * X * X;
      currentStep.setTimeSquare(t.toc());
  
      currentStep.setNnzX(X.nnz());
      currentStep.setNnzX2(X2.nnz());
      currentStep.computeEigAccLoss();
      currentStep.computeExactValues(X, X2);
      currentStep.setTraceX(X.trace());
      currentStep.setTraceX2(X2.trace());

      /* Now we are about to compute the Euclidean norm of (X -
	 X2). Previously this operation needed lots of memory. Now,
	 however, we hope that the memory usage is smaller because the
	 difference matrix is never explicitly computed. */
      
      typename Tmatrix::VectorType * eigVecPtr = new typename Tmatrix::VectorType;
      Treal diffAcc;
      Interval<Treal> XmX2EN;
      t.tic();
      if (info.ShouldComputeXmX2EuclNormAccurately(diffAcc)) {
	if (normXmX2 == euclNorm)
	  XmX2EN = Tmatrix::euclDiffIfSmall(X, X2,  
					    diffAcc, 
					    XmX2ENIsSmallValue,
					    eigVecPtr);	  
	else
	  XmX2EN = Tmatrix::diffIfSmall(X, X2,  
					normXmX2, diffAcc, 
					XmX2ENIsSmallValue);
      }
      else {
	XmX2EN = Tmatrix::diffIfSmall(X, X2,  
				      frobNorm, diffAcc, 
				      XmX2ENIsSmallValue);
      }
      currentStep.setTimeXmX2Norm(t.toc());

      if (!eigVecPtr->is_empty())
	currentStep.setEigVecPtr(eigVecPtr);
      else
	delete eigVecPtr;

      XmX2EN.increase(currentStep.getEigAccLoss());
      Interval<Treal> zeroInt(0.0, XmX2EN.upp());
      XmX2EN.intersect(zeroInt);
      // FIXME: 
      // Consider improving lanczos so that if to good accuracy is requested
      // the best possible accuracy is returned
      currentStep.setXmX2EuclNorm(XmX2EN);
      currentStep.checkIntervals("Purification::stepComputeInfo after setXmX2EuclNorm.");

      Treal tmpVal = template_blas_sqrt(1 + 4 * XmX2EN.upp());
      currentStep.improveEigInterval
	(Interval<Treal>((1 - tmpVal) / 2, (1 + tmpVal) / 2));

      currentStep.checkIntervals("Purification::stepComputeInfo B");
      
      info.improveInfo();
      if (currentStep.getEigInterval().length() > 1.5)
	throw Failure("Purification<Treal, Tmatrix, TdebugPolicy>"
		      "::step() : "
		      "Eigenvalues moved to far from [0, 1] interval.");
      /* Now decide which polynomial to use for first step. */
      currentStep.setPoly();
      currentStep.checkIntervals("Purification::stepComputeInfo C");
    }




} /* end namespace mat */
#undef __ID__
#endif
