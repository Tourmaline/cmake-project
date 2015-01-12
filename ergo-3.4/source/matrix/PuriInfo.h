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

/** @file PuriInfo.h PuriInfo class
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2007
 *
 */
#ifndef MAT_PURIINFO
#define MAT_PURIINFO
#include <math.h>
#include <iomanip>
#include "PuriStepInfo.h"
#define __ID__ "$Id$"
namespace mat {
  /** Contains information about a purification process
   *
   */
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    class PuriInfo : public TdebugPolicy {
  public:
    PuriInfo(int nn, int noc, 
	     Interval<Treal> eigFInt,
	     Interval<Treal> hoF,
	     Interval<Treal> luF,
	     Treal toleratedEigenvalError,
	     Treal toleratedSubspaceError,
	     normType normForTruncation_,
	     int maxS = 100)
      : n(nn),nocc(noc), step(new PuriStepInfo<Treal, Tvector, TdebugPolicy>[maxS]),
      maxSteps(maxS), nSteps(0), correct_occupation_was_forced_flag(false), 
      eigFInterval(eigFInt), homoF(hoF), lumoF(luF), 
      tolSubspaceError(toleratedSubspaceError),
      tolEigenvalError(toleratedEigenvalError), 
      normForTruncation(normForTruncation_) {
      for (int ind = 0; ind < maxSteps; ++ind)
	step[ind] = 
	  PuriStepInfo<Treal, Tvector, TdebugPolicy>(n, nocc, tolEigenvalError); 
    }
      virtual ~PuriInfo() {
	delete[] step;
      }
      /** Set the correctOccupation flag in the current step to 1. */
      void forceCorrectOccupation();
      /** Improves the correct occupation indicator 
       *  Call AFTER convergence and ONLY if it is known that
       *  the occupation is correct at convergence.
       */
      void improveCorrectOccupation();
      /** Improve homo / lumo values in each step.
       *  Call after call to improveCorrectOccupation()
       */
      void improveInfo();
      inline Interval<Treal> getEigFInterval() const {
	return eigFInterval;
      }
      inline PuriStepInfo<Treal, Tvector, TdebugPolicy> & getNext() {
	nSteps++;
	ASSERTALWAYS(nSteps - 1 < maxSteps);
	return step[nSteps - 1];
      }
      inline PuriStepInfo<Treal, Tvector, TdebugPolicy> & operator()(int const ind) {
	assert(ind >= 0);
	assert(ind < nSteps);
	return step[ind];
      }


      /** Returns the subspace error introduced so far. */
      inline Treal subspaceError() const {
	return subspaceError(nSteps);
      }
      /** Estimates the number of steps (upper bound) needed for 
       *  convergence based on homo/lumo information and desired accuracy
       *  in subspace and eigenvalues. Also computes optimal thresh-value 
       *  for the current step.
       */
      void estimateStepsLeft(int& stepsLeft, Treal& thresh) const;
      Treal getOptimalThresh() const;
      Treal getThreshIncreasingGap(Interval<Treal> const & middleGap) const;
      bool ShouldComputeXmX2EuclNormAccurately(Treal & howAccurate) const;
      /** Returns the best interval containing the homo eigenvalue
       *  (not transformed to [0, 1])
       */
      inline Interval<Treal> getHomoF() const {return homoF;}
      /** Returns the best interval containing the lumo eigenvalue
       *  (not transformed to [0, 1])
       */
      inline Interval<Treal> getLumoF() const {return lumoF;}
      inline int getMaxSteps() const {return maxSteps;}
      inline int getNSteps() const {return nSteps;}
      /* Returns a vector or 0 and 1 giving the used poly sequence. */
      void getPolys(std::vector<int> & resultVector) {
	resultVector.resize(nSteps);
	for(int i = 0; i < nSteps; i++)
	  resultVector[i] = step[i].getPoly();
      }
      /* Returns a vector of chosen threshold values. */
      void getThreshValues(std::vector<Treal> & resultVector) {
	resultVector.resize(nSteps);
	for(int i = 0; i < nSteps; i++)
	  resultVector[i] = step[i].getChosenThresh();
      }
      /* Tries to improve the homoF and lumoF eigenvalues.
       * Called after the purification process.
       */
      bool correct_occupation_was_forced() const 
      {return correct_occupation_was_forced_flag;}
      void improveHomoLumoF();

      inline bool converged() {return step[nSteps - 1].converged();}

      double getAccumulatedTimeSquare() const;
      double getAccumulatedTimeThresh() const;
      double getAccumulatedTimeXmX2Norm() const;
      double getAccumulatedTimeTotal() const;

      void mTimings(std::ostream & file) const;
      void mInfo(std::ostream & file) const;
      void mMemUsage(std::ostream & file) const;
      /** HOMO estimation is considered to be accurate
       *  if the error is smaller than tolSubspaceError / 100 in
       *  some purification step.
       */
      bool homoIsAccuratelyKnown() const {
	bool res = false;
	for(int ind = 0; ind < nSteps; ++ind)
	  res = res || step[ind].homoIsAccuratelyKnown(tolSubspaceError / 100);
	return res;
      }
      /** LUMO estimation is considered to be accurate
       *  if the error is smaller than tolSubspaceError / 100 in
       *  some purification step.
       */      
      bool lumoIsAccuratelyKnown() const {
	bool res = false;
	for(int ind = 0; ind < nSteps; ++ind)
	  res = res || step[ind].lumoIsAccuratelyKnown(tolSubspaceError / 100);
	return res;
      }

      normType getNormForTruncation() const {
	return normForTruncation;
      }

      void getHOMOandLUMOeigVecs(Tvector & eigVecLUMO, 
				 Tvector & eigVecHOMO) const;
      
  protected:
      int n; /**< System size */
      int nocc; /**< Number of occupied orbitals. */  

      PuriStepInfo<Treal, Tvector, TdebugPolicy>* step;
      int maxSteps; /**< Capacity of step array. */
      
      int nSteps; /**< Number of taken steps. Number of used elements in step 
		   *   array.
		   */
      bool correct_occupation_was_forced_flag; /**< Correct occupation was
						  assumed, not guaranteed. */
      Interval<Treal> const eigFInterval; 
      /**< Interval containing all eigenvalues before transformation 
       *   to the [0, 1] interval. 
       *   Also, these bounds will be used for the initial transformation.
       */
      Interval<Treal> homoF; 
      /**< Interval containing the HOMO eigenvalue before 
       *   transformation to [0, 1].
       *   Hopefully improved after purification. 
       */
      Interval<Treal> lumoF; 
      /**< Interval containing the LUMO eigenvalue before 
       *   transformation to [0, 1]. 
       *   Hopefully improved after purification. 
       */
      Treal const tolSubspaceError; 
      /**< Tolerated error in angle between correct and computed
       *   subspace.
       */
      Treal const tolEigenvalError; 
      /**< Tolerated error in eigenvalues at convergence.  */

      normType const normForTruncation;
      /**< Norm used for truncation of small matrix elements. */

      /** Returns the subspace error introduced until step end. */
      Treal subspaceError(int end) const;
      
  private:
  };
  
  
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::forceCorrectOccupation() {
    if ( step[nSteps-1].getCorrectOccupation() ) 
      return;
    step[nSteps-1].setCorrectOccupation();
    correct_occupation_was_forced_flag = true;
    return;
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::improveCorrectOccupation() {
    if (step[nSteps-1].getCorrectOccupation()) {
      Treal distance;
      Interval<Treal> middleInt;
      Interval<Treal> zeroOneInt(0.0,1.0);
      for (int ind = 0; ind < nSteps; ind++) {
	Treal XmX2Eucl = step[ind].getXmX2EuclNorm().upp();
	if ( XmX2Eucl < 1 / (Treal)4) {
	  distance = (1 - template_blas_sqrt(1 - 4 * XmX2Eucl)) / 2;
	  middleInt = Interval<Treal>(distance, 1.0 - distance); 
	  int i = ind;
	  while (!middleInt.empty() && i < nSteps - 1) {
	    middleInt.puriStep(step[i].getPoly());
	    middleInt.decrease(step[i+1].getActualThresh());
	    /* Make sure we stay in [0, 1] */
	    middleInt.intersect(zeroOneInt);
	    ++i;
	  } /* end while */
	  if (middleInt.cover(0.5))
	    step[ind].setCorrectOccupation();
	} /* end if */
      } /* end for */
    }
  }
  
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::improveInfo() {
    for (int ind = nSteps - 2; ind >= 0; ind--) {
      step[ind].exchangeInfoWithNext(step[ind + 1]);
    }
    for (int ind = 0; ind < nSteps - 1; ind++) {
      step[ind].exchangeInfoWithNext(step[ind + 1]);
    }  
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    Treal PuriInfo<Treal, Tvector, TdebugPolicy>::subspaceError(int end) const {
    Treal error = 0;
    for (int ind = 0; ind < end; ind++) {
      error += step[ind].subspaceError();
    }
    return error;
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::
    estimateStepsLeft(int& stepsLeft, Treal& thresh) const {
    stepsLeft = -1;
    thresh = 0;
    Interval<Treal> initialGap;
    Interval<Treal> gap;
    Treal tolError = tolSubspaceError - subspaceError(nSteps - 1);
    Treal initialAltThresh = 1;
    if (tolError <= 0.0) 
      return;
    
    /* Compute number of steps needed to converge. */
    /* Compute initial interval */
    /* nSteps == 0 means that a Purification object has not been 
     * associated to this instance yet.
     * nSteps == 1 means that a Purification object has been associated
     * to this instance but no steps have been taken yet.
     */
    Treal lastThresh = 0;
    if (nSteps == 0 || nSteps == 1) {
      Treal lmax = eigFInterval.upp();
      Treal lmin = eigFInterval.low();
      /* Compute transformed homo and lumo eigenvalues. */
      initialGap = Interval<Treal>((lumoF.low() - lmax) / (lmin - lmax),
				   (homoF.upp() - lmax) / (lmin - lmax));
    }
    else {
      initialGap = Interval<Treal>(step[nSteps - 2].getLumo().upp(),
				   step[nSteps - 2].getHomo().low());
      initialAltThresh = getThreshIncreasingGap(initialGap);
      lastThresh = step[nSteps - 2].getChosenThresh();
      initialAltThresh = initialAltThresh > lastThresh / 10 ? 
	initialAltThresh : lastThresh / 10;
      initialGap.puriStep(step[nSteps - 2].getPoly());
    }    
    if (initialGap.empty())
      return;
#if 0
    /* Already converged? */
    if (1.0 - initialGap.upp() < tolEigenvalError && 
	initialGap.low() < tolEigenvalError) {
      stepsLeft = 0;
      thresh = 0;
      return;
    }
#endif
    Treal thetaPerStep = 0.0; /* Tolerated subspace error per step. */
    Treal altThresh = 0; /* Maximum threshold that guarantees increased 
			  * gap.
			  */
    Treal currThresh = 0; /* Chosen threshold. */
    int stepsLeftPrev = -2;
    
    while (stepsLeft > stepsLeftPrev) {
      gap = initialGap;
      altThresh = initialAltThresh;
      currThresh = thetaPerStep * gap.length() / (1 + thetaPerStep);
      currThresh = currThresh < altThresh ? currThresh : altThresh;
      gap.decrease(currThresh);
      lastThresh = currThresh; 
      stepsLeftPrev = stepsLeft;
      stepsLeft = 0;
      while (!gap.empty() && 
	     (1.0 - gap.upp() > tolEigenvalError || 
	      gap.low() > tolEigenvalError)) {
	altThresh = getThreshIncreasingGap(gap); 
	altThresh = altThresh > lastThresh / 10 ? altThresh : lastThresh / 10;
	
	gap.puriStep(gap.upp() + gap.low() < 1);
	
	currThresh = thetaPerStep * gap.length() / (1 + thetaPerStep);
	currThresh = currThresh < altThresh ? currThresh : altThresh;
	gap.decrease(currThresh);
	lastThresh = currThresh; 
	++stepsLeft;
      }
      thetaPerStep = tolError / (stepsLeft + 1);
    }
    
    /* Compute optimal threshold for convergence of subspace. */
    thetaPerStep = tolError / (stepsLeftPrev + 1);
    thresh = thetaPerStep * initialGap.length() / (1 + thetaPerStep);
    return;
  }


  /* FIXME !! */
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    Treal PuriInfo<Treal, Tvector, TdebugPolicy>::getOptimalThresh() const {
    int stepsLeft = -1;
    Treal thresh = 0.0;
    estimateStepsLeft(stepsLeft, thresh);
    if (nSteps > 0)
      step[nSteps - 1].setEstimatedStepsLeft(stepsLeft);
    if (stepsLeft < 0)
      thresh = tolSubspaceError / 100; /* No reason */
    if (nSteps > 1) {
      Interval<Treal> middleGap = 
	Interval<Treal>(step[nSteps - 2].getLumo().upp(),
			step[nSteps - 2].getHomo().low());
      if (!middleGap.empty()) {
	/* Get thresh that definitely gives increasing gap. */
	Treal altThresh = getThreshIncreasingGap(middleGap);
	/* Allow thresh to decrease only to a ten times smaller threshold */
	Treal lastThresh = step[nSteps - 2].getChosenThresh();
	altThresh = altThresh > lastThresh / 10 ? altThresh : lastThresh / 10;
	thresh = thresh < altThresh ? thresh : altThresh;
#if 0
	Interval<Treal> xmX2Eucl = step[nSteps - 2].getXmX2EuclNorm();
	Treal tmp = 1 - xmX2Eucl.upp() * 4.0;
	if (tmp > 0) {
	  Treal altThresh = 0.25 * (1 - template_blas_sqrt(tmp)) / 2; 
	  thresh = thresh < altThresh ? thresh : altThresh;
	}
#endif
      }
    }
    //std::cout<<"nsteps left, thresh: "<<stepsLeft<<" , "<<thresh<<std::endl;
    ASSERTALWAYS(thresh >= 0.0);
    return thresh;
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    Treal PuriInfo<Treal, Tvector, TdebugPolicy>::
    getThreshIncreasingGap(Interval<Treal> const & middleGap) const {
    Treal x = middleGap.midPoint();
    Treal delta = middleGap.upp() - x;
    Treal thresh;
#if 0
    thresh = delta * template_blas_fabs(2 * x - 1) / 10;
#else
    /* After a purification step, truncation is done. 
     * The threshold t (measured by the Euclidean norm of the error matrix)
     * has to satisfy t <= | 1 - 2x | * d / 2 where x and d are the midpoint 
     * and the length of the HOMO-LUMO gap respectively. In this way the 
     * gap is guaranteed to increase to the next step. Note that x = 0.5 and
     * d = 0 gives zero threshold. x -> 0.5 as the process converges. 
     */
    if (delta > 0.4)
      /* This choice gives much better estimation of the remaining 
       * number of iterations.
       */
      /* Close to convergence we chose quadratical dependence on the 
       * midpoint since we expect the midpoint to go to 0.5 quadratically
       * at convergence.
       */
      thresh = delta * template_blas_fabs(2 * x - 1) * template_blas_fabs(2 * x - 1);
    else
      thresh = delta * template_blas_fabs(2 * x - 1) / 2;
    //  thresh = thresh > 1e-7 ? thresh : 1e-7;
#endif
    /**************************************************************
       ELIAS NOTE 2010-11-18: I got assertion failure when using gcc
       in Fedora 14 [g++ (GCC) 4.5.1 20100924 (Red Hat 4.5.1-4)] and
       it turned out to be because the fabs call returned zero while
       thresh was something like 1e-30 for single precision. Therefore
       I added std::numeric_limits<Treal>::epsilon() in the assertion,
       which seemed to solve the problem.
    ***************************************************************/
    ASSERTALWAYS(thresh <= delta * template_blas_fabs(2 * x - 1) + std::numeric_limits<Treal>::epsilon());
    ASSERTALWAYS(thresh >= 0.0);
    return thresh;
  }
  
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    bool PuriInfo<Treal, Tvector, TdebugPolicy>::
    ShouldComputeXmX2EuclNormAccurately(Treal & howAccurate) const {
    
    if (nSteps == 0 || nSteps == 1) {
      howAccurate = 0;
      return false;
    }
    Treal ep = 0.207106781186547; 
    /* = (sqrt(2) - 1) / 2 
     * This value is obtained by noting that x = 1 / sqrt(2) -> x^2 = 1 / 2.
     * If an interval ]1 - 1 / sqrt(2), 1 / sqrt(2)[ is empty from eigenvalues,
     * and if the occ. count is correct, then the occupation 
     * count is guaranteed to be correct after one step as well.
     * This transforms to the value (sqrt(2) - 1) / 2 for the 
     * ||X - X^2||_2 norm via theorem 1.
     */
    Interval<Treal> homoCopy = step[nSteps - 2].getHomo();
    Interval<Treal> lumoCopy = step[nSteps - 2].getLumo();
    homoCopy.puriStep(step[nSteps - 2].getPoly());
    lumoCopy.puriStep(step[nSteps - 2].getPoly());
    /* Note: we changed this from getActualThresh() to
       getChosenThresh() to avoid ridiculously small values for
       howAccurate, as happened earlier for cases when no truncation
       occurred, i.e. when the matrices were small. */
    howAccurate = step[nSteps - 1].getChosenThresh() / 100;
    Treal highestPossibleAccuracy = 10.0 * step[nSteps - 1].getEigAccLoss();
    ASSERTALWAYS(howAccurate >= 0);
    ASSERTALWAYS(highestPossibleAccuracy >= 0);
    howAccurate = howAccurate > highestPossibleAccuracy ?
      howAccurate : highestPossibleAccuracy;
    
    if (homoCopy.length() > 0.2 || lumoCopy.length() > 0.2) {
      /* Base decision on n0 and n1 and XmX2Norm from previous step */
      if (step[nSteps - 2].getN0() / (n - nocc) > 0.5 &&
	  step[nSteps - 2].getN1() / nocc > 0.5 &&
	  step[nSteps - 2].getXmX2EuclNorm().midPoint() < ep)
	return true;
      else
	return false;
    }
    else {
      /* Decision can probably be made from homo and lumo estimates */
      bool computeHomo = true; /* Do we want to try to compute homo */
      bool computeLumo = true; /* Do we want to try to compute lumo */
      if (homoIsAccuratelyKnown() || 
	  homoCopy.upp() < 0.5 ||
	  template_blas_fabs(lumoCopy.low() - 0.5) < template_blas_fabs(homoCopy.low() - 0.5))
	computeHomo = false;
      if (lumoIsAccuratelyKnown() ||
	  lumoCopy.low() > 0.5 ||
	  template_blas_fabs(homoCopy.upp() - 0.5) < template_blas_fabs(lumoCopy.upp() - 0.5))
	computeLumo = false;
      return computeHomo || computeLumo; 
    }
  }

  
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::
    improveHomoLumoF() {
    Treal lmax = eigFInterval.upp();
    Treal lmin = eigFInterval.low();
    Interval<Treal> altHomo(step[0].getHomo() * (lmin - lmax) + lmax);
    Interval<Treal> altLumo(step[0].getLumo() * (lmin - lmax) + lmax);
    homoF.intersect(altHomo);
    lumoF.intersect(altLumo);
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    double PuriInfo<Treal, Tvector, TdebugPolicy>::getAccumulatedTimeSquare() const {
    double accTime = 0;
    for (int ind = 0; ind < nSteps; ++ind) 
      accTime += (double)step[ind].getTimeSquare();
    return accTime;
  }
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    double PuriInfo<Treal, Tvector, TdebugPolicy>::getAccumulatedTimeThresh() const {
    double accTime = 0;
    for (int ind = 0; ind < nSteps; ++ind) 
      accTime += (double)step[ind].getTimeThresh();
    return accTime;    
  }
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    double PuriInfo<Treal, Tvector, TdebugPolicy>::getAccumulatedTimeXmX2Norm() const {
    double accTime = 0;
    for (int ind = 0; ind < nSteps; ++ind) 
      accTime += (double)step[ind].getTimeXmX2Norm();
    return accTime;    
  }
  template<typename Treal, typename Tvector, typename TdebugPolicy>
    double PuriInfo<Treal, Tvector, TdebugPolicy>::getAccumulatedTimeTotal() const {
    double accTime = 0;
    for (int ind = 0; ind < nSteps; ++ind) 
      accTime += (double)step[ind].getTimeTotal();
    return accTime;    
  }


  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::
    mTimings(std::ostream & s) const {
    s<<"puriTime = [";
    for (int ind = 0; ind < nSteps; ++ind) {
      s<<
	step[ind].getTimeSquare()<<"  "<<
	step[ind].getTimeThresh()<<"  "<<
	step[ind].getTimeXmX2Norm()<<"  "<<
	step[ind].getTimeTotal()<<"  "<<
	step[ind].getTimeXX2Write()<<"  "<<
	step[ind].getTimeXX2Read()<<"  "<<
	std::endl;
    }
    s<<"];\n"<<"figure; bar(puriTime(:,1:3),'stacked')"<<std::endl<<
      "legend('Matrix Square', 'Truncation', '||X-X^2||'),"<<
      " xlabel('Iteration'), ylabel('Time (seconds)')"<<std::endl;
    s<<"figure; plot(puriTime(:,4),'-x')"<<std::endl;
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::
    mInfo(std::ostream & s) const {
    s<<"% PURIFICATION INFO IN MATLAB/OCTAVE FILE FORMAT"<<std::endl;
    s<<"% Norm for truncation: "<< getNormTypeString(normForTruncation) 
     << std::endl;
    s<<"n = "<<n<<";"<<std::endl;
    s<<"nocc = "<<nocc<<";"<<std::endl;
    s<<"tolSubspaceError = "<<tolSubspaceError<<";"<<std::endl;
    s<<"tolEigenvalError = "<<tolEigenvalError<<";"<<std::endl;
    s<<"correct_occupation_was_forced_flag = "
     <<correct_occupation_was_forced_flag<<";"<<std::endl;
    s<<"% Each row of the following matrix contains purification info \n"
     <<"% for one step.\n"
     <<"% The columns are arranged as: \n"
     <<"% ind, n0, n1, trace(X), trace(X*X), poly, actualThresh, delta, "
     <<"correctOcc, XmX2low, XmX2upp, HOMOmid, LUMOmid, "
     <<"nnz(X), nnz(X*X), chosenThresh, estimatedStepsLeft "
     <<std::endl;
    s<<"puriMat = [";
    for (int ind = 0; ind < nSteps; ++ind) {
      s<<
	ind+1<<"  "<<
	step[ind].getN0()<<"  "<<
	step[ind].getN1()<<"  "<<
	step[ind].getTraceX()<<"  "<<
	step[ind].getTraceX2()<<"  "<<
	step[ind].getPoly()<<"  "<<
	step[ind].getActualThresh()<<"  "<<
	step[ind].getDelta()<<"  "<<
	step[ind].getCorrectOccupation()<<"  "<<
	step[ind].getXmX2EuclNorm().low()<<"  "<<
	step[ind].getXmX2EuclNorm().upp()<<"  "<<
	step[ind].getHomo().midPoint()<<"  "<<
	step[ind].getLumo().midPoint()<<"  "<<
	step[ind].getNnzX()<<"  "<<
	step[ind].getNnzX2()<<"  "<<
	step[ind].getChosenThresh()<<"  "<<
	step[ind].getEstimatedStepsLeft()<<"  "<<
	std::endl;
    }
    s<<"];\n";
    s<<"if(1)\n";
    s<<"ind = puriMat(:,1);\n";
    s<<"figure; \n"
     <<"plot(ind,puriMat(:,12),'x-b',ind,puriMat(:,13),'o-r')\n"
     <<"legend('HOMO','LUMO'),xlabel('Iteration')\n"
     <<"axis([0 max(ind) 0 1])\n";
    s<<"figure; \n"
     <<"plot(ind,100*puriMat(:,2)/(n-nocc),'o-r',...\n"
     <<"ind, 100*puriMat(:,3)/nocc,'x-b')\n"
     <<"legend('N0','N1'),xlabel('Iteration'), ylabel('Percentage')\n"
     <<"axis([0 max(ind) 0 100])\n";
    s<<"figure; \n"
     <<"subplot(211)\n"
     <<"plot(ind,100*puriMat(:,14)/(n*n),'o-r',...\n"
     <<"ind, 100*puriMat(:,15)/(n*n),'x-b')\n"
     <<"legend('nnz(X)','nnz(X^2)'),xlabel('Iteration') \n"
     <<"ylabel('Percentage')\n"
     <<"axis([0 max(ind) 0 100])\n";
    s<<"subplot(212)\n"
     <<"semilogy(ind,puriMat(:,16),'x-r',ind,puriMat(:,7),'o-b')\n"
     <<"xlabel('Iteration'), ylabel('Threshold')\n"
     <<"legend('Chosen threshold', 'Actual threshold')\n"
     <<"axis([0 max(ind) min(puriMat(:,7))/10 max(puriMat(:,16))*10])\n";
    s<<"figure; \n"
     <<"plot(ind,ind(end:-1:1)-1,ind,puriMat(:,17),'x-r')\n"
     <<"legend('Steps left', 'Estimated steps left')\n"
     <<"xlabel('Iteration')\n";

    s<<"end\n";
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::
    mMemUsage(std::ostream & s) const {
    s<<"puriMemUsage = [";
    for (int ind = 0; ind < nSteps; ++ind) {
      MemUsage::Values memUsageBeforeTrunc = step[ind].getMemUsageBeforeTrunc();
      MemUsage::Values memUsageInXmX2Diff  = step[ind].getMemUsageInXmX2Diff();
      s<<
	memUsageBeforeTrunc.res <<"  "<<
	memUsageBeforeTrunc.virt<<"  "<<
	memUsageBeforeTrunc.peak<<"  "<<
	memUsageInXmX2Diff.res <<"  "<<
	memUsageInXmX2Diff.virt<<"  "<<
	memUsageInXmX2Diff.peak<<"  "<<
	std::endl;
    }
    s<<"];\n";
    s<<"figure; \n"
     <<"plot(puriMemUsage(:,1),'x-b')\n"
     << "hold on\n"
     <<"plot(puriMemUsage(:,2),'o-r')\n"
     <<"plot(puriMemUsage(:,3),'^-g')\n"
     <<"legend('Resident','Virtual','Peak'),xlabel('Iteration'),ylabel('Mem usage before trunc [GB]')\n"
     <<"%force y axis to start at 0\n"
     <<"axissaved = axis; axissaved(3) = 0; axis(axissaved);\n"
     << std::endl;
    s<<"figure; \n"
     <<"plot(puriMemUsage(:,4),'x-b')\n"
     << "hold on\n"
     <<"plot(puriMemUsage(:,5),'o-r')\n"
     <<"plot(puriMemUsage(:,6),'^-g')\n"
     <<"legend('Resident','Virtual','Peak'),xlabel('Iteration'),ylabel('Mem usage in XmX2Diff [GB]')\n"
     <<"%force y axis to start at 0\n"
     <<"axissaved = axis; axissaved(3) = 0; axis(axissaved);\n"
     << std::endl;
  }

  template<typename Treal, typename Tvector, typename TdebugPolicy>
    void PuriInfo<Treal, Tvector, TdebugPolicy>::
    getHOMOandLUMOeigVecs(Tvector & eigVecLUMO, 
			  Tvector & eigVecHOMO) const {
    bool haveHOMO = 0;
    bool haveLUMO = 0;
    for (int ind = 0; ind < nSteps; ++ind) {
      if (!haveHOMO && step[ind].getHomoWasComputed()) {
	eigVecHOMO = *step[ind].getEigVecPtr();
	haveHOMO = 1;
      }
      if (!haveLUMO && step[ind].getLumoWasComputed()) {
	eigVecLUMO = *step[ind].getEigVecPtr();
	haveLUMO = 1;
      }
    }
  }

  
} /* end namespace mat */
#undef __ID__
#endif
