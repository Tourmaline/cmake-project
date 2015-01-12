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

/** @file purification_old.h Purification methods
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date 2005-2006
 *
 */
#ifndef MAT_PURIFICATION_OLD
#define MAT_PURIFICATION_OLD
#include <iomanip>
#include <math.h>
namespace mat{


  template<class Treal, class MAT>
    void tc2_extra(MAT& D, const int nshells, 
		   const int nsteps, const Treal frob_trunc = 0) {
    MAT tmp(D);
    for (int step = 0; step < nsteps; step++) {
      Treal tracediff =  tmp.trace() - nshells;
      if (tracediff > 0) {
	tmp = (Treal)1.0 * D * D;
	D = tmp;
	D = (Treal)-1.0 * tmp * tmp + (Treal)2.0 * D;
      }
      else {
	tmp = D;
	tmp = (Treal)-1.0 * D * D + (Treal)2.0 * tmp;
	D = (Treal)1.0 * tmp * tmp;	
      }      
      D.frob_thresh(frob_trunc);
    }      
  } 

  template<class Treal, class MAT>
    void tc2_extra_auto(MAT& D, const int nshells, 
			int& nsteps, const Treal frob_trunc,
			const int maxiter = 50) {
    MAT tmp(D);
    Treal froberror;
    Treal tracediff =  tmp.trace() - nshells;
    nsteps = 0;
    do { 
      if (tracediff > 0) {
	tmp = (Treal)1.0 * D * D;
	D = tmp;
	D = (Treal)-1.0 * tmp * tmp + (Treal)2.0 * D;
      }
      else {
	tmp = D;
	tmp = (Treal)-1.0 * D * D + (Treal)2.0 * tmp;
	D = (Treal)1.0 * tmp * tmp;	
      }
      froberror = MAT::frob_diff(D, tmp);
      D.frob_thresh(frob_trunc);
      tracediff = D.trace() - nshells;
#if 0
      std::cout<<"Iteration:"<<nsteps<<"  Froberror: "
	       <<std::setprecision(8)<<froberror<<
	"  Tracediff: "<<tracediff<<std::endl;
#endif
      nsteps++;
      if (nsteps > maxiter)
	throw AcceptableMaxIter("Extra Auto Purification reached maxiter" 
				" without convergence", maxiter);
    } while (froberror > frob_trunc);
  } 

  
  
  template<class Treal, class MAT>
    void tc2(MAT& D, int& iterations,
	     const MAT& H, const int nshells, 
	     const Treal trace_error_limit = 1e-2,
	     const Treal frob_trunc = 0,
	     int* polys= NULL,
	     const int maxiter = 100) {
    MAT tmp(H);
    D = H;
    iterations = 0;
    Treal tracediff =  tmp.trace() - nshells;
    Treal tracediff_old = 2.0 * trace_error_limit;
    
    while (template_blas_fabs(tracediff) > trace_error_limit && 
	   template_blas_fabs(tracediff_old) > trace_error_limit) {
      //      std::cout<<"Iteration:"<<iterations
      //	       <<"   Tracediff ="<<tracediff<<std::endl;
      if (tracediff > 0) {
	D = (Treal)1.0 * tmp * tmp;
	if (polys)
	  polys[iterations] = 0;
      }
      else {
	D = (Treal)-1.0 * tmp * tmp + (Treal)2.0 * D;
	if (polys)
	  polys[iterations] = 1;
      }
      D.frob_thresh(frob_trunc);
      tmp = D;
      tracediff_old = tracediff;
      tracediff =  D.trace() - nshells;
      iterations++;
      if (iterations > maxiter)
	throw AcceptableMaxIter("Purification reached maxiter without convergence", maxiter);
    }
  }

#if 1
  template<class Treal, class MAT>
    void tc2_auto(MAT& D, int& iterations,
		  const MAT& H, const int nocc, 
		  const Treal frob_trunc = 0,
		  int* polys= NULL,
		  const int maxiter = 100) {
    assert(frob_trunc >= 0);
    assert(nocc >= 0);
    assert(maxiter >= 0);
    MAT tmp(H);
    D = H;
    iterations = 0;
    Treal tracediff = 2;
    Treal newtracediff = tmp.trace() - nocc;
    Treal froberror = 1;
    Treal tracechange = 0; /* Initial value has no relevance */
    /* Need check for too loose truncation? */
    while (2 * froberror > (3 - template_blas_sqrt((Treal)5)) / 2 - frob_trunc || 
	   template_blas_fabs(tracediff) > 1 || 
	   template_blas_fabs(tracechange) > (1 - 2 * froberror) * (1 - template_blas_fabs(tracediff))) {
      //      std::cout<<"Iteration:"<<iterations
      //	       <<"   Tracediff ="<<tracediff<<std::endl;
      tracediff = newtracediff;
      if (tracediff > 0) {
	D = (Treal)1.0 * tmp * tmp;
	if (polys)
	  polys[iterations] = 0;
      }
      else {
	D = (Treal)-1.0 * tmp * tmp + (Treal)2.0 * D;
	if (polys)
	  polys[iterations] = 1;
      }
      D.frob_thresh(frob_trunc);
      newtracediff = D.trace() - nocc;
      tracechange = newtracediff - tracediff;
      froberror = MAT::frob_diff(D, tmp);
      tmp = D;
#if 0
      std::cout<<"Iteration:"<<iterations<<"  epsilon: "
	       <<std::setprecision(8)<<std::setw(12)<<froberror<<
	"  delta: "<<std::setw(12)<<tracediff<<
	"  tracechange: "<<std::setw(12)<<tracechange<<std::endl;
#endif
      iterations++;
      if (iterations > maxiter)
	throw AcceptableMaxIter("tc2_auto::Purification reached maxiter"
				" without convergence", maxiter);
    }
    
    /* Take one step to make sure the eigenvalues stays in */
    /* { [ 0 , 2 * epsilon [ , ] 1 - 2 * epsilon , 1] }    */
    if (tracediff > 0) {
      D = (Treal)-1.0 * tmp * tmp + (Treal)2.0 * D;
      if (polys)
	polys[iterations] = 1;
    }
    else {
      D = (Treal)1.0 * tmp * tmp;
      if (polys)
	polys[iterations] = 0;
    }
    iterations++;

    /* Use second order convergence polynomials to converge completely */
    D.frob_thresh(frob_trunc);
    tracediff = D.trace() - nocc;
    do { 
      if (tracediff > 0) {
	tmp = (Treal)1.0 * D * D;
	D = tmp;
	D = (Treal)-1.0 * tmp * tmp + (Treal)2.0 * D;
	if (polys) {
	  polys[iterations] = 0;
	  polys[iterations + 1] = 1;
	}
      }
      else {
	tmp = D;
	tmp = (Treal)-1.0 * D * D + (Treal)2.0 * tmp;
	D = (Treal)1.0 * tmp * tmp;
	if (polys) {
	  polys[iterations] = 1;
	  polys[iterations + 1] = 0;
	}
      }
      froberror = MAT::frob_diff(D, tmp);
      D.frob_thresh(frob_trunc);
      tracediff = D.trace() - nocc;
#if 0
      std::cout<<"Iteration:"<<iterations<<"  Froberror: "
	       <<std::setprecision(8)<<froberror<<
	"  Tracediff: "<<tracediff<<std::endl;
#endif
      iterations += 2;
      if (iterations > maxiter)
	throw AcceptableMaxIter("tc2_auto::Purification reached maxiter" 
				" in extra part without convergence", maxiter);
    } while (froberror > frob_trunc);
    

  }

#endif



#if 0
  iterations = 0;
  /* D = (F - lmax * I) / (lmin - lmax)                                  */
  D = F;
  D.add_identity(-lmax);       
    D *= (1.0 / (lmin - lmax));
    /*****/ clock_t start;
    /*****/ float syrktime = 0;
    /*****/ float threshtime = 0;
    MAT tmp;     
    Treal tracediff;
    for (int steps = 0; steps < nextra_steps; steps++) {
      tmp = D;
      tracediff = tmp.trace() - nshells;
      while (template_blas_fabs(tracediff) > trace_error_limit) {
	iterations++;
	/*****/ start = clock();
	if (tracediff > 0) 
	  MAT::syrk('U', false, 1.0, tmp, 0.0, D);
	else 
	  MAT::syrk('U', false, -1.0, tmp, 2.0, D);
        /*****/	syrktime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
	D.sym_to_nosym();
	/*****/ start = clock();
	D.frob_thresh(frob_trunc);
        /*****/	threshtime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
	tmp = D;
	tracediff = tmp.trace() - nshells;
      } /* end while */
      
      /* extra steps */
      iterations += 2;
      if (tracediff > 0) {
	/*****/ start = clock();
	MAT::syrk('U', false, 1.0, tmp, 0.0, D);
        /*****/	syrktime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
	D.sym_to_nosym();
	tmp = D;
	/*****/ start = clock();
	MAT::syrk('U', false, -1.0, tmp, 2.0, D);
        /*****/	syrktime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
      }
      else {
	/*****/ start = clock();
	MAT::syrk('U', false, -1.0, tmp, 2.0, D);
        /*****/	syrktime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
	D.sym_to_nosym();
	tmp = D;
	/*****/ start = clock();
	MAT::syrk('U', false, 1.0, tmp, 0.0, D);
        /*****/	syrktime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
      }      
      D.sym_to_nosym();
	/*****/ start = clock();
      D.frob_thresh(frob_trunc);
        /*****/	threshtime += ((float)(clock()-start))/(CLOCKS_PER_SEC);
    } /* end for */
    /** Continue here */
    std::cout<<std::endl<<" total syrktime in tcp ="<<std::setw(5)
             <<syrktime<<std::endl;
    std::cout<<" total threshtime in tcp ="<<std::setw(5)
             <<threshtime<<std::endl;

  }
#endif  
} /* end namespace mat */

#endif
