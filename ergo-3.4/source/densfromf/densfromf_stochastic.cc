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

/** @file densfromf_stochastic.cc

    \brief Routine get_dens_from_fock_stochastic() for getting density matrix from a given Fock matrix using stochastic orbitals.

    @author: Elias Rudberg <em>responsible</em>. 
*/
#include "densfromf_stochastic.h"
#include "utilities.h"
#include "AllocatorManager.h"
#include <math.h>


static ergo_real get_scalar_product(const std::vector<ergo_real> & a, const std::vector<ergo_real> & b) {
  if(a.size() != b.size())
    throw "Error in get_scalar_product: (a.size() != b.size()).";
  ergo_real sum = 0;
  for(size_t i = 0; i < a.size(); i++)
    sum += a[i] * b[i];
  return sum;
}


static int verify_that_vector_is_normalized(const std::vector<ergo_real> & v) {
  ergo_real sqSum = 0;
  for(size_t i = 0; i < v.size(); i++)
    sqSum += v[i] * v[i];
  ergo_real diff_from_one = std::fabs(sqSum-1);
  if(diff_from_one > 1e-11) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "Error: verify_that_vector_is_normalized failure; diff_from_one = %8.4g", (double)diff_from_one);
    return -1;
  }
  return 0;
}


int get_dens_from_fock_stochastic(int n,
				  int noOfOccupiedOrbs,
				  symmMatrix & resultDens,
				  ergo_real factor,
				  symmMatrix const & Finput,
				  triangMatrix const & invCholFactor,
				  mat::SizesAndBlocks const & matrixSizesAndBlocks,
				  const std::vector< std::vector<ergo_real> > stochastic_orbitals) {
  Util::TimeMeter timeMeterTot;
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_stochastic() start!");

  symmMatrix F(Finput);
  F.readFromFile();
  std::string allocStatsStr3 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "After F.readFromFile(): %s", allocStatsStr3.c_str());
 
  {
    triangMatrix invCholFactor_tmp(invCholFactor);
    invCholFactor_tmp.readFromFile();
    output_current_memory_usage(LOG_AREA_DENSFROMF, "In get_dens_from_fock_sparse, before F = tr(Z) * F * Z");
    Util::TimeMeter timeMeterFortTransf;
    F = transpose(invCholFactor_tmp) * F * invCholFactor_tmp;
    timeMeterFortTransf.print(LOG_AREA_DENSFROMF, "F = transpose(invCholFactor) * F * invCholFactor");
    output_current_memory_usage(LOG_AREA_DENSFROMF, "In get_dens_from_fock_sparse,  after F = tr(Z) * F * Z");
  } // invCholFactor_tmp goes out of scope here

  // Now F contains F_ort. 

  // Diagonalize F.
  std::vector<ergo_real> F_full(n*n);
  {
    // Create full matrix version of F
    normalMatrix* tmpMat;
    tmpMat = new normalMatrix(F);
    tmpMat->fullMatrix(F_full);
    delete tmpMat;
  }

  std::vector<ergo_real> eigVals(n);
  std::vector< std::vector<ergo_real> > eigVecs(n);
  {
    std::vector<ergo_real> A(n*n);
    for(int k = 0; k < n*n; k++)
      A[k] = F_full[k];
    int info = -1;
    int lwork = 3*n-1;
    std::vector<ergo_real> work(lwork);
    Util::TimeMeter timeMeter_syev;
    mat::syev("V", 
	      "U", 
	      &n, 
	      &A[0], 
	      &n, 
	      &eigVals[0], 
	      &work[0], 
	      &lwork, 
	      &info);
    timeMeter_syev.print(LOG_AREA_DENSFROMF, "syev diagonalization");
    if(info != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "Error in sygv");
      return -1;
    }
    for(int i = 0; i < n; i++) {
      eigVecs[i].resize(n);
      for(int k = 0; k < n; k++)
	eigVecs[i][k] = A[i*n+k];
    }
  }
  // Test that computed eigenpairs seem OK.
  for(int i = 0; i < n; i++) {
    std::vector<ergo_real> Fv(n);
    for(int k = 0; k < n; k++) {
      ergo_real sum = 0;
      for(int m = 0; m < n; m++)
	sum += F_full[k*n+m]*eigVecs[i][m];
      Fv[k] = sum;
    }
    ergo_real maxabsdiff = 0;
    for(int k = 0; k < n; k++) {
      ergo_real absdiff = std::fabs(Fv[k] - eigVals[i] * eigVecs[i][k]);
      if(absdiff > maxabsdiff)
	maxabsdiff = absdiff;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxabsdiff for eigenpair %6d: %9.4g", i, (double)maxabsdiff);
    if(maxabsdiff > 1e-8) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "Error: too large diff when checking result after sygv");
      return -1;
    }
  }
  // Verify that list of eigenpairs is sorted.
  for(int i = 0; i < n-1; i++) {
    ergo_real diff = eigVals[i+1] - eigVals[i];
    if(diff < 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "Error: eigenvalues not sorted in ascending order.");
      return -1;
    }
  }
  // Verify that all eigenvectors are normalized
  for(int i = 0; i < n; i++) {
    if(verify_that_vector_is_normalized(eigVecs[i]) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "Error: eigenvectors not normalized.");
      return -1;
    }
  }

  int no_of_stochastic_orbitals = stochastic_orbitals.size();
  std::vector< std::vector<ergo_real> > projected_orbitals(no_of_stochastic_orbitals);
  ergo_real stochastic_orbital_scale_factor = ( (ergo_real)n / no_of_stochastic_orbitals);

  // Verify that all stochastic orbital vectors are normalized
  for(int i = 0; i < no_of_stochastic_orbitals; i++) {
    if(verify_that_vector_is_normalized(stochastic_orbitals[i]) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "Error: stochastic orbital vectors not normalized.");
      return -1;
    }
  }

  ergo_real chemical_potential_lo = -100; // FIXME use better limits here
  ergo_real chemical_potential_hi =  100; // FIXME use better limits here
  const ergo_real beta = 100; // FIXME use other value here
  while(1) {
    ergo_real chemical_potential = ( chemical_potential_lo + chemical_potential_hi ) / 2;
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "=============================================");
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "chemical_potential = %15.8f", chemical_potential);
    // Set up vector of occupation numbers.
    std::vector<ergo_real> occupation_numbers(n);
    for(int i = 0; i < n; i++) {
      occupation_numbers[i] = 0.5 * erfc(-beta*(chemical_potential-eigVals[i]));
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "occupation_numbers[%5d] = %14.9f", i, (double)occupation_numbers[i]);
    }
    // Now compute the projection of each stochastic vector onto the occupied subspace.
    for(int i = 0; i < no_of_stochastic_orbitals; i++) {
      // Now take care of stochastic vector i
      projected_orbitals[i].resize(n);
      for(int k = 0; k < n; k++)
	projected_orbitals[i][k] = 0;
      for(int k = 0; k < n; k++) {
	ergo_real scalarProd = get_scalar_product(stochastic_orbitals[i], eigVecs[k]);
	for(int m = 0; m < n; m++)
	  projected_orbitals[i][m] += occupation_numbers[k] * scalarProd * eigVecs[k][m];
      }
    } // end for i
    // Now compute trace of "density matrix" from stochastic orbitals, one stochastic orbital at a time.
    double occupationCount = 0;
    for(int i = 0; i < no_of_stochastic_orbitals; i++) {
      // Now take care of stochastic vector i
      ergo_real sqSum = 0;
      for(int k = 0; k < n; k++)
	sqSum += projected_orbitals[i][k] * projected_orbitals[i][k];
      occupationCount += sqSum * stochastic_orbital_scale_factor;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "occupationCount = %15.8f, noOfOccupiedOrbs = %d", occupationCount, noOfOccupiedOrbs);
    if(occupationCount > noOfOccupiedOrbs) {
      // Too high occupation. Then we know the chemical potential should be smaller than the value used now, so we can move down the upper bound.
      chemical_potential_hi = chemical_potential;
    }
    else {
      chemical_potential_lo = chemical_potential;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "chemical_potential_lo chemical_potential_hi = %15.4f %15.4f", chemical_potential_lo, chemical_potential_hi);
    if(chemical_potential_hi - chemical_potential_lo < 1e-8)
      break;
  } // end while trying different chemical potential values

  // OK, now contruct new density matrix from the projected stochastic orbitals.
  std::vector<ergo_real> D_full(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      D_full[i*n+j] = 0;
  // Add contributions to density matrix for one stochastic orbital at a time.
  for(int i = 0; i < no_of_stochastic_orbitals; i++) {
    for(int k = 0; k < n; k++)
      for(int m = 0; m < n; m++)
	D_full[k*n+m] += stochastic_orbital_scale_factor * projected_orbitals[i][k] * projected_orbitals[i][m];
  }
  // OK, now D_full is the new density matrix computed from the projected stochastic orbitals.

  resultDens.assignFromFull(D_full);  
  {
    triangMatrix invCholFactor_tmp(invCholFactor);
    invCholFactor_tmp.readFromFile();
    resultDens = invCholFactor_tmp * resultDens * transpose(invCholFactor_tmp);
  } // invCholFactor_tmp goes out of scope here
  
  resultDens *= factor;
  
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_stochastic ending OK");
  timeMeterTot.print(LOG_AREA_DENSFROMF, "get_dens_from_fock_stochastic");
  
  return 0;
}

