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

#include "TestMatrix.h"
#include "Purification_scaled.h"
int main() {
  typedef double real;
  int n = 10;
  int nocc = 5;
  real * elements = new real[n];
  for (int i = 0; i < n; ++i) 
    elements[i] = i;

  TestMatrix<real> F_and_D(n, elements);
  mat::Interval<real> eigFInt( F_and_D.min(), F_and_D.max() );
  mat::Interval<real> hoF( elements[nocc-1]-0.1, elements[nocc-1]+0.01 );
  mat::Interval<real> luF( elements[nocc]-0.01,   elements[nocc]+0.1 );
  real const toleratedEigenvalError = 1e-6;
  real const toleratedSubspaceError = 1e-6;
  int const max_steps = 100;
  mat::normType normForTruncation = mat::euclNorm; 
  bool use_scaling = 1;
  pur::Purification_scaled<TestMatrix<real> > myPuri( F_and_D,  
						      eigFInt,
						      hoF, luF,
						      toleratedEigenvalError,
						      toleratedSubspaceError,
						      max_steps,
						      normForTruncation,
						      use_scaling );
  myPuri.purify();

  F_and_D.get_diag(elements);
  for (int i = 0; i < n; ++i) 
    std::cout << elements[i] << "  ";
  std::cout << std::endl;

  mat::Interval<real> hoF_new;
  mat::Interval<real> luF_new;
  myPuri.get_homo_lumo_intervals(hoF_new, luF_new);
  std::cout << "NEW HOMO: " << hoF_new << std::endl;
  std::cout << "NEW LUMO: " << luF_new << std::endl;

  std::ofstream ff("puriInfo.m");
  myPuri.mInfo(ff);
  ff.close();
  //  myPuri.mInfo( std::cout );

  std::ofstream ff2("puriTime.m");
  myPuri.mTime(ff2);
  ff2.close();

}
