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

#include "matInclude.h"
#include "mat_gblas.h"
#include "matrix_proxy.h"
#include "Interval.h"
template<typename Treal>
class TestMatrix {
 public:
  typedef Treal real;
  TestMatrix() :elements(0),n(0) {}  
  ~TestMatrix() {delete[] elements;}
  TestMatrix(TestMatrix<real> const & other);
  TestMatrix<real>& operator=(TestMatrix<real> const & other);
  TestMatrix(int const n, real const * const elem);
  void get_diag(real * elem) const;
  real thresh(real chosen_thresh, mat::normType normTruncation);
  TestMatrix<real>& operator*=(real const alpha);
  void add_identity(real const alpha);
  TestMatrix<real>& operator=(mat::XYZ<real, TestMatrix<real>, TestMatrix<real> > const & sm2);
  TestMatrix<real>& operator+=(mat::XY<real, TestMatrix<real> > const & sm);
  real trace() const;
  static mat::Interval<real> diffIfSmall( TestMatrix<real> const & A, 
					  TestMatrix<real> const & B, 
					  mat::normType const norm, 
					  real const reqAcc, 
					  real const maxAbsVal );
  static real mixed_diff( TestMatrix<real> const & A, 
			  TestMatrix<real> const & B, 
			  real const reqAcc );

  void transfer(TestMatrix<real> & dest);

  
  real min() const {
    real min_val = elements[0];
    for (int ind = 1; ind < n; ++ind) {
      min_val = min_val < elements[ind] ? min_val : elements[ind];
    }
    return min_val;
  }
  real max() const {
    real max_val = elements[0];
    for (int ind = 1; ind < n; ++ind) {
      max_val = max_val > elements[ind] ? max_val : elements[ind];
    }
    return max_val;
  }

  size_t nnz() const {return n;}
  
  int get_nrows() const {return n;}
 private:
  real * elements; // Diagonal matrix, length of elements vector: n
  int n;
};

template<typename Treal>
TestMatrix<Treal>::TestMatrix(TestMatrix<Treal> const & other) 
:n(other.n)
{
  elements = new real[n];
  for (int i = 0; i < n; ++i) 
    elements[i] = other.elements[i];  
}

template<typename Treal>
TestMatrix<Treal>& TestMatrix<Treal>::operator=(TestMatrix<Treal> const & other) {
  delete[] elements;
  n = other.n;
  elements = new real[n];
  for (int i = 0; i < n; ++i) 
    elements[i] = other.elements[i];    
}

template<typename Treal>
TestMatrix<Treal>::TestMatrix(int const n, real const * const elem)
:n(n) {
  elements = new real[n];
  for (int i = 0; i < n; ++i) 
    elements[i] = elem[i];
} 

template<typename Treal>
void TestMatrix<Treal>::get_diag(real * elem) const {
  for (int i = 0; i < n; ++i) 
    elem[i] = elements[i];  
}


template<typename Treal>
Treal TestMatrix<Treal>::thresh(real chosen_thresh, 
			      mat::normType normTruncation) {
  for (int ind = 0; ind < n; ind++)
    elements[ind] += chosen_thresh;
  return chosen_thresh;
}

template<typename Treal>
TestMatrix<Treal>& TestMatrix<Treal>::operator*=(real const alpha) {
  for (int ind = 0; ind < n; ind++)
    elements[ind] *= alpha;
  return *this;
}

template<typename Treal>
void TestMatrix<Treal>::add_identity(real const alpha) {
  for (int ind = 0; ind < n; ind++)
    elements[ind] += alpha;
}

template<typename Treal>
TestMatrix<Treal>& TestMatrix<Treal>::operator=(mat::XYZ<real, TestMatrix<Treal>, TestMatrix<Treal> > const & sm2) {
  assert(this != &sm2.B);
  assert(this != &sm2.C);    
  delete[] elements;
  n = sm2.B.n;
  elements = new real[n];
  for (int ind = 0; ind < n; ind++)
    elements[ind] = sm2.A * sm2.B.elements[ind] * sm2.C.elements[ind];
  return *this;
}
  
template<typename Treal>
TestMatrix<Treal>& TestMatrix<Treal>::
operator+=(mat::XY<real, TestMatrix<Treal> > const & sm) {
  for (int ind = 0; ind < n; ind++)
    elements[ind] += sm.A * sm.B.elements[ind];
  return *this;  
}

template<typename Treal>
Treal TestMatrix<Treal>::trace() const {
  real tr = 0;
  for (int ind = 0; ind < n; ind++)
    tr += elements[ind];
  return tr;
}

template<typename Treal>
mat::Interval<Treal> TestMatrix<Treal>::
diffIfSmall( TestMatrix<Treal> const & A, 
	     TestMatrix<Treal> const & B, 
	     mat::normType const norm, 
	     real const reqAcc, 
	     real const maxAbsVal ) {
  real diff = 0;
  for (int ind = 0; ind < A.n; ind++) {
    real tmp_val = fabs( A.elements[ind] - B.elements[ind] ); 
    diff = diff > tmp_val ? diff : tmp_val;
  }
  return mat::Interval<Treal>(diff,diff);
}

template<typename Treal>
Treal TestMatrix<Treal>::mixed_diff( TestMatrix<Treal> const & A, 
				     TestMatrix<Treal> const & B, 
				     real const reqAcc ) {
  real diff = 0;
  for (int ind = 0; ind < A.n; ind++) {
    real tmp_val = fabs( A.elements[ind] - B.elements[ind] ); 
    diff = diff > tmp_val ? diff : tmp_val;
  }
  return diff;
}


template<typename Treal>
void TestMatrix<Treal>::
transfer(TestMatrix<real> & dest) {
  delete[] dest.elements;
  dest.n = n;
  dest.elements = elements;
  elements = 0;
  n = 0;
}
