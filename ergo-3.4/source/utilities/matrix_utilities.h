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

#ifndef MATRIX_UTILITIES_HEADER
#define MATRIX_UTILITIES_HEADER

#include "matrix_typedefs.h"
#include "basisinfo.h"

#if 0
/** prepare_matrix_permutation creates a perm object as required by
    further matrix manipulation routines. */
Perm* prepare_matrix_permutation(const BasisInfoStruct& basisInfo,
				 int sparse_block_size,
				 int factor1, int factor2, int factor3);
#else

mat::SizesAndBlocks prepareMatrixSizesAndBlocks(int n_basis_functions,
						int sparse_block_size,
						int factor1, 
						int factor2, 
						int factor3);

void getMatrixPermutation(const BasisInfoStruct& basisInfo,
			  int sparse_block_size,
			  int factor1, 
			  int factor2, 
			  int factor3,
			  std::vector<int> & permutation);
void getMatrixPermutation(const BasisInfoStruct& basisInfo,
			  int sparse_block_size,
			  int factor1, 
			  int factor2, 
			  int factor3,
			  std::vector<int> & permutation,
			  std::vector<int> & inversePermutation);
void getMatrixPermutationOnlyFactor2(const std::vector<ergo_real> & xcoords,
				     const std::vector<ergo_real> & ycoords,
				     const std::vector<ergo_real> & zcoords,
				     int sparse_block_size_lowest,
				     int first_factor, // this factor may be different from 2, all other factors are always 2.
				     std::vector<int> & permutation,
				     std::vector<int> & inversePermutation);
void getMatrixPermutationOnlyFactor2(const BasisInfoStruct& basisInfo,
				     int sparse_block_size_lowest,
				     int first_factor, // this factor may be different from 2, all other factors are always 2.
				     std::vector<int> & permutation,
				     std::vector<int> & inversePermutation);

#endif
void fill_matrix_with_random_numbers(int n, symmMatrix & M);

void add_random_diag_perturbation(int n, 
				  symmMatrix & M, 
				  ergo_real eps);

bool check_if_matrix_contains_strange_elements(const symmMatrix & M,
                                               std::vector<int> const & inversePermutationHML);

void output_matrix(int n, const ergo_real* matrix, const char* matrixName);

template<class Tmatrix>
ergo_real compute_maxabs_sparse(const Tmatrix & M)
{
  return M.maxAbsValue();

}

template<typename RandomAccessIterator>
struct matrix_utilities_CompareClass {
  RandomAccessIterator first;
  explicit matrix_utilities_CompareClass(RandomAccessIterator firstel)
    : first(firstel){}
  bool operator() (int i,int j) { return (*(first + i) < *(first + j));}
};

template<typename Tmatrix>
void get_all_nonzeros( Tmatrix const & A, 
		       std::vector<int> const & inversePermutation,
		       std::vector<int> & rowind,
		       std::vector<int> & colind,
		       std::vector<ergo_real> & values) {
  rowind.resize(0);
  colind.resize(0);
  values.resize(0);
  size_t nvalues = 0;
  size_t nvalues_tmp = A.nvalues();
  std::vector<int>       rowind_tmp; rowind_tmp.reserve(nvalues_tmp);
  std::vector<int>       colind_tmp; colind_tmp.reserve(nvalues_tmp);
  std::vector<ergo_real> values_tmp; values_tmp.reserve(nvalues_tmp);
  A.get_all_values(rowind_tmp,
		   colind_tmp,
		   values_tmp,
		   inversePermutation,
		   inversePermutation);
  // Count the number of nonzeros
  for(size_t i = 0; i < nvalues_tmp; i++) {
    nvalues += ( values_tmp[i] != 0 );
  }
  rowind.reserve(nvalues);
  colind.reserve(nvalues);
  values.reserve(nvalues);
  // Extract all nonzeros
  for(size_t i = 0; i < nvalues_tmp; i++) {
    if ( values_tmp[i] != 0 ) {
      rowind.push_back( rowind_tmp[i] );
      colind.push_back( colind_tmp[i] );
      values.push_back( values_tmp[i] );	
    }
  }
} // end get_all_nonzeros(...)



template<typename Tmatrix>
void output_distance_vs_magnitude( BasisInfoStruct const & basisInfo, /**< Info about basis set. */
				   Tmatrix const & A, /**< The matrix. */
				   std::vector<int> const & inversePermutation, /**< Permutation to be used when accessing matrix elements. */
				   std::string name, /**< File name. */
				   int resolution_r, /**< Resolution in r-direction, r is the distance between atoms. */
				   int resolution_m  /**< Resolution m-direction, m is the magnitude of matrix elements. */
				   ) {
  std::string m_name = name + ".m";
  std::ofstream os(m_name.c_str());
  // Get xyz coords for all indices
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> x(n);
  std::vector<ergo_real> y(n);
  std::vector<ergo_real> z(n);
  for(int i = 0; i < n; i++) {
    x[i] = basisInfo.basisFuncList[i].centerCoords[0];
    y[i] = basisInfo.basisFuncList[i].centerCoords[1];
    z[i] = basisInfo.basisFuncList[i].centerCoords[2];
  }
  
  size_t number_of_stored_zeros = 0;
  ergo_real minAbsValue = 1e22;
  ergo_real maxAbsValue = 0;

  // Get all matrix elements
  size_t nvalues = 0;
  std::vector<int>       rowind;
  std::vector<int>       colind;
  std::vector<ergo_real> values;
  {
    std::vector<int>       rowind_tmp;
    std::vector<int>       colind_tmp;
    std::vector<ergo_real> values_tmp;
    get_all_nonzeros( A, inversePermutation, rowind_tmp, colind_tmp, values_tmp);
    
    bool matrixIsSymmetric = (A.obj_type_id() == "MatrixSymmetric");
    if (matrixIsSymmetric) {
      // Also include lower triangle
      size_t nvalues_tmp = values_tmp.size();
      rowind.reserve(nvalues_tmp*2);
      colind.reserve(nvalues_tmp*2);
      values.reserve(nvalues_tmp*2);
      for(size_t i = 0; i < nvalues_tmp; i++) {
	rowind.push_back( rowind_tmp[i] );
	colind.push_back( colind_tmp[i] );
	values.push_back( values_tmp[i] );	
	if ( rowind_tmp[i] != colind_tmp[i] ) {
	  rowind.push_back( colind_tmp[i] );
	  colind.push_back( rowind_tmp[i] );
	  values.push_back( values_tmp[i] );	
	}
      } // end for
    } // end if
    else {  
      rowind = rowind_tmp;
      colind = colind_tmp;
      values = values_tmp;
    } // end else
    
    nvalues = values.size();
    // Take absolute value
    for(size_t i = 0; i < nvalues; i++) {
      ergo_real fabsVal = fabs( values[i] );
      values[i] = fabsVal;
      minAbsValue = fabsVal < minAbsValue ? fabsVal : minAbsValue;
      maxAbsValue = fabsVal > maxAbsValue ? fabsVal : maxAbsValue;      
    }
  }
 
  os << "%% Run for example like this: matlab -nosplash -nodesktop -r " <<  name << std::endl;
  os << "number_of_stored_zeros    = " << number_of_stored_zeros << ";" << std::endl;
  os << "number_of_stored_nonzeros = " << nvalues << ";" << std::endl;
  
  // Get distances for all matrix elements  
  std::vector<ergo_real> distances(nvalues);
  for(size_t i = 0; i < nvalues; i++) {
    ergo_real diff_x = x[ rowind[i] ] - x[ colind[i] ];
    ergo_real diff_y = y[ rowind[i] ] - y[ colind[i] ];
    ergo_real diff_z = z[ rowind[i] ] - z[ colind[i] ];
    distances[i] = std::sqrt( diff_x * diff_x + diff_y * diff_y + diff_z * diff_z );
  }

  // Index vector   
  std::vector<size_t> index(nvalues);
  for ( size_t ind = 0; ind < index.size(); ++ind ) {
    index[ind] = ind;
  }
  
  // Sort based on distance
  matrix_utilities_CompareClass<typename std::vector<ergo_real>::const_iterator> 
    compareDist( distances.begin() );
  std::sort ( index.begin(), index.end(), compareDist );
  
  // Min and max distances
  ergo_real minDistance = *std::min_element( distances.begin(), distances.end() );
  ergo_real maxDistance = *std::max_element( distances.begin(), distances.end() );
  // Size of box in r direction
  ergo_real rbox_length = (maxDistance - minDistance) / resolution_r;

  // Get max absolute value of A
  ergo_real maxMagLog10 = std::log10(maxAbsValue);
  ergo_real minMagLog10 = std::log10(minAbsValue) > -20 ? std::log10(minAbsValue) : -20;
  // Size of box in m direction
  ergo_real mbox_length = (maxMagLog10 - minMagLog10) / resolution_m;
  
  os << "A = [ " << std::endl; 
  // Loop over r boxes
  size_t start_ind = 0;
  ergo_real r_low = minDistance;
  for ( int rbox = 0; rbox < resolution_r; rbox++ ) {
    ergo_real r_upp = r_low + rbox_length;
    // Find end index
    size_t end_ind = start_ind;
    while ( end_ind < nvalues && distances[index[end_ind]] < r_upp )
      end_ind++;
    // Now we have the bounds for box in r-direction 
    // Sort based on magnitude
    matrix_utilities_CompareClass<typename std::vector<ergo_real>::const_iterator> 
      compareMagnitude( values.begin() );
    std::sort ( index.begin() + start_ind, index.begin() + end_ind, compareMagnitude );
    // Loop over m boxes
    ergo_real m_low = minMagLog10;
    size_t ind_m = start_ind;

    // Skip very small values
    while ( ind_m < end_ind && std::log10( values[index[ind_m]] ) < m_low ) 
      ind_m++;
    size_t skipped_small = ind_m - start_ind;
    os << r_low << "  " 
       << r_upp << "  " 
       << 0 << "  " 
       << std::pow(10,m_low) << "  " 
       << skipped_small 
       << std::endl;
    
    for ( int mbox = 0; mbox < resolution_m; mbox++ ) {
      ergo_real m_upp = m_low + mbox_length;
      size_t count = 0;
      while ( ind_m < end_ind && std::log10( values[index[ind_m]] ) < m_upp ) {
	ind_m++;
	count++;
      }
      // Now we have r_low r_upp m_low m_upp count
      // Write to stream
      os << r_low << "  " 
		<< r_upp << "  " 
		<< std::pow(10,m_low) << "  " 
		<< std::pow(10,m_upp) << "  " 
		<<  count 
		<< std::endl;
      m_low = m_upp;
    }
    
    r_low = r_upp;
    start_ind = end_ind;
  }
  os << "];" << std::endl; 
  os << "B=[];" << std::endl;
  os << "for ind = 1 : size(A,1)" << std::endl;
  os << "  if (A(ind,3) ~= 0)" << std::endl;
  os << "    B = [B; A(ind,:)];" << std::endl;
  os << "  end" << std::endl;
  os << "end" << std::endl;
  os << "%col = jet(101);" << std::endl;
  os << "col = gray(101);col=col(end:-1:1,:);" << std::endl;
  os << "maxCount = max(B(:,5));" << std::endl;
  os << "ax = [0 30 1e-12 1e3]" << std::endl;

  os << "fighandle = figure;" << std::endl;
  os << "for ind = 1 : size(B,1)" << std::endl;
  os << "  rmin = B(ind, 1); rmax = B(ind, 2);" << std::endl;
  os << "  mmin = B(ind, 3); mmax = B(ind, 4);" << std::endl;
  os << "  colind = round( 1+100 * B(ind,5) / maxCount);" << std::endl;
  os << "  fill([rmin rmin rmax rmax rmin], [mmin mmax mmax mmin mmin], col(colind,:), 'EdgeColor', col(colind,:) )" << std::endl;
  os << "  hold on" << std::endl;
  os << "end" << std::endl;
  os << "set(gca,'YScale','log')" << std::endl;
  os << "axis(ax)" << std::endl;
  os << "set(gca,'FontSize',16)" << std::endl;
  os << "xlabel('Distance')" << std::endl;
  os << "ylabel('Magnitude')" << std::endl;  
  os << "print( fighandle, '-depsc2', '" << name << "')" << std::endl;  

  os << "fighandle = figure;" << std::endl;  
  os << "for ind = 1 : size(B,1)" << std::endl;  
  os << "  if (B(ind,5) ~= 0)" << std::endl;  
  os << "    rmin = B(ind, 1); rmax = B(ind, 2);" << std::endl;
  os << "    mmin = B(ind, 3); mmax = B(ind, 4);" << std::endl;
  os << "    msize = 3+1*ceil(20 * B(ind,5) / maxCount);" << std::endl;
  os << "    plot((rmin+rmax)/2,(mmin+mmax)/2,'k.','MarkerSize',msize)" << std::endl;
  os << "    hold on" << std::endl;
  os << "  end" << std::endl;
  os << "end" << std::endl;
  os << "set(gca,'YScale','log')" << std::endl;
  os << "axis(ax)" << std::endl;
  os << "set(gca,'FontSize',16)" << std::endl;
  os << "xlabel('Distance')" << std::endl;
  os << "ylabel('Magnitude')" << std::endl;
  os << "print( fighandle, '-depsc2', '" << name << "_dots')" << std::endl;
  os << "exit(0);" << std::endl;  
  os.close();
} // end output_distance_vs_magnitude(...)

template<typename Tmatrix>
void output_magnitude_histogram( Tmatrix const & A, /**< The matrix. */
				 std::string name, /**< File name. */
				 int resolution_m  /**< Resolution m-direction, m is the magnitude of matrix elements. */
				 ) {
  std::string m_name = name + ".m";
  std::ofstream os(m_name.c_str());

  size_t number_of_stored_zeros = 0;
  ergo_real minAbsValue = 1e22;
  ergo_real maxAbsValue = 0;

  // Get all matrix elements
  size_t nvalues = 0;
  std::vector<int>       rowind;
  std::vector<int>       colind;
  std::vector<ergo_real> values;
  {
    // Get all nonzeros
    rowind.resize(0);
    colind.resize(0);
    values.resize(0);
    size_t nvalues_tmp = A.nvalues();
    std::vector<int>       rowind_tmp; rowind_tmp.reserve(nvalues_tmp);
    std::vector<int>       colind_tmp; colind_tmp.reserve(nvalues_tmp);
    std::vector<ergo_real> values_tmp; values_tmp.reserve(nvalues_tmp);
    A.get_all_values(rowind_tmp,
		     colind_tmp,
		     values_tmp);
    // Count the number of nonzeros
    for(size_t i = 0; i < nvalues_tmp; i++) {
      nvalues += ( values_tmp[i] != 0 );
    }

    bool matrixIsSymmetric = (A.obj_type_id() == "MatrixSymmetric");
    if (matrixIsSymmetric) {
      // Also include lower triangle
      rowind.reserve(nvalues*2);
      colind.reserve(nvalues*2);
      values.reserve(nvalues*2);
      // Extract all nonzeros
      for(size_t i = 0; i < nvalues_tmp; i++) {
	if ( values_tmp[i] != 0 ) {
	  rowind.push_back( rowind_tmp[i] );
	  colind.push_back( colind_tmp[i] );
	  values.push_back( values_tmp[i] );	
	  if ( rowind_tmp[i] != colind_tmp[i] ) {
	    rowind.push_back( colind_tmp[i] );
	    colind.push_back( rowind_tmp[i] );
	    values.push_back( values_tmp[i] );	
	  }
	}
      }
      nvalues = values.size();
    } // end if
    else {      
      rowind.reserve(nvalues);
      colind.reserve(nvalues);
      values.reserve(nvalues);
      // Extract all nonzeros
      for(size_t i = 0; i < nvalues_tmp; i++) {
	if ( values_tmp[i] != 0 ) {
	  rowind.push_back( rowind_tmp[i] );
	  colind.push_back( colind_tmp[i] );
	  values.push_back( values_tmp[i] );	
	}
      }
      assert( nvalues == values.size() );
    } // end else
    // Take absolute value
    for(size_t i = 0; i < nvalues; i++) {
      ergo_real fabsVal = fabs( values[i] );
      values[i] = fabsVal;
      minAbsValue = fabsVal < minAbsValue ? fabsVal : minAbsValue;
      maxAbsValue = fabsVal > maxAbsValue ? fabsVal : maxAbsValue;      
    }
  }
  
  os << "%% Run for example like this: matlab -nosplash -nodesktop -r " <<  name << std::endl;
  os << "number_of_stored_zeros    = " << number_of_stored_zeros << ";" << std::endl;
  os << "number_of_stored_nonzeros = " << nvalues << ";" << std::endl;

  // Index vector   
  std::vector<size_t> index(nvalues);
  for ( size_t ind = 0; ind < index.size(); ++ind ) {
    index[ind] = ind;
  }
  
  // Get min and max absolute value of A
  ergo_real maxMagLog10 = std::log10(maxAbsValue);
  ergo_real minMagLog10 = std::log10(minAbsValue) > -20 ? std::log10(minAbsValue) : -20;
  // Size of box in m direction
  ergo_real mbox_length = (maxMagLog10 - minMagLog10) / resolution_m;

  os << "A = [ " << std::endl; 
  // Sort based on magnitude
  matrix_utilities_CompareClass<typename std::vector<ergo_real>::const_iterator> 
    compareMagnitude( values.begin() );
  std::sort ( index.begin(), index.end(), compareMagnitude );
  // Loop over m boxes
  ergo_real m_low = minMagLog10;
  size_t ind_m = 0;

  // Skip very small values
  while ( ind_m < nvalues && std::log10( values[index[ind_m]] ) < m_low ) 
    ind_m++;
  size_t skipped_small = ind_m;
  os   << 0 << "  " 
       << std::pow(10,m_low) << "  " 
       << skipped_small 
       << std::endl;
    
  for ( int mbox = 0; mbox < resolution_m; mbox++ ) {
    ergo_real m_upp = m_low + mbox_length;
    size_t count = 0;
    while ( ind_m < nvalues && std::log10( values[index[ind_m]] ) < m_upp ) {
      ind_m++;
      count++;
    }
    // Now we have m_low m_upp count
    // Write to stream
    os << std::pow(10,m_low) << "  " 
       << std::pow(10,m_upp) << "  " 
       <<  count 
       << std::endl;
    m_low = m_upp;
  }
  os << "];" << std::endl; 
  
  os.close();
} // end output_magnitude_histogram

template<typename Tmatrix>
void write_matrix_in_matrix_market_format( Tmatrix const & A, 
					   std::vector<int> const & inversePermutation,
					   std::string filename,
					   std::string identifier,
					   std::string method_and_basis)
{

  // Get all matrix elements
  size_t nvalues = 0;
  std::vector<int>       rowind;
  std::vector<int>       colind;
  std::vector<ergo_real> values;
  get_all_nonzeros( A, inversePermutation, rowind, colind, values);
  nvalues = values.size();
  // Now we have all matrix elements
  // Open file stream
  std::string mtx_filename = filename + ".mtx";
  std::ofstream os(mtx_filename.c_str());

  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  std::string matrix_market_matrix_type = "general";
  bool matrixIsSymmetric = (A.obj_type_id() == "MatrixSymmetric");
  if (matrixIsSymmetric)
    matrix_market_matrix_type = "symmetric";
  os << "%%MatrixMarket matrix coordinate real " << matrix_market_matrix_type << std::endl
     << "%===============================================================================" << std::endl
     << "% Generated by the Ergo quantum chemistry program version " << VERSION << " (www.ergoscf.org)" << std::endl
     << "% Date      : " << asctime (timeinfo) // newline added by asctime
     << "% ID-string : " << identifier << std::endl
     << "% Method    : " << method_and_basis << std::endl
     << "%" << std::endl
     << "% MatrixMarket file format:" << std::endl
     << "% +-----------------" << std::endl
     << "% | % comments" << std::endl
     << "% | nrows ncols nentries" << std::endl
     << "% | i_1        j_1        A(i_1,j_1)" << std::endl
     << "% | i_2        j_2        A(i_1,j_1)" << std::endl
     << "% | ..." << std::endl
     << "% | i_nentries j_nentries A(i_nentries,j_nentries) " << std::endl
     << "% +----------------" << std::endl
     << "% Note that indices are 1-based, i.e. A(1,1) is the first element." << std::endl
     << "%" << std::endl
     << "%===============================================================================" << std::endl;
  os << A.get_nrows() << "  " << A.get_ncols() << "  " << nvalues << std::endl;
  if (matrixIsSymmetric)
    for(size_t i = 0; i < nvalues; i++) {
      // Output lower triangle
      if ( rowind[i] < colind[i] )
	os << colind[i]+1 << "  " << rowind[i]+1 << "  " << std::setprecision(10) << values[i] << std::endl;
      else
	os << rowind[i]+1 << "  " << colind[i]+1 << "  " << std::setprecision(10) << values[i] << std::endl;
    }      
  else
    for(size_t i = 0; i < nvalues; i++) {
      os << rowind[i]+1 << "  " << colind[i]+1 << "  " << std::setprecision(10) << values[i] << std::endl;
    }  
  os.close();
} // end write_matrix_in_matrix_market_format(...)


#endif
