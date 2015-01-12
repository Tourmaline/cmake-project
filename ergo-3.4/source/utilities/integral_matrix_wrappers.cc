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

#include <string.h>

#include "integral_matrix_wrappers.h"
#include "output.h"
#include "basis_func_pair_list_1el.h"
#include "integrals_1el_kinetic.h"
#include "integrals_1el_potential.h"
#include "operator_matrix.h"
#include "basis_func_pair_list.h"
#include "integrals_2el_boxed.h"
#include "integrals_2el_coulomb.h"
#include "integrals_2el_exchange.h"
#include "integrals_2el_exchange_prep.h"
#include "utilities.h"
#include "memorymanag.h"


static ergo_real get_max_charge(const Molecule& molecule) {
  ergo_real maxCharge = 0;
  for(int i = 0; i < molecule.getNoOfAtoms(); i++) {
    ergo_real currCharge = molecule.getAtom(i).charge;
    if(currCharge > maxCharge)
      maxCharge = currCharge;
  }
  return maxCharge;
}


int
compute_V_sparse(const BasisInfoStruct& basisInfo,
		 const IntegralInfo& integralInfo,
		 const Molecule& molecule,
		 ergo_real threshold,
		 ergo_real boxSize,
		 symmMatrix & V,
		 std::vector<int> const & permutationHML)
{
  int n = basisInfo.noOfBasisFuncs;

  ergo_real maxCharge = get_max_charge(molecule);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "compute_V_sparse, no of atoms = %5i, maxCharge = %6.2f", 
	    molecule.getNoOfAtoms(), (double)maxCharge);

  int noOfBasisFuncIndexPairs = get_basis_func_pair_list_1el(basisInfo,
							     threshold,
							     maxCharge,
							     NULL,
							     2000000000);
  if(noOfBasisFuncIndexPairs <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_basis_func_pair_list_1el, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
      return -1;
    }
  std::vector<basis_func_index_pair_struct_1el> basisFuncIndexPairList(noOfBasisFuncIndexPairs);
  std::vector<ergo_real> V_list(noOfBasisFuncIndexPairs);
  noOfBasisFuncIndexPairs = get_basis_func_pair_list_1el(basisInfo,
							 threshold,
							 maxCharge,
							 &basisFuncIndexPairList[0],
							 noOfBasisFuncIndexPairs);
  if(noOfBasisFuncIndexPairs <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_basis_func_pair_list_1el, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
      return -1;
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "noOfBasisFuncIndexPairs = %i ==> storing %6.2f %% of a full matrix", 
	    noOfBasisFuncIndexPairs, (double)100*noOfBasisFuncIndexPairs/((double)n*n));

  if(compute_V_and_gradient_linear(basisInfo,
				   integralInfo,
				   molecule,
				   threshold,
				   boxSize,
				   &basisFuncIndexPairList[0],
				   &V_list[0],
				   noOfBasisFuncIndexPairs,
				   false, NULL, NULL) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_V_and_gradient_linear");
      return -1;
    }
  
  // Now transfer result from V_list to V
  {
    std::vector<int> rowind(noOfBasisFuncIndexPairs);
    std::vector<int> colind(noOfBasisFuncIndexPairs);
    std::vector<ergo_real> V_list_Tmp(&V_list[0], &V_list[noOfBasisFuncIndexPairs]);
    for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
      {
	rowind[i] = basisFuncIndexPairList[i].index_1;
	colind[i] = basisFuncIndexPairList[i].index_2;
      }
    V.assign_from_sparse(rowind,
			 colind,
			 V_list_Tmp,
			 permutationHML,
			 permutationHML);
  }

  return 0;
}


int
compute_gradient_of_nucl_and_trDV(const BasisInfoStruct& basisInfo,
				  const IntegralInfo& integralInfo,
				  const Molecule& molecule,
				  ergo_real threshold,
				  ergo_real boxSize,
				  const symmMatrix & densityMatrix_sparse,
				  std::vector<int> const & permutationHML,
				  ergo_real* result_gradient_list) {
  int n = basisInfo.noOfBasisFuncs;

  ergo_real maxCharge = get_max_charge(molecule);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "compute_gradient_of_tr_D_V, no of atoms = %5i, maxCharge = %6.2f", 
	    molecule.getNoOfAtoms(), (double)maxCharge);

  int noOfBasisFuncIndexPairs = get_basis_func_pair_list_1el(basisInfo, threshold, maxCharge, NULL, 2000000000);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_basis_func_pair_list_1el, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
    return -1;
  }
  std::vector<basis_func_index_pair_struct_1el> basisFuncIndexPairList(noOfBasisFuncIndexPairs);
  std::vector<ergo_real> V_list(noOfBasisFuncIndexPairs);
  noOfBasisFuncIndexPairs = get_basis_func_pair_list_1el(basisInfo, threshold, maxCharge, &basisFuncIndexPairList[0], noOfBasisFuncIndexPairs);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_basis_func_pair_list_1el, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
    return -1;
  }
  
  // Get density matrix elements
  std::vector<ergo_real> D_list(noOfBasisFuncIndexPairs);
  {
    std::vector<int> rowind(noOfBasisFuncIndexPairs);
    std::vector<int> colind(noOfBasisFuncIndexPairs);
    for(int i = 0; i < noOfBasisFuncIndexPairs; i++) {
      rowind[i] = basisFuncIndexPairList[i].index_1;
      colind[i] = basisFuncIndexPairList[i].index_2;
    }
    densityMatrix_sparse.get_values(rowind, colind, D_list, permutationHML, permutationHML);
  }
  
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "noOfBasisFuncIndexPairs = %i ==> storing %6.2f %% of a full matrix", 
	    noOfBasisFuncIndexPairs, (double)100*noOfBasisFuncIndexPairs/((double)n*n));

  if(compute_V_and_gradient_linear(basisInfo,
				   integralInfo,
				   molecule,
				   threshold,
				   boxSize,
				   &basisFuncIndexPairList[0],
				   &V_list[0],
				   noOfBasisFuncIndexPairs,
				   true, &D_list[0],
				   result_gradient_list) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_V_and_gradient_linear");
      return -1;
    }
  
  return 0;
}


ergo_real
get_electron_nuclear_attraction_energy(const IntegralInfo& integralInfo,
				       const Molecule& molecule,
				       const BasisInfoStruct& basisInfo,
				       const symmMatrix & D,
				       ergo_real threshold_integrals_1el,
				       mat::SizesAndBlocks const & matrix_size_block_info,
				       std::vector<int> const & permutationHML)
{
  // Get V
  symmMatrix V;
  V.resetSizesAndBlocks(matrix_size_block_info, matrix_size_block_info);
  ergo_real boxSize = 10.0;
  if(compute_V_sparse(basisInfo,
		      integralInfo, 
		      molecule,
		      threshold_integrals_1el,
		      boxSize,
		      V,
		      permutationHML) != 0)
    throw "error in compute_V_sparse";
  return symmMatrix::trace_ab(D, V);
}


int
compute_T_sparse(const BasisInfoStruct& basisInfo,
		 const IntegralInfo& integralInfo,
		 ergo_real threshold,
		 symmMatrix & T,
		 std::vector<int> const & permutationHML) {
  int n = basisInfo.noOfBasisFuncs;
  std::vector<int> nvaluesList(n);
  std::vector<int*> colindList(n);
  std::vector<ergo_real*> valuesList(n);
  if(compute_T_matrix_sparse(basisInfo,
			     threshold,
			     n,
			     &nvaluesList[0],
			     &colindList[0],
			     &valuesList[0]) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_T_matrix_sparse");
    return -1;
  }
  // Now convert result to three vectors so prepare for HML storage.
  int nvalues = 0;
  for(int i = 0; i < n; i++)
    nvalues += nvaluesList[i];
  // allocate vectors
  std::vector<int> rowind(nvalues);
  std::vector<int> colind(nvalues);
  std::vector<ergo_real> values(nvalues);
  // populate vectors
  int count = 0;
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < nvaluesList[i]; j++)
	{
	  rowind[count] = i;
	  colind[count] = colindList[i][j];
	  values[count] = valuesList[i][j];
	  count++;
	} // END FOR j
    } // END FOR i
  // Now the information is in rowind colind values.
  // free memory allocated by compute_operator_matrix_sparse.
  for(int i = 0; i < n; i++) {
    ergo_free(colindList[i]);
    ergo_free(valuesList[i]);
  }
  T.assign_from_sparse(rowind,
		       colind,
		       values,
		       permutationHML,
		       permutationHML);
  return 0;
}



static int
check_diagonal_elements_of_overlap_matrix(int n, const symmMatrix & S_symm)
{
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  std::vector<ergo_real> values(n);
  for(int i = 0; i < n; i++)
    {
      rowind[i] = i;
      colind[i] = i;
    }
  S_symm.get_values(rowind,
		    colind,
		    values);
  ergo_real maxabsdiff = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real absdiff = std::fabs(values[i] - 1);
      if(absdiff > maxabsdiff)
	maxabsdiff = absdiff;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "diag elements of ovl matrix differ from 1 by at most %9.6g",
            (double)maxabsdiff);
  if(maxabsdiff > 1e-3)
    {
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "WARNING: bad overlap matrix: diag elements of ovl matrix differ from 1 by more than safe limit.");
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "WARNING: bad overlap matrix: diag elements of ovl matrix differ from 1 by more than safe limit.");
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "WARNING: bad overlap matrix: diag elements of ovl matrix differ from 1 by more than safe limit.");
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "WARNING: this means that basis functions are not completely normlized.");
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "WARNING: this means that basis functions are not completely normlized.");
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "WARNING: this means that basis functions are not completely normlized.");
    }
  return 0;
}


int
compute_overlap_matrix_sparse(const BasisInfoStruct& basisInfo,
			      symmMatrix & S_symm,
			      std::vector<int> const & permutationHML)
{
  if(compute_operator_matrix_sparse_symm(basisInfo,
					 0,
					 0,
					 0,
					 S_symm,
					 permutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_operator_matrix_sparse_symm");
      return -1;
    }
  
  // check diagonal elements of overlap matrix
  int n = basisInfo.noOfBasisFuncs;
  check_diagonal_elements_of_overlap_matrix(n, S_symm);
  return 0;
}



typedef int* int_ptr;
typedef ergo_real* ergo_real_ptr;

int
compute_operator_matrix_sparse_symm(const BasisInfoStruct& basisInfo,
				    int pow_x,
				    int pow_y,
				    int pow_z,
				    symmMatrix & A_symm,
				    std::vector<int> const & permutationHML)
{
  int n = basisInfo.noOfBasisFuncs;

  std::vector<int> nvaluesList(n);
  std::vector< std::vector<int> > colindList(n);
  std::vector< std::vector<ergo_real> > valuesList(n);
  if(compute_operator_matrix_sparse(basisInfo,
				    basisInfo,
				    pow_x,
				    pow_y,
				    pow_z,
				    n,
				    n,
				    nvaluesList,
				    colindList,
				    valuesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_operator_matrix_sparse");
      return -1;
    }
  
  // Now convert result to three vectors so prepare for HML storage.
  int nvalues = 0;
  for(int i = 0; i < n; i++)
    nvalues += nvaluesList[i];

  // allocate vectors
  std::vector<int> rowind(nvalues);
  std::vector<int> colind(nvalues);
  std::vector<ergo_real> values(nvalues);
  // populate vectors
  int count = 0;
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < nvaluesList[i]; j++)
	{
	  rowind[count] = i;
	  colind[count] = colindList[i][j];
	  values[count] = valuesList[i][j];
	  count++;
	} // END FOR j
    } // END FOR i

  normalMatrix A_norm(A_symm);

  A_norm.assign_from_sparse(rowind,
			    colind,
			    values,
			    permutationHML,
			    permutationHML);

  // Convert to symmetric form.
  A_symm = A_norm;
  A_norm.clear();

  // Earlier, the matrix was truncated here but we changed this so it is truncated outside this routine.

  // Write to file and read again to reduce memory fragmentation.
  A_symm.writeToFile();
  A_symm.readFromFile();

  return 0;
}




int
compute_J_by_boxes_sparse(const BasisInfoStruct& basisInfo,
			  const IntegralInfo& integralInfo, 
			  const JK::Params& J_K_params,
			  symmMatrix & J,
			  const symmMatrix & densityMatrix_sparse,
			  std::vector<int> const & permutationHML)
{
  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "compute_J_by_boxes_sparse() start.");
  Util::TimeMeter timeMeterWriteAndReadAll;
  std::string sizesStr = mat::FileWritable::writeAndReadAll();
  timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll: '" + sizesStr).c_str());

  int n = basisInfo.noOfBasisFuncs;

  ergo_real maxDensityMatrixElement = densityMatrix_sparse.maxAbsValue();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_J_by_boxes_sparse, maxDensityMatrixElement = %22.11f", 
	    (double)maxDensityMatrixElement);

  std::vector<basis_func_index_pair_struct> basisFuncIndexPairList;
  int noOfBasisFuncIndexPairs = get_basis_func_pair_list_2el(basisInfo,
							     integralInfo,
							     J_K_params.threshold_J,
							     maxDensityMatrixElement,
							     basisFuncIndexPairList);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, 
	      "error in get_basis_func_pair_list, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
    return -1;
  }
  
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "compute_J_by_boxes_sparse: noOfBasisFuncIndexPairs = %i ==> storing %6.2f %% of a full matrix", 
	    noOfBasisFuncIndexPairs, (double)100*noOfBasisFuncIndexPairs/((double)n*n));
  
  // Setup D_list
  std::vector<ergo_real> D_list(noOfBasisFuncIndexPairs);
  {
    std::vector<int> rowind(noOfBasisFuncIndexPairs);
    std::vector<int> colind(noOfBasisFuncIndexPairs);
    for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
      {
	rowind[i] = basisFuncIndexPairList[i].index_1;
	colind[i] = basisFuncIndexPairList[i].index_2;
      }
    densityMatrix_sparse.get_values(rowind,
				    colind,
				    D_list,
				    permutationHML,
				    permutationHML);
  }

  std::vector<ergo_real> J_list(noOfBasisFuncIndexPairs);

  output_current_memory_usage(LOG_AREA_SCF, "Before calling compute_J_by_boxes_linear");
  if(compute_J_by_boxes_linear(basisInfo,
			       integralInfo,
			       J_K_params,
			       &basisFuncIndexPairList[0],
			       noOfBasisFuncIndexPairs,
			       &D_list[0],
			       &J_list[0],
			       noOfBasisFuncIndexPairs) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_J_by_boxes_linear");
      return -1;
    }
  output_current_memory_usage(LOG_AREA_SCF, "After calling compute_J_by_boxes_linear");

  // Now transfer result from J_list to J
  {
    std::vector<int> rowind(noOfBasisFuncIndexPairs);
    std::vector<int> colind(noOfBasisFuncIndexPairs);
    for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
      {
	rowind[i] = basisFuncIndexPairList[i].index_1;
	colind[i] = basisFuncIndexPairList[i].index_2;
      }
    J.assign_from_sparse(rowind,
			 colind,
			 J_list,
			 permutationHML,
			 permutationHML);
    output_current_memory_usage(LOG_AREA_SCF, "After J.assign_from_sparse");
  }

  timeMeter.print(LOG_AREA_SCF, "compute_J_by_boxes_sparse total");
  return 0;
}



static int 
get_CSR_from_symmMatrix(int n,
			const symmMatrix & A,
			std::vector<int> const & inversePermutationHML,
			csr_matrix_struct & CSR)
{
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  A.get_all_values(rowind,
		   colind,
		   values,
		   inversePermutationHML,
		   inversePermutationHML);
  int nvalues = values.size();      
  // switch rows and columns if necessary
  for(int i = 0; i < nvalues; i++)
    {
      int row = rowind[i];
      int col = colind[i];
      if(row > col)
	{
	  rowind[i] = col;
	  colind[i] = row;
	}
    }      
  if(ergo_CSR_create(&CSR, 
		     1,
		     n,
		     nvalues,
		     &rowind[0],
		     &colind[0]) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_create.");
      return -1;
    }
  for(int i = 0; i < nvalues; i++)
    {
      ergo_CSR_add_to_element(&CSR, 
			      rowind[i],
			      colind[i],
			      values[i]);
    }
  // check that CSR matrix is correct.
  for(int i = 0; i < nvalues; i++)
    {
      if(ergo_CSR_get_element(&CSR, rowind[i], colind[i]) != values[i])
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error checking CSR matrix.");
	  return -1;
	}
    }
  return 0;
}


static int 
get_CSR_from_normalMatrix(int n,
			  const normalMatrix & A,
			  std::vector<int> const & inversePermutationHML,
			  csr_matrix_struct & CSR)
{
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  A.get_all_values(rowind,
		   colind,
		   values,
		   inversePermutationHML,
		   inversePermutationHML);
  int nvalues = values.size();      
  if(ergo_CSR_create(&CSR, 
		     0,
		     n,
		     nvalues,
		     &rowind[0],
		     &colind[0]) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_create.");
      return -1;
    }
  for(int i = 0; i < nvalues; i++) {
    ergo_CSR_add_to_element(&CSR, 
			    rowind[i],
			    colind[i],
			    values[i]);
  }
  // check that CSR matrix is correct.
  for(int i = 0; i < nvalues; i++) {
    if(ergo_CSR_get_element(&CSR, rowind[i], colind[i]) != values[i])
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "Error in get_CSR_from_normalMatrix, error checking CSR matrix.");
	return -1;
      }
  }
  return 0;
}



/** Returns the exchange matrix multiplied by 0.5.
 * To get the correct value multiply K by 2.
 */
int
compute_K_by_boxes_sparse(const BasisInfoStruct& basisInfo,
			  const IntegralInfo& integralInfo, 
			  const JK::ExchWeights & CAM_params,
			  const JK::Params& J_K_params,
			  symmMatrix & K,
			  symmMatrix & densityMatrix_sparse,
			  std::vector<int> const & permutationHML,
			  std::vector<int> const & inversePermutationHML)
{
  Util::TimeMeter timeMeter;
  int n = basisInfo.noOfBasisFuncs;
  output_current_memory_usage(LOG_AREA_SCF, "Beginning of compute_K_by_boxes_sparse");
  
  csr_matrix_struct dens_CSR;
  memset(&dens_CSR, 0, sizeof(csr_matrix_struct));

  csr_matrix_struct K_CSR;
  memset(&K_CSR, 0, sizeof(csr_matrix_struct));

  if(get_CSR_from_symmMatrix(n, 
			     densityMatrix_sparse,
			     inversePermutationHML,
			     dens_CSR) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_CSR_from_symmMatrix");
      return -1;
    }

  // set up CSR K matrix for result
  // note that this call will allocate memory, we must call ergo_CSR_destroy later.
  if(create_CSR_for_K(basisInfo,
		      integralInfo,
		      J_K_params,
		      &dens_CSR,
		      &K_CSR,
		      1) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes");
      return -1;
    }

  densityMatrix_sparse.writeToFile();
  
  output_current_memory_usage(LOG_AREA_SCF, "Before calling compute_K_by_boxes");
  if(compute_K_by_boxes(basisInfo,
			integralInfo,
			CAM_params,
			J_K_params,
			NULL,
			&K_CSR,
			NULL,
			&dens_CSR,
			1) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes");
      return -1;
    }
  output_current_memory_usage(LOG_AREA_SCF, "After calling compute_K_by_boxes");

  // collect result
  // Convert K matrix from CSR to HML format.
  int nvalues = ergo_CSR_get_nvalues(&K_CSR);
  std::vector<int> rowind(nvalues);
  std::vector<int> colind(nvalues);
  std::vector<ergo_real> values(nvalues);
  if(ergo_CSR_get_values(&K_CSR,
			 &rowind[0], 
			 &colind[0], 
			 &values[0],
			 nvalues) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_get_values.");
      return -1;
    }
  // now the information is in the vectors rowind colind values, we can free memory for K_CSR.

  ergo_CSR_destroy(&K_CSR);
  K.assign_from_sparse(rowind,
		       colind,
		       values,
		       permutationHML,
		       permutationHML);

  ergo_CSR_destroy(&dens_CSR);
  
  densityMatrix_sparse.readFromFile();
  
  timeMeter.print(LOG_AREA_SCF, "compute_K_by_boxes_sparse total");
  return 0;
}



#if 0
static void
get_fullmatrix_from_normalMatrix(int n,
                                 const normalMatrix & A,
                                 std::vector<int> const & inversePermutationHML,
                                 ergo_real* full)
{
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  A.get_all_values(rowind,
                   colind,
                   values,
                   inversePermutationHML,
                   inversePermutationHML);
  int nvalues = values.size();
  for(int i = 0; i < n*n; i++)
    full[i] = 0;
  for(int i = 0; i < nvalues; i++) {
    int p = rowind[i];
    int q = colind[i];
    full[p*n+q] = values[i];
  }
}
#endif


int
compute_K_by_boxes_sparse_nosymm(const BasisInfoStruct& basisInfo,
				 const IntegralInfo& integralInfo, 
				 const JK::ExchWeights & CAM_params,
				 const JK::Params& J_K_params,
				 normalMatrix & K,
				 normalMatrix & densityMatrix_sparse,
				 std::vector<int> const & permutationHML,
				 std::vector<int> const & inversePermutationHML)
{
  Util::TimeMeter timeMeter;
  int n = basisInfo.noOfBasisFuncs;
  output_current_memory_usage(LOG_AREA_SCF, "Beginning of compute_K_by_boxes_sparse");
  
  csr_matrix_struct dens_CSR;
  memset(&dens_CSR, 0, sizeof(csr_matrix_struct));

  csr_matrix_struct K_CSR;
  memset(&K_CSR, 0, sizeof(csr_matrix_struct));

  csr_matrix_struct K_CSR_2;
  memset(&K_CSR, 0, sizeof(csr_matrix_struct));

  if(get_CSR_from_normalMatrix(n, 
			       densityMatrix_sparse,
			       inversePermutationHML,
			       dens_CSR) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_CSR_from_normalMatrix");
      return -1;
    }

  // set up CSR K matrix for result
  // note that this call will allocate memory, we must call ergo_CSR_destroy later.
  if(create_CSR_for_K(basisInfo,
		      integralInfo,
		      J_K_params,
		      &dens_CSR,
		      &K_CSR,
		      0) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes");
      return -1;
    }
  if(create_CSR_for_K(basisInfo,
		      integralInfo,
		      J_K_params,
		      &dens_CSR,
		      &K_CSR_2,
		      0) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes");
      return -1;
    }

  densityMatrix_sparse.writeToFile();

  output_current_memory_usage(LOG_AREA_SCF, "Before calling compute_K_by_boxes");
  if(compute_K_by_boxes(basisInfo,
			integralInfo,
			CAM_params,
			J_K_params,
			NULL,
			&K_CSR,
			NULL,
			&dens_CSR,
			0) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes");
      return -1;
    }
  output_current_memory_usage(LOG_AREA_SCF, "After calling compute_K_by_boxes");

  // collect result
  // Convert K matrix from CSR to HML format.
  int nvalues = ergo_CSR_get_nvalues(&K_CSR);
  std::vector<int> rowind(nvalues);
  std::vector<int> colind(nvalues);
  std::vector<ergo_real> values(nvalues);
  if(ergo_CSR_get_values(&K_CSR,
			 &rowind[0], 
			 &colind[0], 
			 &values[0],
			 nvalues) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in ergo_CSR_get_values.");
      return -1;
    }
  // now the information is in the vectors rowind colind values, we can free memory for K_CSR.

  ergo_CSR_destroy(&K_CSR);
  K.assign_from_sparse(rowind,
		       colind,
		       values,
		       permutationHML,
		       permutationHML);

  ergo_CSR_destroy(&dens_CSR);
  
  densityMatrix_sparse.readFromFile();
  
  timeMeter.print(LOG_AREA_SCF, "compute_K_by_boxes_sparse total");
  return 0;
}
