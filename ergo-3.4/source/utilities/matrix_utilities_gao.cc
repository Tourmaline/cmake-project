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
#include <stdio.h>
#include <stdlib.h>
#include "matrix_utilities.h"
#include "output.h"


typedef struct
{
  int row;
  int col;
  ergo_real value;
} CSR_element_struct;

static int
csr_compare_index_pairs_for_qsort(const void* p1, const void* p2)
{
  CSR_element_struct* pair_1 = (CSR_element_struct*)p1;
  CSR_element_struct* pair_2 = (CSR_element_struct*)p2;
  if(pair_1->row > pair_2->row)
    return 1;
  if(pair_1->row < pair_2->row)
    return -1;
  if(pair_1->col > pair_2->col)
    return 1;
  if(pair_1->col < pair_2->col)
    return -1;
  return 0;
}

void write_CSR_matrix(int n,
		      const symmMatrix & M, 
		      const char* fileName, 
		      ergo_real threshold,
		      std::vector<int> const & inversePermutationHML)
{
  // get general matrix
  normalMatrix MM(M);

  // Get all nonzero elements
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  MM.get_all_values(rowind,
		    colind,
		    values,
		    inversePermutationHML,
		    inversePermutationHML);
  int nvalues = values.size();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "write_CSR_matrix, threshold = %7.3g, nvalues = %9i", 
	    (double)threshold, nvalues);


  CSR_element_struct* list = new CSR_element_struct[nvalues];
  int i;
  int count = 0;
  for(i = 0; i < nvalues; i++)
    {
      if(std::fabs(values[i]) > threshold)
	{
	  list[count].row   = rowind[i];
	  list[count].col   = colind[i];
	  list[count].value = values[i];
	  count++;
	}
    }

  // Now all info we need is stored in list
  // Sort list by rows and columns
  qsort(list, count, sizeof(CSR_element_struct), csr_compare_index_pairs_for_qsort);

  int nnz_per_row_vector[n];
  i = 0;
  int row;
  for(row = 0; row < n; row++)
    {
      int nnz_curr_row = 0;
      while(i < count)
	{
	  if(list[i].row == row)
	    {
	      nnz_curr_row++;
	      i++;
	    }
	  else
	    break;
	} // END WHILE
      nnz_per_row_vector[row] = nnz_curr_row;
    }
  if(i != count)
    throw "Error in write_CSR_matrix: (i != count)";

  colind.assign(count, 0);
  rowind.assign(count, 0);
  values.assign(count, 0);
  
  for(i = 0; i < count; i++)
    {
      colind[i] = list[i].col;
      rowind[i] = list[i].row;
      values[i] = list[i].value;
    }

  FILE* f_bin = fopen(fileName, "wb");
  if(f_bin == NULL)
    throw "error in write_CSR_matrix: (f_bin == NULL)";

  char fileNameTXT[888];
  sprintf(fileNameTXT, "%s.txt", fileName);
  FILE* f_txt = fopen(fileNameTXT, "wt");
  if(f_txt == NULL)
    throw "error in write_CSR_matrix: (f_txt == NULL)";

  // Write matrix dimension n
  if(fwrite(&n, sizeof(int), 1, f_bin) != 1)
    throw "error in write_CSR_matrix, in fwrite.";
  fprintf(f_txt, "matrix dimension: %6i\n", n);

  // Write number of nonzeros
  if(fwrite(&count, sizeof(int), 1, f_bin) != 1)
    throw "error in write_CSR_matrix, in fwrite.";
  fprintf(f_txt, "number of nonzero elements: %9i\n", count);

  // Write nnz_per_row_vector
  if(fwrite(nnz_per_row_vector, sizeof(int), n, f_bin) != unsigned(n))
    throw "error in write_CSR_matrix, in fwrite.";
  fprintf(f_txt, "number of nonzeros in each row:\n");
  for(i = 0; i < n; i++)
    fprintf(f_txt, "%6i\n", nnz_per_row_vector[i]);

  // Write column index vector
  if(fwrite(&colind[0], sizeof(int), count, f_bin) != unsigned(count))
    throw "error in write_CSR_matrix, in fwrite.";
  fprintf(f_txt, "column index vector:\n");
  for(i = 0; i < count; i++)
    fprintf(f_txt, "%6i\n", colind[i]);
  
  // Write element value vector
  if(fwrite(&values[0], sizeof(ergo_real), count, f_bin) != unsigned(count))
    throw "error in write_CSR_matrix, in fwrite.";
  fprintf(f_txt, "element values vector:\n");
  for(i = 0; i < count; i++)
    fprintf(f_txt, "%15.8f\n", (double)values[i]);

  fprintf(f_txt, "End of CSR text file.\n");

  fclose(f_bin);
  fclose(f_txt);

  delete list;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "write_CSR_matrix, file '%s' written OK, count = %i", 
	    fileName,
	    count);
}

