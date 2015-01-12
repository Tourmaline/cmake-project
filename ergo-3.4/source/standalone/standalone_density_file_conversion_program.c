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

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

typedef double ergo_real;


#define MAX_NO_OF_CONTR_GAUSSIANS 20

#define ergo_malloc malloc
#define ergo_free free



const int DENSITY_FILE_VERSION_NUMBER = 10001;


// Matrix storage types
#define MATRIX_STORAGE_TYPE_FULL      1
#define MATRIX_STORAGE_TYPE_TRIANGLE  2



struct ShellSpecStruct_old_{
  int noOfContr;
  ergo_real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real sizeList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real extent;
  ergo_real centerCoords[3]; /* x0, y0, z0 */
  int shellType;
  int shell_ID;
  int noOfBasisFuncs;
  int startIndexInMatrix; /* start index in density matrix  */
};
typedef struct ShellSpecStruct_old_ ShellSpecStruct_old;

typedef struct
{
  int noOfShells;
  int noOfBasisFuncs;
  int noOfDensityMatrices;
  int fileSizeInBytes;
} densityFileHeaderStruct_old;


struct ShellSpecStruct_new_{
  ergo_real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real sizeList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real extent;
  ergo_real centerCoords[3]; /* x0, y0, z0 */
  int noOfContr;
  int shellType;
  int shell_ID;
  int noOfBasisFuncs;
  int startIndexInMatrix; /* start index in density matrix  */
  int dummy; /* padding to make sure the size of this structure is a multiple of 8 bytes */
};
typedef struct ShellSpecStruct_new_ ShellSpecStruct_new;

typedef struct
{
  int densityFileVersion;
  int typeOfMatrixStorage;
  int noOfShells;
  int noOfBasisFuncs;
  int noOfDensityMatrices;
  int matrixSize_1;
  int matrixSize_2;
  int fileSizeInBytes;
} densityFileHeaderStruct_new;


static void
ddf_store_triangular_matrix(char* p, int n, const ergo_real* matrix)
{
  ergo_real* resultPtr = (ergo_real*)p;
  int i, j;
  int count = 0;
  for(i = 0; i < n; i++)
    for(j = i; j < n; j++)
      {
	resultPtr[count] = matrix[i*n+j];
	count++;
      }
}


int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      printf("usage: ./x old_density_file new_density_file\n");
      return 0;
    }

  printf("standalone_density_file_conversion_program version 1.0\n");

  char* fileName_old = argv[1];
  char* fileName_new = argv[2];

  // open old density file
  FILE* f_in = fopen(fileName_old, "rb");
  if(f_in == NULL)
    {
      printf("error opening file '%s'\n", fileName_old);
      return -1;
    }

  // check file size
  fseeko(f_in, 0L, SEEK_END);
  off_t fileSize_old = ftello(f_in);
  if(fileSize_old < sizeof(densityFileHeaderStruct_old))
    {
      printf("error: (fileSize_old < sizeof(densityFileHeaderStruct_old))\n");
      return -1;
    }
  fseeko(f_in, 0L, SEEK_SET);
  char* buffer_old = (char*)ergo_malloc(fileSize_old);
  if(buffer_old == NULL)
    {
      printf("error in malloc.\n");
      return -1;
    }
  size_t readsz_old = fread(buffer_old, 1, fileSize_old, f_in);
  fclose(f_in);
  if(readsz_old != fileSize_old)
    {
      printf("error reading file '%s'\n", fileName_old);
      ergo_free(buffer_old);
      return -1;
    }
  
  // OK, old density file contents are now in buffer_old.
  densityFileHeaderStruct_old* headerPtr_old = (densityFileHeaderStruct_old*)buffer_old;
  if(headerPtr_old->fileSizeInBytes != fileSize_old)
    {
      printf("error: (headerPtr_old->fileSizeInBytes != fileSize_old)\n");
      return -1;
    }
  int noOfDensityMatrices = headerPtr_old->noOfDensityMatrices;
  if(noOfDensityMatrices != 1 && noOfDensityMatrices != 2)
    {
      printf("error: (noOfDensityMatrices != 1 && noOfDensityMatrices != 2)\n");
      return -1;
    }

  int n = headerPtr_old->noOfBasisFuncs;
  int nShells = headerPtr_old->noOfShells;

  printf("nShells = %i, n = %i\n", nShells, n);

  if(n <= 0 || nShells <= 0)
    {
      printf("error: (n <= 0 || nShells <= 0)\n");
      return -1;
    }

  int fileSize_expected_old = 
    sizeof(densityFileHeaderStruct_old) + 
    nShells * sizeof(ShellSpecStruct_old) +
    noOfDensityMatrices * n * n * sizeof(ergo_real);
  if(fileSize_expected_old != fileSize_old)
    {
      printf("error: (fileSize_expected_old != fileSize_old)");
      return -1;
    }

  printf("OK, file '%s' seems to contain a proper old-style density description, noOfDensityMatrices = %i\n", 
	 fileName_old, noOfDensityMatrices);
  printf("Now attempting to create new density description file '%s'\n", fileName_new);

  char* p = buffer_old + sizeof(densityFileHeaderStruct_old);

  ShellSpecStruct_old* shellList_old = (ShellSpecStruct_old*)p;
  p += nShells * sizeof(ShellSpecStruct_old);

  ergo_real* matrixPtr_1_old = NULL;
  ergo_real* matrixPtr_2_old = NULL;
  matrixPtr_1_old = (ergo_real*)p;
  p += n*n*sizeof(ergo_real);
  if(noOfDensityMatrices == 2)
    {
      matrixPtr_2_old = (ergo_real*)p;
      p += n*n*sizeof(ergo_real);
    }
  if((p - buffer_old) != fileSize_old)
    {
      printf("error: ((p - buffer_old) != fileSize_old)\n");
      return -1;
    }

  int tringularMatrixStorageSize = sizeof(ergo_real) * n * (n+1) / 2;
  int fileSize_new = 
    sizeof(densityFileHeaderStruct_new) +  
    nShells * sizeof(ShellSpecStruct_new) + 
    noOfDensityMatrices * tringularMatrixStorageSize;
  
  // allocate buffer for new density description
  char* buffer_new = (char*)malloc(fileSize_new);
  if(buffer_new == NULL)
    {
      printf("error in malloc.\n");
      return -1;
    }
  memset(buffer_new, 0, fileSize_new);
  p = buffer_new;

  densityFileHeaderStruct_new densityFileHeader_new;
  densityFileHeader_new.densityFileVersion = DENSITY_FILE_VERSION_NUMBER;
  densityFileHeader_new.typeOfMatrixStorage = MATRIX_STORAGE_TYPE_TRIANGLE;
  densityFileHeader_new.noOfShells = nShells;
  densityFileHeader_new.noOfBasisFuncs = n;
  densityFileHeader_new.noOfDensityMatrices = noOfDensityMatrices;
  densityFileHeader_new.matrixSize_1 = tringularMatrixStorageSize;
  if(noOfDensityMatrices == 2)
    densityFileHeader_new.matrixSize_2 = tringularMatrixStorageSize;
  else
    densityFileHeader_new.matrixSize_2 = 0;
  densityFileHeader_new.fileSizeInBytes = fileSize_new;

  // do header
  memcpy(p, &densityFileHeader_new, sizeof(densityFileHeaderStruct_new));
  p += sizeof(densityFileHeaderStruct_new);

  // do shell list
  ShellSpecStruct_new* shellList_new = (ShellSpecStruct_new*)p;
  int i;
  for(i = 0; i < nShells; i++)
    {
      int j;
      for(j = 0; j < MAX_NO_OF_CONTR_GAUSSIANS; j++)	
	{
	  shellList_new[i].coeffList[j]    = shellList_old[i].coeffList[j];
	  shellList_new[i].exponentList[j] = shellList_old[i].exponentList[j];
	  shellList_new[i].sizeList[j]     = shellList_old[i].sizeList[j];
	}
      shellList_new[i].extent              = shellList_old[i].extent;

      for(j = 0; j < 3; j++)
	shellList_new[i].centerCoords[j]   = shellList_old[i].centerCoords[j];

      shellList_new[i].noOfContr           = shellList_old[i].noOfContr;
      shellList_new[i].shellType           = shellList_old[i].shellType;
      shellList_new[i].shell_ID            = shellList_old[i].shell_ID;
      shellList_new[i].noOfBasisFuncs      = shellList_old[i].noOfBasisFuncs;
      shellList_new[i].startIndexInMatrix  = shellList_old[i].startIndexInMatrix;
      shellList_new[i].dummy = 0;
    } // END FOR i
  p += nShells*sizeof(ShellSpecStruct_new);

  // do density matrices
  ddf_store_triangular_matrix(p, n, matrixPtr_1_old);
  p += tringularMatrixStorageSize;
  if(noOfDensityMatrices == 2)
    {
      ddf_store_triangular_matrix(p, n, matrixPtr_2_old);
      p += tringularMatrixStorageSize;
    }
  
  if((p - buffer_new) != fileSize_new)
    {
      printf("error: ((p - buffer_new) != fileSize_new)\n");
      return -1;
    }

  // OK, buffer_new is complete, time to write it to file.
  
  // write to file
  FILE* f_out = fopen(fileName_new, "wb");
  if(f_out == NULL)
    {
      printf("error opening file '%s' for writing\n", fileName_new);
      return -1;
    }

  int noOfBytesWritten = (int)fwrite(buffer_new, 1, fileSize_new, f_out);
  fclose(f_out);
  if(noOfBytesWritten != fileSize_new)
    {
      printf("error in fwrite, fileSize_new = %i, noOfBytesWritten = %i",
	     fileSize_new, noOfBytesWritten);
      return -1;
    }

  ergo_free(buffer_old);
  ergo_free(buffer_new);

  printf("OK, density description file '%s' created, %i bytes.\n", fileName_new, fileSize_new);
  
  return 0;
}

