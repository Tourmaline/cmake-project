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


typedef double ergo_real;


#define MAX_NO_OF_CONTR_GAUSSIANS 20

#define ergo_malloc malloc
#define ergo_free free



const int DENSITY_FILE_VERSION_NUMBER_OLD = 10001;
const int DENSITY_FILE_VERSION_NUMBER_NEW = 10002;


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
} densityFileHeaderStruct_old;


typedef struct
{
  int densityFileVersion;
  int typeOfMatrixStorage;
  int noOfShells;
  int noOfBasisFuncs;
  int noOfDensityMatrices;
  int dummy_int;
  double matrixSize_1;
  double matrixSize_2;
  double fileSizeInBytes;
} densityFileHeaderStruct_new;





int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      printf("usage: ./x old_density_file new_density_file\n");
      return 0;
    }

  printf("standalone_density_file_conversion_program_2 version 1.1\n");

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
  
  printf("OK, file '%s' seems to contain a proper old-style density description, noOfDensityMatrices = %i\n", 
	 fileName_old, noOfDensityMatrices);
  printf("Now attempting to create new density description file '%s'\n", fileName_new);

  char* p = buffer_old + sizeof(densityFileHeaderStruct_old);

  densityFileHeaderStruct_new densityFileHeader_new;
  densityFileHeader_new.densityFileVersion = DENSITY_FILE_VERSION_NUMBER_NEW;
  densityFileHeader_new.typeOfMatrixStorage = headerPtr_old->typeOfMatrixStorage;
  densityFileHeader_new.noOfShells = headerPtr_old->noOfShells;
  densityFileHeader_new.noOfBasisFuncs = headerPtr_old->noOfBasisFuncs;
  densityFileHeader_new.noOfDensityMatrices = headerPtr_old->noOfDensityMatrices;
  densityFileHeader_new.matrixSize_1 = (double)headerPtr_old->matrixSize_1;
  densityFileHeader_new.matrixSize_2 = (double)headerPtr_old->matrixSize_2;
  densityFileHeader_new.fileSizeInBytes = (double)(headerPtr_old->fileSizeInBytes + sizeof(densityFileHeaderStruct_new) - sizeof(densityFileHeaderStruct_old));
  
  // write to file
  FILE* f_out = fopen(fileName_new, "wb");
  if(f_out == NULL)
    {
      printf("error opening file '%s' for writing\n", fileName_new);
      return -1;
    }
  
  size_t noOfBytesWritten_1 = fwrite(&densityFileHeader_new, 1, sizeof(densityFileHeaderStruct_new), f_out);
  if(noOfBytesWritten_1 != sizeof(densityFileHeaderStruct_new))
    {
      printf("error in fwrite\n");
      return -1;
    }
  int restSize = densityFileHeader_new.fileSizeInBytes - sizeof(densityFileHeaderStruct_new);
  size_t noOfBytesWritten_2 = fwrite(p, 1, restSize, f_out);
  if(noOfBytesWritten_2 != restSize)
    {
      printf("error in fwrite\n");
      return -1;
    }
  fclose(f_out);

  ergo_free(buffer_old);
  
  printf("OK, density description file '%s' created, %i bytes.\n", fileName_new, (int)densityFileHeader_new.fileSizeInBytes);
  
  return 0;
}

