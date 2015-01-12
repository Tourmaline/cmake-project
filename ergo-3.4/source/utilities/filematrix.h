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

#ifndef FILE_MATRIX_HEADER
#define FILE_MATRIX_HEADER

#include "matrix_typedefs.h"
#include "realtype.h"

#include <exception>
#include "Failure.h"

template<typename Treal, typename TMatrixType>
  class FileMatrix 
{
 public:
  FileMatrix();
  ~FileMatrix();
  TMatrixType M;
  void AssignFileName(const char* fileName);
  void ReadFromFile();
  void WriteToFile();
 private:
  char* FileName;
};

template<typename Treal, typename TMatrixType>
  FileMatrix<Treal, TMatrixType>::FileMatrix()
{
  FileName = NULL;
}

template<typename Treal, typename TMatrixType>
  FileMatrix<Treal, TMatrixType>::~FileMatrix()
{
  if(FileName != NULL)
    {
      // Remove file if it exists
      unlink(FileName);
      // free memory used for file name string.
      delete FileName;
      FileName = NULL;
    }
}

template<typename Treal, typename TMatrixType>
  void FileMatrix<Treal, TMatrixType>::AssignFileName(const char* fileName)
{
  if(fileName != NULL)
    {
      FileName = new char[strlen(fileName)+1];
      strcpy(FileName, fileName);
    }
  else
    throw "Error in FileMatrix::AssignFileName: fileName == NULL";
}

template<typename Treal, typename TMatrixType>
  void FileMatrix<Treal, TMatrixType>::WriteToFile()
{
  if(FileName != NULL)
    {
      int noOfBytes = 0;
      M.write_to_buffer_count(noOfBytes);
      Treal* buffer = (ergo_real*)new char[noOfBytes];
      M.write_to_buffer(buffer, noOfBytes);
      FILE* f = fopen(FileName, "wb");
      if(f == NULL)
	throw "error in FileMatrix::WriteToFile, in fopen";
      if(fwrite(buffer, sizeof(char), noOfBytes, f) != (unsigned int)noOfBytes)
	throw "error in FileMatrix::WriteToFile, in fwrite";
      fclose(f);
      delete buffer;
      // Free memory used by matrix.
      M.clear();
    }
}

template<typename Treal, typename TMatrixType>
  void FileMatrix<Treal, TMatrixType>::ReadFromFile()
{
  if(FileName != NULL)
    {
      // open file
      FILE* f = fopen(FileName, "rb");
      if(f == NULL)
	throw "error in FileMatrix::ReadFromFile, in fopen";
      // get file size
      fseek(f, 0L, SEEK_END);
      int fileSize = ftell(f);
      fseek(f, 0L, SEEK_SET);
      if(fileSize <= 0)
	throw "error in FileMatrix::ReadFromFile, (fileSize <= 0)";
      // allocate buffer
      char* buffer = new char[fileSize];
      // read file
      if(fread(buffer, sizeof(char), fileSize, f) != (unsigned int)fileSize)
	throw "error in FileMatrix::ReadFromFile, in fread";
      // close file
      fclose(f);
      // Create matrix
      M.read_from_buffer(buffer, fileSize);
      delete buffer;
    }
}


#endif

