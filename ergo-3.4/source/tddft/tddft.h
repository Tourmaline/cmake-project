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

/** @file tddft.h Provides integral evaluation and export routines.
    
The main usage is is to perform the explicitly time-dependent HF/DFT
calculations.

*/

#include "molecule.h"
#include "matrix_typedefs.h"

#if !defined(BEGIN_NAMESPACE)
#define BEGIN_NAMESPACE(x) namespace x {
#define END_NAMESPACE(x)   }; /* x */
#endif

BEGIN_NAMESPACE(TDDFT);

int writeMatlab(FILE *f, const ergo_real *mat, int n, const char *matName);

int savePotential(const Molecule& m, const BasisInfoStruct& bis,
		  const IntegralInfo& ii, FILE *f);

int saveKinetic(const BasisInfoStruct& bis, FILE *f);
int saveOverlap(const BasisInfoStruct& bis, FILE *f);
int saveDipole(const BasisInfoStruct& bis, FILE *f);

int saveCoulomb(const BasisInfoStruct& bis,
                const IntegralInfo& ii, FILE *f);

int saveXC(const Molecule& m, const BasisInfoStruct& bis,
           const ergo_real *densityMatrix_full, FILE *f);

END_NAMESPACE(TDDFT);
