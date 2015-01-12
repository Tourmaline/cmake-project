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

#ifndef DENSFROMFFULLHEADER
#define DENSFROMFFULLHEADER

#include "realtype.h"

int get_F_orbs(int n, 
	       const ergo_real* F, 
	       const ergo_real* ovl, 
	       ergo_real* cmo, 
	       ergo_real* eigv);

int get_dens_from_fock_full(int n, 
			    int noOfOccupiedOrbs, 
			    ergo_real* result_P, 
			    const ergo_real* F, 
			    const ergo_real* ovl, 
			    ergo_real factor,
			    ergo_real electronicTemperature,
			    ergo_real & resultEntropyTerm,
			    ergo_real * const lumoVec = 0,
			    ergo_real * const homoVec = 0);

#endif
