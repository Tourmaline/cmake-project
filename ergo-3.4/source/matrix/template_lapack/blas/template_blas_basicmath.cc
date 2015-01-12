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
 
 /* This file belongs to the template_lapack part of the Ergo source 
  * code. The source files in the template_lapack directory are modified
  * versions of files originally distributed as CLAPACK, see the
  * Copyright/license notice in the file template_lapack/COPYING.
  */
 

/* We need to include config.h to get macros HAVE_SQRTL etc. */
#include "config.h"

#include <cmath>
#include <math.h>
#include <stdio.h>

#include "template_blas_basicmath.h"

/* fabs function */

#ifdef HAVE_FABSF
template<> float template_blas_fabs<float>(float x) { return fabsf(x); }
#else
template<> float template_blas_fabs<float>(float x) { return fabs(x); }
#endif

template<> double template_blas_fabs<double>(double x) { return fabs(x); }

#ifdef HAVE_FABSL
template<> long double template_blas_fabs<long double>(long double x) { return fabsl(x); }
#else
template<> long double template_blas_fabs<long double>(long double x) { return fabs(x); }
#endif



/* sqrt function */

#ifdef HAVE_SQRTF
template<> float template_blas_sqrt<float>(float x) { return sqrtf(x); }
#else
template<> float template_blas_sqrt<float>(float x) { return sqrt(x); }
#endif

template<> double template_blas_sqrt<double>(double x) { return sqrt(x); }

#ifdef HAVE_SQRTL
template<> long double template_blas_sqrt<long double>(long double x) { return sqrtl(x); }
#else
template<> long double template_blas_sqrt<long double>(long double x) { return sqrt(x); }
#endif



/* exp function */

#ifdef HAVE_EXPF
template<> float template_blas_exp<float>(float x) { return expf(x); }
#else
template<> float template_blas_exp<float>(float x) { return exp(x); }
#endif

template<> double template_blas_exp<double>(double x) { return exp(x); }

#ifdef HAVE_EXPL
template<> long double template_blas_exp<long double>(long double x) { return expl(x); }
#else
template<> long double template_blas_exp<long double>(long double x) { return exp(x); }
#endif



/* log function */

#ifdef HAVE_LOGF
template<> float template_blas_log<float>(float x) { return logf(x); }
#else
template<> float template_blas_log<float>(float x) { return log(x); }
#endif

template<> double template_blas_log<double>(double x) { return log(x); }

#ifdef HAVE_LOGL
template<> long double template_blas_log<long double>(long double x) { return logl(x); }
#else
template<> long double template_blas_log<long double>(long double x) { return log(x); }
#endif



/* error function erf */

#ifdef HAVE_ERFF
template<> float template_blas_erf<float>(float x) { return erff(x); }
#else
template<> float template_blas_erf<float>(float x) { return erf(x); }
#endif

template<> double template_blas_erf<double>(double x) { return erf(x); }

#ifdef HAVE_ERFL
template<> long double template_blas_erf<long double>(long double x) { return erfl(x); }
#else
template<> long double template_blas_erf<long double>(long double x) { return erf(x); }
#endif



/* complementary error function erfc */

#ifdef HAVE_ERFCF
template<> float template_blas_erfc<float>(float x) { return erfcf(x); }
#else
template<> float template_blas_erfc<float>(float x) { return erfc(x); }
#endif

template<> double template_blas_erfc<double>(double x) { return erfc(x); }

#ifdef HAVE_ERFCL
template<> long double template_blas_erfc<long double>(long double x) { return erfcl(x); }
#else
template<> long double template_blas_erfc<long double>(long double x) { return erfc(x); }
#endif




