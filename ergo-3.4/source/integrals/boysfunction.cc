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

#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "boysfunction.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"


static double
semiFactorial(int n)
{
  switch(n)
    {
    case -1:
      return 1;
    case 0:
      return 1;
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 3;
    case 4:
      return 4*2;
    case 5:
      return 5*3;
    case 6:
      return 6*4*2;
    case 7:
      return 7*5*3;
    case 8:
      return 8*6*4*2;
    case 9:
      return 9*7*5*3;
    case 10:
      return 10*8*6*4*2;
    case 11:
      return 11*9*7*5*3;
    case 12:
      return 12*10*8*6*4*2;
    case 13:
      return 13*11*9*7*5*3;
    case 14:
      return 14*12*10*8*6*4*2;
    case 15:
      return 15*13*11*9*7*5*3;
    case 16:
      return 16*14*12*10*8*6*4*2;
    case 17:
      return 17*15*13*11*9*7*5*3;
    case 18:
      return 18*16*14*12*10*8*6*4*2;
    case 19:
      return 19*(double)17*15*13*11*9*7*5*3;
    case 20:
      return 20*(double)18*16*14*12*10*8*6*4*2;
    case 21:
      return 21*19*(double)17*15*13*11*9*7*5*3;
    case 22:
      return 22*20*(double)18*16*14*12*10*8*6*4*2;
    case 23:
      return 23*21*19*(double)17*15*13*11*9*7*5*3;
    case 24:
      return 24*22*20*(double)18*16*14*12*10*8*6*4*2;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: semiFactorial not implemented for n > 24");
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "n = %i", n);
      exit(0);
      return 0;
    }
}


#if 0
static ergo_real
BoysFunction_raw(int n, ergo_real x)
{
  const int N = 400;
  int i;
  ergo_real h = (ergo_real)0.5 / N;
  ergo_real sum = 0;
  for(i = 0; i < N; i++)
    {
      ergo_real t = (ergo_real)i / N + h;
      sum += 2 * h * std::exp(-x*t*t) * std::pow(t, 2*n);
    }
  return sum;
}
#endif


static ergo_real
BoysFunction_raw_simpson(int n, ergo_real x)
{
  const int N = 100;
  ergo_real h = (ergo_real)0.5 / N;
  ergo_real sum = 0;
  for(int k = 0; k <= 2*N; k++)
    {
      ergo_real tk = (ergo_real)k / (2*N);
      // Compute f(tk) = exp(-x*tk*tk) * pow(tk, 2*n)
      ergo_real foftk = std::exp(-x*tk*tk);
      if(n != 0)
	{
	  if(k != 0)
	    foftk *= std::pow(tk, 2*n);
	  else
	    foftk = 0;
	}
      // OK, foftk done, now add to sum.
      if(k == 0 || k == 2*N)
	{
	  sum += foftk;
	  continue;
	}
      if(k % 2 == 1)
	{
	  sum += 4 * foftk;
	  continue;
	}
      sum += 2 * foftk;
    }
  return (h/3) * sum;
}




void
BoysFunctionManager::init(void)
{
  if(Boys_init_flag == 1)
    return;
  Util::TimeMeter timeMeter;
  ergo_real halfstep, kfactorial, BoysFuncRawResult, Ak, midx;
  halfstep = (ergo_real)BOYS_X_MAX / BOYS_NO_OF_INTERVALS * 0.5;
  for(int n = 0; n < BOYS_N_MAX; n++)
    {
      for(int j = 0; j < BOYS_NO_OF_INTERVALS; j++)
	{
	  midx = (ergo_real)BOYS_X_MAX * j / BOYS_NO_OF_INTERVALS + halfstep;
	  Boys_list[n][j].midx = midx;
	  kfactorial = 1;
	  int minusOneToPowk = 1;
	  for(int k = 0; k < BOYS_TAB_DEGREE; k++)
	    {
	      BoysFuncRawResult = BoysFunction_raw_simpson(n+k, midx);
	      Ak = minusOneToPowk * BoysFuncRawResult / kfactorial;
	      Boys_list[n][j].A[k] = Ak;
	      kfactorial *= k+1;
	      minusOneToPowk *= -1;
	    } /* END FOR k */
	} /* END FOR j */
    } /* END FOR n */
  Boys_init_flag = 1;
  timeMeter.print(LOG_AREA_INTEGRALS, "BoysFunctionManager::init");
}

ergo_real
BoysFunctionManager::BoysFunction_pretabulated(int n, ergo_real x) const
{
  const BoysFuncIntervalStruct* interval;
  ergo_real intervalWidth, count, sum, deltax, deltaxtopowk;
  int intervalIndex, k;
  if(Boys_init_flag != 1)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (Boys_init_flag != 1).");
  if(x < 0)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (x < 0).");
  if(n >= BOYS_N_MAX)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (n >= BOYS_N_MAX).");
  if(x >= BOYS_X_MAX) {
    /* use "large x formula" */
    return (semiFactorial(2*n-1) / std::pow((ergo_real)2, n+1)) * std::sqrt(pi / std::pow(x, 2*n+1));
  }
  /* choose which interval to use */
  intervalWidth = (ergo_real)BOYS_X_MAX / BOYS_NO_OF_INTERVALS;
  count = x / intervalWidth;
  intervalIndex = (int)std::floor(count);
  if((intervalIndex < 0) || (intervalIndex >= BOYS_NO_OF_INTERVALS))
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: bad intervalIndex.");
  interval = &Boys_list[n][intervalIndex];
  sum = 0;
  deltax = x - interval->midx;
  deltaxtopowk = 1;
  for(k = 0; k < BOYS_TAB_DEGREE; k++)
    {
      ergo_real Ak = interval->A[k];
      sum += Ak * deltaxtopowk;
      deltaxtopowk *= deltax;
    }
  return sum;
}

ergo_real
BoysFunctionManager::BoysFunction(int n, ergo_real x) const {
  return BoysFunction_pretabulated(n, x);
}

/** Function needed for Chunks&Tasks usage. */
void BoysFunctionManager::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error: bufferSize too small.");
  // Boys_list
  memcpy(p, Boys_list, sizeof(Boys_list));
  p += sizeof(Boys_list);
  // Boys_init_flag
  memcpy(p, &Boys_init_flag, sizeof(int));
  p += sizeof(int);
  // DONE!  
}

/** Function needed for Chunks&Tasks usage. */
size_t BoysFunctionManager::get_size() const {
  return sizeof(Boys_list) + sizeof(int);
}

/** Function needed for Chunks&Tasks usage. */
void BoysFunctionManager::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  // Boys_list
  memcpy(Boys_list, p, sizeof(Boys_list));
  p += sizeof(Boys_list);
  // Boys_init_flag
  memcpy(&Boys_init_flag, p, sizeof(int));
  p += sizeof(int);
  // DONE!
  if(static_cast<size_t>(p-dataBuffer) > bufferSize)
    throw std::runtime_error("Error: (p > bufferSize).");  
}
