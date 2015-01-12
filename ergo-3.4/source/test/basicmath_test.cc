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

/** @file basicmath_test.cc Tests some basic math functions
    such as template_blas_sqrt() template_blas_log() etc to
    see that they are working properly and that they deliver
    the expected accuracy. */

#include <stdio.h>
#include <stdlib.h>
#include <limits>

#include "realtype.h"
#include "template_blas_common.h"

int main(int argc, char *argv[])
{
  int failed = 0;
  int verbose = getenv("VERBOSE") != NULL;
  ergo_real machine_epsilon = std::numeric_limits<ergo_real>::epsilon();
  
  printf("machine_epsilon = %g Run with env VERBOSE for more info.\n",
         (double)machine_epsilon);


  /* Test sqrt function for a set of random numbers. */
  ergo_real maxabsdiff_sqrt = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real x2 = x * x;
      ergo_real sqrt_of_x2 = template_blas_sqrt(x2);
      ergo_real diff = sqrt_of_x2 - x;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      if(absdiff > maxabsdiff_sqrt)
	maxabsdiff_sqrt = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_sqrt: %g\n",
           (double)maxabsdiff_sqrt);
  ergo_real maxabsdiff_sqrt_requested = machine_epsilon;
  if(maxabsdiff_sqrt > maxabsdiff_sqrt_requested)
    {
      printf("ERROR: template_blas_sqrt() not accurate enough!\n");
      printf("maxabsdiff_sqrt: %g, requested: %g\n", (double)maxabsdiff_sqrt, (double)maxabsdiff_sqrt_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_sqrt() accuracy OK.\n");
  }




  /* Test exp function by computing exp(a)*exp(b) and comparing to exp(a+b) for a list of random pairs (a,b) */
  ergo_real maxabsdiff_exp = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real a = (ergo_real)rand() / RAND_MAX;
      ergo_real b = (ergo_real)rand() / RAND_MAX;
      ergo_real product_of_exps = template_blas_exp(a) * template_blas_exp(b);
      ergo_real exp_of_sum = template_blas_exp(a + b);
      ergo_real diff = product_of_exps - exp_of_sum;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      if(absdiff > maxabsdiff_exp)
	maxabsdiff_exp = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_exp: %g\n", (double)maxabsdiff_exp);
  ergo_real maxabsdiff_exp_requested = machine_epsilon * 15;
  if(maxabsdiff_exp > maxabsdiff_exp_requested)
    {
      printf("ERROR: template_blas_exp() not accurate enough!\n");
      printf("maxabsdiff_exp: %g, requested: %g\n", (double)maxabsdiff_exp, (double)maxabsdiff_exp_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_exp() accuracy OK.\n");
  }





  /* Test log function by computing log(a) + log(b) and comparing to log(a*b) for a list of random pairs (a,b) */
  ergo_real maxabsdiff_log = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real a = (ergo_real)rand() / RAND_MAX;
      ergo_real b = (ergo_real)rand() / RAND_MAX;
      ergo_real sum_of_logs = template_blas_log(a) + template_blas_log(b);
      ergo_real log_of_product = template_blas_log(a * b);
      ergo_real diff = sum_of_logs - log_of_product;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      if(absdiff > maxabsdiff_log)
	maxabsdiff_log = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_log: %g\n", (double)maxabsdiff_log);
  ergo_real maxabsdiff_log_requested = machine_epsilon * 10;
  if(maxabsdiff_log > maxabsdiff_log_requested)
    {
      printf("ERROR: template_blas_log() not accurate enough!\n");
      printf("maxabsdiff_log: %g, requested: %g\n", (double)maxabsdiff_log, (double)maxabsdiff_log_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_log() accuracy OK.\n");
  }




  /* Test erf function by comparing with a series expression */
  ergo_real piBBP = template_blas_compute_pi_BBP((ergo_real)0);
  ergo_real maxabsdiff_erf = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real minus_1_to_pow_n = 1;
      ergo_real n_factorial = 1;
      ergo_real x_to_pow_2n_plus_1 = x;
      int n = 0;
      ergo_real sum = 0;
      while(((ergo_real)1 / n_factorial) > machine_epsilon)
	{
	  sum += (minus_1_to_pow_n / ( n_factorial * (ergo_real)( 2 * n + 1) )) * x_to_pow_2n_plus_1;
	  n++;
	  minus_1_to_pow_n *= -1;
	  n_factorial *= n;
	  x_to_pow_2n_plus_1 *= x * x;
	}
      ergo_real series_result = ((ergo_real)2 / template_blas_sqrt(piBBP)) * sum;
      ergo_real erf_value = template_blas_erf(x);
      ergo_real diff = series_result - erf_value;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      if(absdiff > maxabsdiff_erf)
	maxabsdiff_erf = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_erf: %g\n", (double)maxabsdiff_erf);
  ergo_real maxabsdiff_erf_requested = machine_epsilon * 5;
  if(maxabsdiff_erf > maxabsdiff_erf_requested)
    {
      printf("ERROR: template_blas_erf() not accurate enough!\n");
      printf("maxabsdiff_erf: %g, requested: %g\n", (double)maxabsdiff_erf, (double)maxabsdiff_erf_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_erf() accuracy OK.\n");
  }




  /* Test erfc function by computing erf(x) + erfc(x) and comparing to 1 */
  ergo_real maxabsdiff_erfc = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real erf_of_x = template_blas_erf(x);
      ergo_real erfc_of_x = template_blas_erfc(x);
      ergo_real sum = erf_of_x + erfc_of_x;
      ergo_real diff = sum - (ergo_real)1.0;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      if(absdiff > maxabsdiff_erfc)
	maxabsdiff_erfc = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_erfc: %g\n", (double)maxabsdiff_erfc);
  ergo_real maxabsdiff_erfc_requested = machine_epsilon * 1;
  if(maxabsdiff_erfc > maxabsdiff_erfc_requested)
    {
      printf("ERROR: template_blas_erfc() not accurate enough!\n");
      printf("maxabsdiff_erfc: %g, requested: %g\n", (double)maxabsdiff_erfc, (double)maxabsdiff_erfc_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_erfc() accuracy OK.\n");
  }

  

  if (!failed)
    puts("Basic math functions test succeeded.");
  else
    puts("Basic math functions test FAILED.");

  return failed;
}
