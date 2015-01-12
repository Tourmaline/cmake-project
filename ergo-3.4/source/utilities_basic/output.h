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

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER


#include <stdarg.h>



/* Log categories */
#define LOG_CAT_UNDEFINED   0
#define LOG_CAT_ERROR       1
#define LOG_CAT_WARNING     2
#define LOG_CAT_INFO        3
#define LOG_CAT_EXTRAINFO   4
#define LOG_CAT_RESULTS     5
#define LOG_CAT_TIMINGS     6
#define LOG_CAT_MEMUSAGE    7

/* Log areas */
#define LOG_AREA_UNDEFINED  0
#define LOG_AREA_MAIN       1
#define LOG_AREA_SCF        2
#define LOG_AREA_LR         3
#define LOG_AREA_INTEGRALS  4
#define LOG_AREA_DENSFROMF  5
#define LOG_AREA_DFT        6
#define LOG_AREA_LOWLEVEL   7
#define LOG_AREA_CI         8

/* output functions */
void do_output(int logCategory, int logArea, const char* format, ...);
int do_voutput(int logCategory, int logArea, const char* format, va_list v);
void do_output_time(int logCategory, int logArea, const char* s);
void output_current_memory_usage(int logArea, const char* contextString);
void enable_memory_usage_output(void);
void enable_output();


#endif
