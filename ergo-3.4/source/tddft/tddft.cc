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

/** Interface from ERGO to TD-DFT routines. */

#include <string.h>

#include "dft_common.h"
#include "integrator.h"
#include "integrals_1el_kinetic.h"
#include "integrals_1el_potential.h"
#include "integrals_2el_explicit.h"
#include "operator_matrix.h"
#include "tddft.h"
#include "grid_matrix.h"

BEGIN_NAMESPACE(TDDFT);

static const ergo_real THRESHOLD = 1e-15;

/** Writes specified quadratic matrix to specified file in matlab
    format.  Returns 0 on success, -1 on failure. */
int
writeMatlab(FILE *f, const ergo_real *mat, int n, const char *matName)
{
  if(fprintf(f, "%s = [\n", matName) < 1) return -1;
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      if(fprintf(f, " %lg", (double)mat[j + i*n]) < 1) return -1;
    if(fputs(";\n", f) == EOF) return -1;
  }
  if(fputs("];\n", f) == EOF) return -1;
  return 0;
}

#if 0
/** Writes specified quadratic matrix to specified file in binary
    format. Returns 0 on success, other number on failure. */
static int
writeBinary(FILE *f, const ergo_real *mat, int n)
{
  return fwrite(mat, sizeof(ergo_real), n*n, f) - n*n;
}
#endif

/** Saves one-electron part of the KS matrix to given file. */
int
savePotential(const Molecule& m, const BasisInfoStruct& bis,
	      const IntegralInfo& ii, FILE *f)
{
  int retval;
  int n = bis.noOfBasisFuncs;
  ergo_real *mat= new ergo_real[n*n];

  memset(mat, 0, n*n*sizeof(ergo_real));
  int res = compute_V_matrix_full(bis, ii, m.getNoOfAtoms(),
                                  m.getAtomListPtr(),
                                  THRESHOLD, mat);
  if(res != 0)
    throw "Error in tddft savePotential, in compute_V_matrix_full.";

  /* Save matrix here. */
  retval = writeMatlab(f, mat, n, "potential");

  delete []mat;
  return retval;
}

/** Saves the kinetic energy matrix. */
int
saveKinetic(const BasisInfoStruct& bis, FILE *f)
{
  int retval;
  int n = bis.noOfBasisFuncs;
  ergo_real *mat= new ergo_real[n*n];

  memset(mat, 0, n*n*sizeof(ergo_real));
  int res = compute_T_matrix_full(bis, THRESHOLD, mat);
  if(res != 0)
    throw "Error in tddft saveKinetic, in compute_T_matrix_full.";

  /* Save matrix here. */
  retval = writeMatlab(f, mat, n, "kinetic");

  delete []mat;
  return retval;
}

/** Saves the overlap matrix. */

int
saveOverlap(const BasisInfoStruct& bis,  FILE *f)
{
  int retval = -1;
  unsigned n = bis.noOfBasisFuncs;
  ergo_real *mat= new ergo_real[n*n];

  if( (retval = compute_overlap_matrix(bis, bis, mat)) == 0) {
    retval = writeMatlab(f, mat, n, "overlap");
  }
  delete []mat;

  return retval;
}

/** Saves the dipole matrix to specified file.  */

int
saveDipole(const BasisInfoStruct& bis, FILE *f)
{
  int comp;
  unsigned n = bis.noOfBasisFuncs;
  
  ergo_real *mat= new ergo_real[n*n];

  for(comp=0; comp<3; comp++) {
    int d[3];
    d[0] = d[1] = d[2]  = 0;
    d[comp] = 1;

    if(compute_operator_matrix_full(bis, bis, d[0], d[1], d[2], mat))
      break;
    char matName[200];
    sprintf(matName, "dipole(1:%d,1:%d,%d)", n, n, comp+1);
    if(writeMatlab(f, mat, n, matName) != 0)
      break;
  }

  delete []mat;

  return comp == 3 ? 0 : -1;
}


int
saveCoulomb(const BasisInfoStruct& bis,
            const IntegralInfo& ii,  FILE *f)
{
  unsigned n = bis.noOfBasisFuncs;
  unsigned p, q, r, s;
  
  ergo_real *mat = new ergo_real[n*n];

  for(p = 0; p < n; p++)
    for(q = 0; q < n; q++) {
      for(r = 0; r < n; r++)
	for(s = 0; s < n; s++)
          mat[s + r*n] = do_2e_integral(p, q, r, s, 
                                        bis, ii);

      char matName[200];
      sprintf(matName, "g(1:%d,1:%d,%d,%d)", n, n, p+1, q+1);
      if(writeMatlab(f, mat, n, matName) != 0) {
        delete []mat;
        return -1;
      }
    }

  delete []mat;

  return 0;
}

/* =================================================================== */
/* Exchange-correlation section.                                       */
/* Start from LDA type-functional handling...                          */
static void
hessianCb(DftIntegratorBl* grid, real *tmp,
          int bllen, int blstart, int blend,
          void* cb_data)
{
  int p, q, r, s, k;
  ergo_real *hessian = (ergo_real*)cb_data;
  FunDensProp dp = { 0 };

  for(k=blstart; k<blend; k++) {
    
    SecondDrv vxc;
    real weight = grid->weight[grid->curr_point+k];
    dp.rhoa = dp.rhob = 0.5*grid->r.rho[k];
    dftpot1_(&vxc, &weight, &dp, &ZEROI);
    tmp[k] = vxc.fRR*2;
  }
  static const int SYMMETRY = 0;
  const ergo_real *aos = grid->atv;

  int n = grid->nbast;
  int (*const blocks)[2] = BASBLOCK(grid,SYMMETRY);
  int blCnt = grid->bas_bl_cnt[SYMMETRY];
            
  for(int pBl=0; pBl<blCnt; pBl++)
    for(p=blocks[pBl][0]; p<blocks[pBl][1]; p++) { 
      const ergo_real *pOrbs = aos + p*bllen;

      for(int qBl=0; qBl<blCnt; qBl++)
        for(q=blocks[qBl][0]; q<blocks[qBl][1]; q++) { 
          const ergo_real *qOrbs = aos + q*bllen;

          for(int rBl=0; rBl<blCnt; rBl++)
            for(r=blocks[rBl][0]; r<blocks[rBl][1]; r++) { 
              const ergo_real *rOrbs = aos + r*bllen;

              for(int sBl=0; sBl<blCnt; sBl++)
                for(s=blocks[sBl][0]; s<blocks[sBl][1]; s++) { 
                  const ergo_real *sOrbs = aos + s*bllen;
                  ergo_real *hessianPQRS =
                    hessian + s + n*(r + n*(q + n*p));
                  for(k=blstart; k<blend; k++)
                    *hessianPQRS +=
                      pOrbs[k]*qOrbs[k]*rOrbs[k]*sOrbs[k]*tmp[k];
                }
            }
        }
    }
}

int saveXC(const Molecule& molecule, const BasisInfoStruct& bis,
           const ergo_real* dMat,  FILE *f)
{
  int n = bis.noOfBasisFuncs;
  ergo_real *hessian = new ergo_real[n*n*n*n];
  int ret = 0;
  Dft::GridParams gss;
  Dft::FullMatrix density(dMat, n);
  const Dft::FullMatrix *densPtr = &density;
  memset(hessian, 0, n*n*n*n*sizeof(ergo_real));

  ergo_real nElectrons =
    Dft::integrate(1, &densPtr, bis, molecule, gss, 1, hessianCb, hessian);

  fprintf(stderr, 
          "Hessian integration got %lf electrons in ground state density\n",
          (double)nElectrons);

  for(int p=0; p<n; p++)
    for(int q=0; q<n; q++) {
      char matName[200];
      sprintf(matName, "xc(1:%d,1:%d,%d,%d)", n, n, p+1, q+1);
      if(writeMatlab(f, hessian + n*n*(p + n*q), n, matName) != 0) {
	ret = -1;
	break;
      }
    }

  delete []hessian;
  return ret;
}


END_NAMESPACE(TDDFT);
