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

/** @file tddft_test.cc Tests the TDDFT interface.

It has to produce number of files for given molecule and basis set:

a. one electron integral matrix (T+V)
b. g_pqrs
c. V_pqrs^{xc}(rho0)

Example usage is:

 source/tddft/tddft_test mol/h2.mol STO-3G "Combine hf=0"


*/

#include <stdio.h>
#include <memory>

#include "Matrix.h"
#include "SCF_restricted.h"
#include "basisinfo.h"
#include "density_description_file.h"
#include "dft_common.h"
#include "integral_info.h"
#include "integrals_2el.h"
#include "matrix_utilities.h"
#include "memorymanag.h"
#include "molecule.h"
#include "scf.h"
#include "tddft.h"

static char usage[] =
  "Usage: tddft_test MOLFILE BASIS FUNCTIONAL\n"
  "\n"
  "Generates the relevant integrals for the TDDFT calculation.\n"
  "The data is saved in output.m containing data in octave format.\n";

int
main(int argc, char *argv[])
{
  if(argc<3) {
    fputs(usage, stderr);
    return 1;
  }

  Molecule molecule;
  Molecule extraCharges; /* Not used in this test. */

  if(molecule.setFromMoleculeFile(argv[1], 0, NULL)) {
    fprintf(stderr, "Reading molecule from %s failed.\n", argv[1]);
    return 1;
  }
  const char *basisFileName = argv[2];


  IntegralInfo integralInfo(true);
  BasisInfoStruct basisInfo;
  if(basisInfo.addBasisfuncsForMolecule(molecule, 
                                        basisFileName,
                                        0, NULL,
                                        integralInfo, 
                                        false,
                                        true, true) != 0) {
    fprintf(stderr, "Error in BasisInfoStruct::add_basisfuncs_for_molecule "
            "for main basis set, Basis='%s'",
            basisFileName);
    return -1;
  }

  SCF::Options scfOptions;  /* Defaults */
  scfOptions.use_dft  = strcmp(argv[3], "HF") != 0;
  if(scfOptions.use_dft) {
    if(dft_setfunc(argv[3]) == 0) {
      fprintf(stderr, "Error in functional definition '%s'.", argv[3]);
      return 1;
    }
  }

  try {
    static const ergo_real THRESHOLD_1EL    = 1e-12;
    JK::Params   jkOptions;   /* Defaults */
    SCF::MatOptions matOpts;  /* Defaults */
    matOpts.prepare(basisInfo);

    BasisInfoStruct basisInfoDensFit;
    Dft::GridParams gridParams;
  
    SCF_restricted scf(molecule,
		       extraCharges,
                       basisInfo, 
                       basisInfoDensFit,
                       integralInfo,
                       NULL,
                       jkOptions,
		       gridParams,
                       scfOptions,
                       matOpts,
                       THRESHOLD_1EL);
    scf.do_SCF_iterations();
    FILE *f = fopen("output.m", "wt");
    TDDFT::saveOverlap(basisInfo, f);
    TDDFT::saveDipole(basisInfo, f);
    TDDFT::saveKinetic(basisInfo, f);
    TDDFT::savePotential(molecule, basisInfo, integralInfo, f);
    TDDFT::saveCoulomb(basisInfo, integralInfo, f);
    {
      ergo_real       *densMatrix_full = NULL;
      BasisInfoStruct *basis_read  = NULL;
      if(ddf_load_density("density.bin", 1, integralInfo,
			  &basis_read, &densMatrix_full) == 0) {
        TDDFT::saveXC(molecule, basisInfo, densMatrix_full, f);
	size_t n = basis_read->noOfBasisFuncs;
	TDDFT::writeMatlab(f, densMatrix_full, n, "density");
      } else
        fprintf(stderr, "ERROR: Density file disappeared?\n");
      ergo_free(densMatrix_full);
      delete basis_read;
    }
    fclose(f);
  } catch (const std::exception& e) {
    fprintf(stderr, "Exception caught: %s\n", e.what());
    return 1;
  } catch (const char* s) {
    fprintf(stderr, "Exception caught: %s\n", s);
    return 1;
  }
  puts("\nTD-DFT data generation completed.");
  return 0;
}
