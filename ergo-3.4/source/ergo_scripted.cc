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

/** @file ergo_scripted.cc

    \brief The main program for the ergo project.  It enables
    scripting and more complex input forms.

    @author: Pawel Salek <em>responsible</em>. But feel free to modify
    the file if you are humbly convinced your ideas are correct.
*/
/* Copyright(c) Pawel Salek 2006. */

#include <dirent.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory>
#include <string>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/times.h>
#include <errno.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atom_labels.h"
#include "density_description_file.h"
#include "ergo_scripted.h"
#include "grid_reader.h"
#include "dft_common.h"
#include "lin_trans.h"
#include "integrals_2el.h"
#include "integrals_2el_explicit.h"
#include "integrals_2el_boxed.h"
#include "integrals_2el_exchange.h"
#include "integrals_2el_coulomb.h"
#include "integrals_general.h"
#include "operator_matrix.h"
#include "memorymanag.h"
#include "molecule.h"
#include "output.h"
#include "scf.h"
#include "scf_utils.h"
#include "slr.h"
#include "matrix_utilities.h"
#include "SCF_restricted.h"
#include "SCF_unrestricted.h"
#include "units.h"
#include "ci.h"
#include "license.h"
#include "xyz_file_parser.h"
#include "electron_dynamics.h"
#include "tdhf_dynamics.h"

// ELIAS NOTE 2014-07-14: define SKIP_UNOFFICIAL_INPUT_PARAMS for "official" releases of the code, so that some testing/debugging parameters are skipped.
#define SKIP_UNOFFICIAL_INPUT_PARAMS

static void variable_free(struct variable* v);

/** An object representing the state of the input processor.  A way to
 * initialize state and to cleanly shut it down and release memory
 * is provided. */
class Ergo {
public:
  static const int NO_OF_BASIS_SET_RANGES = 3;

  struct variable* var_list;
  struct variable* J_K_params;
  struct variable* lr_params;
  struct variable* ed_params;
  struct variable* mat_params;
  struct variable* scf_params;
  struct variable* XC_params;
  struct variable* output_params;

  Molecule molecule;
  Molecule ghostMolecule;
  Molecule extraChargesMolecule;
  ergo_real moleculeUnit; /**< the distance unit for inline molecule
                             input. */
  enum MolType readingMoleculeClass; /**< tells which inline molecule we are
                                        reading now: main or ghost. */

  JK::Params   jkOptions;
  SCF::Options    scfOptions;
  SCF::MatOptions matOptions;
  ED::Params edOptions; /* Electron dynamics (ED) options. */

  void registerInputVariables();
  char *Basis;            /**< name of the current basis set. */
  char *GhostBasis;       /**< name of the ghost basis set. */

  BasissetNameRange basissetRangeList[NO_OF_BASIS_SET_RANGES];
  BasissetNameRange basissetRangeListGhost[NO_OF_BASIS_SET_RANGES];

  Ergo() : Basis(NULL), GhostBasis(NULL) {
    memset(basissetRangeList, 0, NO_OF_BASIS_SET_RANGES * sizeof(BasissetNameRange));
    memset(basissetRangeListGhost, 0, NO_OF_BASIS_SET_RANGES * sizeof(BasissetNameRange));
  }
  ~Ergo() {
    if(Basis)
      ergo_free(Basis);
    if(GhostBasis)
      ergo_free(GhostBasis);
    variable_free(var_list); /* This one owns the data. Other ones are
                              * just helpers... */
  }
};

static Ergo ergo;

/** Molecule stores geometry of the current molecule. */

static IntegralInfo* ergoIntegralInfo = NULL;
static BasisInfoStruct* Basis_info = NULL;
/* End of static variable block. */

/** Macro for compact expression of recognized keywords. We make some
    effort to convert all the floating-point default values to double
    type so that they can be passed through the stack without
    problem. The only potential problem is a potential loss of
    precision if sizeof(ergo_real) > sizeof(double) but this we can
    hopefully live with for input variables, can we? */
#define KW(kl,vname, type, defval, desc)                          \
  kl = variable_new_ ##type(kl, (#vname), (desc), (type), (defval))
#define variable_new_VAR_STRING variable_new
#define variable_new_VAR_FLOAT(kl,n,h,t,v) variable_new(kl,n,h,t,double(v))
#define variable_new_VAR_INT    variable_new
#define variable_new_VAR_LIST   variable_new

/** creates new variable item. Such variable can be later assigned
    values etc.

    @param tail is a tail of the variable list, allowing easy variable
    list creation.
    @param name is the variable name.
    @param description is a string with a few sentences describing what the variable is for.
    @param type is the variable type (string, int, or float).
*/
static struct variable*
variable_new(struct variable* tail, const char *name, const char *description,
             enum VarType type, ...)
{
  struct variable * v = ergo_new(1, struct variable);
  va_list ap;

  va_start(ap, type);
  v->next = tail;
  v->name = name;
  v->description = description;
  v->type = type;
  switch(type) {
  case VAR_STRING: v->v.str  = strdup(va_arg(ap, char *)); break;
  case VAR_FLOAT:  v->v.num  = va_arg(ap, double);         break;
  case VAR_INT:    v->v.vint = va_arg(ap, int);            break;
  case VAR_LIST:   v->v.list = va_arg(ap, struct variable*); break;
  default: do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "unknown variable type %d\n", type);
  }
  va_end(ap);
  return v;
}

/** release variable data structure and its children. */
static void
variable_free(struct variable* v)
{
  switch(v->type) {
  case VAR_STRING: free(v->v.str); break;
  case VAR_FLOAT:                  break;
  case VAR_INT:                    break;
  case VAR_LIST:   variable_free(v->v.list); break;
  default: do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
                     "unknown variable type %d\n", v->type);
  }
  if(v->next)
    variable_free(v->next);
  ergo_free(v);
}

/** es_assign_num assigns given numerical value to the variable. */
void
es_assign_num(struct variable *v, double val)
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, " %s := %g\n", v->name, val);
  switch(v->type) {
  case VAR_FLOAT:  v->v.num  =      val; break;
  case VAR_INT:    v->v.vint = (int)val; break;
  default: 
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
	      "Assignment of numerical value to nonnumerical "
	      "variable %s ignored.\n", v->name);
  }
}

/** es_assign_str assigns given string to the variable. It
    additionally clears some local variables if a value is assigned to
    one of the "special" variables like "output_basis or "basis". */
void
es_assign_str(struct variable *v, const char *str)
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, " %s := %s\n", v->name, str);
  if(v->v.str) free(v->v.str);
  v->v.str = strdup(str);
  /* FIXME: move the following code to a modify callback. */
  if( Basis_info && (strcmp(v->name, "output_basis") == 0 ||
		     strcmp(v->name, "basis") == 0 ||
		     strcmp(v->name, "ghost_basis") == 0) ) {
    delete Basis_info;
    Basis_info = NULL;
  }
  if(strcmp(v->name, "basis") == 0) {
    if(ergo.Basis) free(ergo.Basis);
    ergo.Basis = strdup(str);
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Basis really set to %s\n",
              ergo.Basis);
  } else if(strcmp(v->name, "ghost_basis") == 0) {
    if(ergo.GhostBasis) free(ergo.GhostBasis);
    ergo.GhostBasis = strdup(str);
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Ghost basis really set to %s\n",
              ergo.GhostBasis);
  }
}

/** Defines a range of atoms that will be assigned a specific basis set.
    The range can be reset by specifying a empty count.

    @param mt chooses main or the ghost molecule.
    @param rangeNo choose the range to be assigned (1 to 3).
    @param start the start index.
    @param cnt the count of atoms in the range.
    @param name the name of the basis set file.
*/
int
es_assign_range(MolType mt, int rangeNo,
                int start, int cnt, const char *name)
{
  if(rangeNo <1 || rangeNo > Ergo::NO_OF_BASIS_SET_RANGES)
    return false;
  BasissetNameRange *bnrs;
  switch(mt) {
  case MOL_MAIN:  bnrs = ergo.basissetRangeList; break;
  case MOL_GHOST: bnrs = ergo.basissetRangeListGhost; break;
  default: return false;
  }
  --rangeNo;
  printf("Assigning range %s %d [%d:%d] = %s\n",
         mt == MOL_MAIN ? "MAIN" : "GHOST",
         rangeNo, start, start+cnt-1, name);
  bnrs[rangeNo].startAtomIndex = start;
  bnrs[rangeNo].count = cnt;
  if (bnrs[rangeNo].basisSetFileName)
    free(bnrs[rangeNo].basisSetFileName);
  if(name && *name)
    bnrs[rangeNo].basisSetFileName = strdup(name);
  else
    bnrs[rangeNo].basisSetFileName = NULL;
  return true;
}

/** finds the variable struct by @param name starting in the specified
    root. @param root must be of type VAR_LIST. */
struct variable*
es_find_var(struct variable *root, const char *name)
{
  struct variable *res;
  if(root) {
    if(root->type == VAR_LIST)
      res = root->v.list;
    else res = root;
  } else res = ergo.var_list;

  const char* dot_pos = strchr(name, '.');
  size_t l = dot_pos ? dot_pos-name : strlen(name);
  
  while(res && (strncmp(res->name, name, l) || l != strlen(res->name)) ) {
    res = res->next;
  }
  if(!res) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "Variable %s not found.\n", name);
  }
  if(res && dot_pos) {
    return (res->type == VAR_LIST)
      ? es_find_var(res->v.list, dot_pos+1) : NULL;
  }
  else return res;
}

static inline int
var_get_int_template(struct variable *root, const char *name) 
{
  struct variable *v = es_find_var(root, name);
  if(v && v->type == VAR_INT)
    return v->v.vint;
  else {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Undefined integer variable %s\n", name);
    return -12345;
  }
}
static inline double
var_get_real_template(struct variable *root, const char *name) 
{
  struct variable *v = es_find_var(root, name);
  if (v && v->type == VAR_FLOAT)
    return v->v.num;
  else {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Undefined real variable: %s\n", name);
    return -123456.0;
  }
}

static inline const char*
var_get_string(struct variable *root, const char *name) 
{
  struct variable *v = es_find_var(root, name);
  return (v && v->type == VAR_STRING) ? v->v.str : "";
}

#define var_get_int(n)    var_get_int_template(NULL,       (n))
#define var_get_intJK(n)  var_get_int_template(ergo.J_K_params, (n))
#define var_get_intLR(n)  var_get_int_template(ergo.lr_params, (n))
#define var_get_intED(n)  var_get_int_template(ergo.ed_params, (n))
#define var_get_intMA(n)  var_get_int_template(ergo.mat_params, (n))
#define var_get_intOU(n)  var_get_int_template(ergo.output_params, (n))
#define var_get_intSCF(n) var_get_int_template(ergo.scf_params, (n))
#define var_get_intXC(n)  var_get_int_template(ergo.XC_params, (n))
#define var_get_real(n)    var_get_real_template(NULL,       (n))
#define var_get_realJK(n)  var_get_real_template(ergo.J_K_params, (n))
#define var_get_realLR(n)  var_get_real_template(ergo.lr_params, (n))
#define var_get_realED(n)  var_get_real_template(ergo.ed_params, (n))
#define var_get_realMA(n)  var_get_real_template(ergo.mat_params, (n))
#define var_get_realOU(n)  var_get_real_template(ergo.output_params, (n))
#define var_get_realSCF(n) var_get_real_template(ergo.scf_params, (n))
#define var_get_realXC(n)  var_get_real_template(ergo.XC_params, (n))
#define var_get_stringLR(n)  var_get_string(ergo.lr_params, (n))
#define var_get_stringED(n)  var_get_string(ergo.ed_params, (n))
#define var_get_stringSCF(n) var_get_string(ergo.scf_params, (n))
#define var_get_stringOU(n)  var_get_string(ergo.output_params, (n))
#define var_get_stringXC(n)  var_get_string(ergo.XC_params, (n))

static void
var_print_tree(struct variable *tree, FILE *f, int indent)
{
  for(;tree; tree = tree->next) {
    for(int i=0; i<indent; i++) fputc(' ', f);

    switch(tree->type) {
    case VAR_STRING:
      fprintf(f, "STRING: %s = \"%s\"\n", tree->name,
              tree->v.str ? tree->v.str : "(empty)");
      break;
    case VAR_FLOAT:
      fprintf(f, "FLOAT : %s = %g\n", tree->name, tree->v.num);
      break;
    case VAR_INT:
      fprintf(f, "INT   : %s = %d\n", tree->name, tree->v.vint);
      break;
    case VAR_LIST:
      fprintf(f, "LIST  : %s\n", tree->name);
      var_print_tree(tree->v.list, f, indent+3);
    }
  }
}
/** starts processing the inline molecule input. Call to this routine
    should be followed by calls to es_add_atom and es_mol_commit.
    @param moleculeClass selects the main molecule (MOL_MAIN) or the ghost
    molecule (MOL_GHOST).
*/
void
es_mol_begin(enum MolType moleculeClass) {
  ergo.moleculeUnit = 1;
  ergo.readingMoleculeClass = moleculeClass;
  switch(moleculeClass) {
  case MOL_MAIN:  ergo.molecule.clear();      break;
  case MOL_GHOST: ergo.ghostMolecule.clear(); break;
  default: assert(0);
  }
  if(Basis_info) {
    delete Basis_info;
    Basis_info = NULL;
  }
}


/** adds single atom at given coordinates and given name. The charge
    is specified currently by the name of the element. */
void
es_add_atom(const char *name, double x, double y, double z)
{

  Molecule *m;
  switch(ergo.readingMoleculeClass) {
  case MOL_MAIN:  m = &ergo.molecule; break;
  case MOL_GHOST: m = &ergo.ghostMolecule; break;
  default: assert(0);
  }
  ergo_real charge = get_charge_int_from_atom_label(name);
  m->addAtom(charge,
	     x*ergo.moleculeUnit,
	     y*ergo.moleculeUnit,
	     z*ergo.moleculeUnit);
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "%s (charge=%3.1f) %f %f %f", name,
	    (double)charge,
            (double)(x*ergo.moleculeUnit), 
	    (double)(y*ergo.moleculeUnit), 
	    (double)(z*ergo.moleculeUnit));
}

/** Finish the inline definition of the molecule. */
void
es_mol_commit(void)
{
  printf("Inline %s molecule with %d atoms.\n", 
         ergo.readingMoleculeClass == MOL_GHOST ? "ghost" : "main",
         ergo.readingMoleculeClass  == MOL_GHOST
         ? ergo.ghostMolecule.getNoOfAtoms() : ergo.molecule.getNoOfAtoms());
}

/** Selects the units for the inline molecule format to be Angtroms,
    as opposed to default atomic units. */
void
es_mol_unit_angstrom(void)
{
  ergo.moleculeUnit = UNIT_one_Angstrom;
}

/** reads molecule data in the MOLECULE.INP (Dalton) or XYZ format.

    @param fname contains the file name to be opened and read. 

    @param moleculeClass determines whether it is the main molecule
    (MOL_MAIN) or the ghost molecule (MOL_GHOST) to be read.
*/
int
es_mol_read_molecule(const char *fname, enum MolType moleculeClass)
{
  char *basissetfile = NULL;
  Molecule *m;
  char **basisFileName;
  switch(moleculeClass) {
  case MOL_MAIN:
    m = &ergo.molecule;      basisFileName = &ergo.Basis;
    break;
  case MOL_GHOST:
    m = &ergo.ghostMolecule; basisFileName = &ergo.GhostBasis;
    break;
  default: assert(0);
  }
  int res = m->setFromMoleculeFile(fname, 
                                   0, /* we are guessing the net charge here */
                                   &basissetfile);
  if(basissetfile) {
    if(!*basisFileName) {
      *basisFileName = basissetfile;
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN,
                "Setting the Basis from the MOLECULE file to %s\n",
                *basisFileName);
    } else {
      ergo_free(basissetfile);
    }
  }
  if(Basis_info) {
    delete Basis_info;
    Basis_info = NULL;
  }
  return res;
}

int
es_set_nthreads(int nThreads)
{
  const char *thread_counters[] = {
    "J_K.threads_J", "J_K.threads_K",
    "mat.threads", "scf.no_of_threads_for_V"
  };
  for(unsigned i=0; i<sizeof(thread_counters)/sizeof(thread_counters[0]); i++) {
    struct variable *var = es_find_var(NULL, thread_counters[i]);
    if(!var) printf("Not found %s\n", thread_counters[i]);
    else es_assign_num(var, nThreads);
  }
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  // If not compiled with OpenMP, the matrix library can only use one
  // thread. Therefore, set mat.threads to 1 in this case.
  struct variable *var = es_find_var(NULL, "mat.threads");
  if(!var) printf("Not found %s\n", "mat.threads");
  else es_assign_num(var, 1);  
#endif
  dft_set_num_threads(nThreads);
  return 0;
}

int
es_set_nthreads_string(const char *str)
{
  bool doDetect = strcasecmp(str, "detect") == 0;
  if(strcasecmp(str, "env") == 0 || doDetect) {
      const char *env = getenv("OMP_NUM_THREADS");
      int defThreads = 1;
      if ( !(env && (defThreads=atoi(env)) > 0) ) {
	if(doDetect) {
	  FILE *f = fopen("/proc/cpuinfo", "rt");
	  if(f) {
	    char line[256];
	    defThreads = 0;
	    while(fgets(line, sizeof(line), f))
	      if(strncmp(line, "processor", 9) == 0)
		defThreads++;
	    fclose(f);
	    /* Protect against case when /proc/cpuinfo exits but
	       contains garbage. Unlikely but possible. */
	    if(defThreads == 0)
	      defThreads = 1;
	  }
	}
      } 
      es_set_nthreads(defThreads);
      return 0;
  } else
    return -1;
}

static void
jkparams_set_from_vars(JK::Params& jkp)
{
  memset(&jkp, 0, sizeof(jkp));
  jkp.threshold_J = var_get_realJK("threshold_2el_J");
  jkp.threshold_K = var_get_realJK("threshold_2el_K");
  jkp.multipole_threshold_factor = var_get_realJK("multipole_threshold_factor");
#ifndef SKIP_UNOFFICIAL_INPUT_PARAMS
  jkp.use_densfit_for_J          = var_get_intJK("use_densfit_for_J");
  jkp.use_differential_density   = var_get_intJK("use_differential_density");
#endif
  jkp.use_fmm                    = var_get_intJK("use_fmm");
  jkp.fmm_box_size               = var_get_realJK("fmm_box_size");
  jkp.fmm_no_of_branches = var_get_intJK("fmm_no_of_branches");
  jkp.fmm_branch_splitter_extent_1 = var_get_realJK("fmm_branch_splitter_extent_1");
  jkp.fmm_branch_splitter_extent_2 = var_get_realJK("fmm_branch_splitter_extent_2");
  jkp.fmm_branch_splitter_extent_3 = var_get_realJK("fmm_branch_splitter_extent_3");
  jkp.fmm_branch_splitter_extent_4 = var_get_realJK("fmm_branch_splitter_extent_4");
  jkp.fmm_branch_splitter_extent_5 = var_get_realJK("fmm_branch_splitter_extent_5");
  jkp.exchange_box_size         = var_get_realJK("exchange_box_size");
  jkp.noOfThreads_J             = var_get_intJK("threads_J");
  jkp.noOfThreads_K             = var_get_intJK("threads_K");
}

void
es_print_help()
{
  var_print_tree(ergo.var_list, stdout, 0);
  printf("\nAvailable commands:\n"
         "help\n"
         "molecule [ghost] \"FILENAME\"\n"
         "molecule_inline [Angstrom]\n"
         "ghost_inline [Angstrom]\n"
         "range NUM = START COUNT \"BASIS-SET\"\n"
         "run \"METHOD\", METHOD=HF or a DFT functional\n"
         "system \"CMD\"\n"
         "warranty\n"
         "quit\n"
         "get_excited_state \"METHOD\" NO_OF_STATES\n"
         "get_polarisability \"METHOD\" \"[XYZ]\" FREQUENCY\n"
         "get_polarisability \"METHOD\" all FREQUENCY\n"
         "set_nthreads(N) where N is a number\n"
         "set_nthreads(\"env\") uses OMP_NUM_THREADS to set the thread count\n"
         "set_nthreads(\"detect\") uses OMP_NUM_THREADS, or hardware info.\n");
}

void
es_print_help_var(const struct variable *var)
{
  printf("%s: %s\n", var->name, var->description);
}

/** Print precision that was selected for building the program. */
void
es_print_precision()
{
#ifdef PRECISION_SINGLE 
  const char *precision = "single";
#elif defined(PRECISION_LONG_DOUBLE)
  const char *precision = "long_double";
#else
  const char *precision = "double";
#endif
  puts(precision);
}

static int
es_rmdir_with_content(const char *dirname)
{
  DIR * dir = opendir(dirname);
  struct dirent *dp;
  if(!dir) return -1;

  std::list<std::string> filesToRemove;
  while ( (dp=readdir(dir)) ) {
    filesToRemove.push_front(dp->d_name);
  }
  closedir(dir);

  for(std::list<std::string>::const_iterator i=filesToRemove.begin();
      i != filesToRemove.end(); ++i) {
    std::string fname(dirname);
    fname.append(1, '/');
    fname.append(*i);
    if (unlink(fname.c_str()) != 0)
      return -1;     
  }
  return rmdir(dirname);
}
 
/** called when an actual calculation is to be commenced.  @param mode
    is the first specified keyword.  Some calculation types - like
    response ones - require the Fock matrix. @param save_pot tells
    whether saving it is required: The save_final_potential
    configuration parameter will be overriden if save_pot is true.
 */
int
es_run(const char *mode, int save_pot)
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "running %s\n", mode);

  if(!ergoIntegralInfo)
    ergoIntegralInfo = new IntegralInfo(true);
  
  if(var_get_int("enable_memory_usage_output"))
    enable_memory_usage_output();

  if(var_get_int("rand_seed")) {
    int seed = var_get_int("rand_seed");
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN,
	      "es_run: calling srand() with seed %9d.",
	      seed);
    srand(seed);
  }

  if (ergo.molecule.getNumberOfElectrons() <= 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
                "es_run: no electrons found. Number of atoms: %d",
                ergo.molecule.getNoOfAtoms());
      return -1;
  }

  if(!Basis_info) {
    int output_basis = var_get_int("output_basis");
    int use_6_d_funcs = var_get_int("use_6_d_functions");
    Basis_info = new BasisInfoStruct(use_6_d_funcs);
    const int do_basis_normalization = 1;
    const int skip_sort_shells = 0;

    if(ergo.Basis != NULL) {
      // We skip adding basis functions here if the special basis set string "none" is given.
      if(strcmp(ergo.Basis, "none") != 0) {
	if(Basis_info->addBasisfuncsForMolecule(ergo.molecule, 
						ergo.Basis,
						Ergo::NO_OF_BASIS_SET_RANGES,
						ergo.basissetRangeList,
						*ergoIntegralInfo, 
						output_basis,
						do_basis_normalization,
						skip_sort_shells) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		      "error in basisInfo->add_basisfuncs_for_molecule "
		      "for main basis set, Basis='%s'",
		      ergo.Basis);
	    return -1;
	  }
      }
    }

    if(ergo.ghostMolecule.getNoOfAtoms() > 0 &&
       Basis_info->addBasisfuncsForMolecule(ergo.ghostMolecule, 
                                            ergo.GhostBasis,
                                            Ergo::NO_OF_BASIS_SET_RANGES,
                                            ergo.basissetRangeListGhost,
                                            *ergoIntegralInfo, 
                                            output_basis,
                                            do_basis_normalization,
                                            skip_sort_shells) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		  "error in basisInfo->add_basisfuncs_for_molecule "
		  "for ghost basis set, Basis='%s'",
                  ergo.GhostBasis);
	return -1;
      }
  } /* else reuse basis info since none of the geometry, Basis
     * has changed. */

  if(Basis_info->noOfBasisFuncs<1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
                "Cannot proceed: No basis functions defined.");
      return -1;
    }
  SCF::Options& scf = ergo.scfOptions;

  scf.calculation_identifier = var_get_stringSCF("calculation_identifier");  

  scf.method_and_basis_set = std::string(mode) + "/" + std::string(ergo.Basis);


  scf.electric_field.v[0] = var_get_realSCF("electric_field_x");
  scf.electric_field.v[1] = var_get_realSCF("electric_field_y");
  scf.electric_field.v[2] = var_get_realSCF("electric_field_z");

  scf.sparse_threshold_for_S = var_get_realSCF("sparse_threshold_for_S");
  scf.sparse_threshold_for_Z = var_get_realSCF("sparse_threshold_for_Z");
  scf.convergence_threshold = var_get_realSCF("convergence_threshold");
  scf.step_length_start     = var_get_realSCF("step_length_start");
  scf.step_length_giveup    = var_get_realSCF("step_length_giveup");
  scf.error_maxabs_for_diis = var_get_realSCF("error_maxabs_for_diis");
  scf.starting_guess_disturbance  = var_get_realSCF("starting_guess_disturbance");
  scf.purification_subspace_err_limit = var_get_realSCF("purification_subspace_err_limit");
  scf.purification_eigvalue_err_limit = var_get_realSCF("purification_eigvalue_err_limit");
  scf.puri_eig_acc_factor_for_guess = var_get_realSCF("puri_eig_acc_factor_for_guess");
  scf.gap_expected_lower_bound = var_get_realSCF("gap_expected_lower_bound");
  scf.shift_using_prev_density_matrix = var_get_realSCF("shift_using_prev_density_matrix");
  scf.electronic_temperature = var_get_realSCF("electronic_temperature");
  scf.purification_truncation_norm = mat::getNormType( var_get_stringSCF("purification_truncation_norm") );  
  scf.break_on_energy_increase    = var_get_intSCF("break_on_energy_increase");
  scf.create_basis_func_coord_file = var_get_intSCF("create_basis_func_coord_file");
  scf.output_homo_and_lumo_eigenvectors = var_get_intSCF("output_homo_and_lumo_eigenvectors");
  scf.create_mtx_file_S           = var_get_intSCF("create_mtx_file_S");
  scf.create_mtx_file_H_core      = var_get_intSCF("create_mtx_file_H_core");
  scf.create_mtx_files_F          = var_get_intSCF("create_mtx_files_F");
  scf.create_mtx_files_D          = var_get_intSCF("create_mtx_files_D");
  scf.create_mtx_files_dipole     = var_get_intSCF("create_mtx_files_dipole");
  scf.create_2el_integral_m_file  = var_get_intSCF("create_2el_integral_m_file");
  scf.force_restricted            = var_get_intSCF("force_restricted");
  scf.force_unrestricted          = var_get_intSCF("force_unrestricted");
  scf.max_no_of_diis_matrices     = var_get_intSCF("max_no_of_diis_matrices");
  scf.max_restart_count           = var_get_intSCF("max_restart_count");
  scf.min_number_of_iterations    = var_get_intSCF("min_number_of_iterations");
  scf.max_number_of_iterations    = var_get_intSCF("max_number_of_iterations");
  scf.no_of_careful_first_scf_steps     = var_get_intSCF("no_of_careful_first_scf_steps");
  scf.do_report_density_diff      = var_get_intSCF("do_report_density_diff");
  scf.no_of_impr_req_for_diis     = var_get_intSCF("no_of_impr_req_for_diis");
  scf.no_of_threads_for_V         = var_get_intSCF("no_of_threads_for_V");
  scf.output_density_at_every_step = var_get_intSCF("output_density_at_every_step");
  scf.save_final_potential        = var_get_intSCF("save_final_potential")
    || save_pot;
  scf.use_diagonalization         = var_get_intSCF("use_diagonalization");
  scf.use_diag_on_error           = var_get_intSCF("use_diag_on_error");
  scf.use_diag_on_error_guess     = var_get_intSCF("use_diag_on_error_guess");
  scf.output_mulliken_pop         = var_get_intSCF("output_mulliken_pop");
  scf.output_density_images       = var_get_intSCF("output_density_images");
  scf.output_density_images_only  = var_get_intSCF("output_density_images_only");
  scf.output_density_images_boxwidth = var_get_realSCF("output_density_images_boxwidth");
  scf.compute_gradient_fixeddens  = var_get_intSCF("compute_gradient_fixeddens");
  scf.verify_gradient_fixeddens   = var_get_intSCF("verify_gradient_fixeddens");
#ifndef SKIP_UNOFFICIAL_INPUT_PARAMS
  scf.do_f_thresh_verification    = var_get_intSCF("do_f_thresh_verification");
  scf.do_comparison_to_simple_purification   = var_get_intSCF("do_comparison_to_simple_purification");
  scf.do_puri_mmul_tests          = var_get_intSCF("do_puri_mmul_tests");
  scf.output_statistics_mfiles    = var_get_intSCF("output_statistics_mfiles");
  scf.do_sparsity_investigation   = var_get_intSCF("do_sparsity_investigation");
  scf.do_sparsity_investigation_reppuri   = var_get_intSCF("do_sparsity_investigation_reppuri");
  scf.sparsity_plots_resolution_r = var_get_intSCF("sparsity_plots_resolution_r");
  scf.sparsity_plots_resolution_m = var_get_intSCF("sparsity_plots_resolution_m");
  scf.do_acc_scan_J               = var_get_intSCF("do_acc_scan_J");
  scf.do_acc_scan_K               = var_get_intSCF("do_acc_scan_K");
  scf.do_acc_scan_Vxc             = var_get_intSCF("do_acc_scan_Vxc");
  scf.scan_no_of_steps            = var_get_intSCF("scan_no_of_steps");
  scf.scan_start_thresh           = var_get_realSCF("scan_start_thresh");
  scf.scan_step_factor            = var_get_realSCF("scan_step_factor");
  scf.write_guess_density_only    = var_get_intSCF("write_guess_density_only");
  scf.compute_core_density        = var_get_intSCF("compute_core_density");
  scf.no_of_core_electrons        = var_get_intSCF("no_of_core_electrons");
  scf.skip_H_core                 = var_get_intSCF("skip_H_core");
  scf.purification_create_m_files = var_get_intSCF("purification_create_m_files");
  scf.purification_use_rand_perturbation_for_alleigsint = var_get_intSCF("purification_use_rand_perturbation_for_alleigsint");
  scf.use_stochastic_orbs               = var_get_intSCF("use_stochastic_orbs");
  scf.stochastic_orbs_no_of_vectors     = var_get_intSCF("stochastic_orbs_no_of_vectors");
  scf.stochastic_orbs_use_unit_vectors  = var_get_intSCF("stochastic_orbs_use_unit_vectors");
  scf.stochastic_orbs_rand_param        = var_get_intSCF("stochastic_orbs_rand_param");
#endif
  scf.write_overlap_matrix        = var_get_intSCF("write_overlap_matrix");
  scf.use_simple_dense_H_core     = var_get_intSCF("use_simple_dense_H_core");
  scf.use_diis_always             = var_get_intSCF("use_diis_always");
  scf.use_simple_starting_guess   = var_get_int("use_simple_starting_guess");
  scf.use_dft  = mode && strcmp(mode, "HF") != 0 
    ? (var_get_intXC("sparse_mode") ? 2 : 1) : 0;

  jkparams_set_from_vars(ergo.jkOptions);

  ED::Params& ed = ergo.edOptions;
  ed.max_time = var_get_realED("max_time");
  ed.timestep = var_get_realED("timestep");
  ed.dc_pulse_strength = var_get_realED("dc_pulse_strength");
  ed.dc_pulse_time = var_get_realED("dc_pulse_time");
  ed.ac_pulse_max = var_get_realED("ac_pulse_max");
  ed.ac_pulse_omega = var_get_realED("ac_pulse_omega");
  ed.field_type = var_get_stringED("field_type");

  ergo.molecule.setNetCharge(var_get_int("charge"));
  int alpha_beta_diff = var_get_int("spin_polarization");
  if( (ergo.molecule.getNumberOfElectrons()-alpha_beta_diff)%2) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
		"Incompatible values for charge and spin_polarization"
		" settings.");
      return -1;
    }
    

  SCF::OutputOptions outputOptions;

  ergo.matOptions.threshold_inch = 
    var_get_realMA("threshold_inch");
  ergo.matOptions.sparse_threshold = 
    var_get_realMA("sparse_threshold");
  ergo.matOptions.sparse_matrix_block_size =
    var_get_intMA("sparse_matrix_block_size");
  ergo.matOptions.sparse_matrix_block_factor_1 =
    var_get_intMA("sparse_matrix_block_factor_1");
  ergo.matOptions.sparse_matrix_block_factor_2 =
    var_get_intMA("sparse_matrix_block_factor_2");
  ergo.matOptions.sparse_matrix_block_factor_3 =
    var_get_intMA("sparse_matrix_block_factor_3");
  ergo.matOptions.threads = var_get_intMA("threads");
  ergo.matOptions.parallelLevel = var_get_intMA("parallelLevel");
  ergo.matOptions.no_of_buffers_per_allocator = var_get_intMA("no_of_buffers_per_allocator");
  ergo.matOptions.prepare(*Basis_info);

  const char *tmpdir = var_get_string(NULL, "tmpdir");
  std::string subdir(tmpdir);
  if(var_get_intMA("write_to_file")) {
    subdir.append("/ergo_");
    char buf[20];
    snprintf(buf, sizeof(buf), "%i", getpid());
    subdir.append(buf);
    const char *matrixDir = subdir.c_str();
    if(mkdir(matrixDir, 0777) != 0 && errno != EEXIST) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
                "Cannot create tmp directory %s: %s", matrixDir,
		strerror(errno));
      return -1;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN,
              "Using directory '%s' for matrix storage.",
              matrixDir);
    static bool initializedMatLib = false;
    if(!initializedMatLib) {
      mat::FileWritable::setPath(matrixDir);
      mat::FileWritable::activate();
      initializedMatLib = true;
    }
  }
  if(scf.use_dft) {
    if(dft_setfunc(mode) == 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in dft_setfunc");
      return -1;
    }
    dftreport_();
  }

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Calling grid_set_tmpdir with tmpdir = '%s'", tmpdir);
  grid_set_tmpdir(tmpdir);
  Dft::GridParams gss(var_get_realXC("radint"),
		      var_get_intXC("angmin"),
		      var_get_intXC("angint"),
		      var_get_realXC("box_size"),
		      var_get_intXC("force_cubic_boxes"),
		      var_get_realXC("hicu_max_error"),
		      var_get_realXC("hicu_box_size"),
		      var_get_realXC("hicu_start_box_size_debug"),
		      var_get_intXC("hicu_use_error_per_volume"),
		      var_get_intXC("hicu_do_double_checking"),
		      var_get_intXC("hicu_compare_to_refined"),
		      var_get_intXC("hicu_use_energy_criterion"),
		      var_get_intXC("hicu_use_energy_criterion_only"),
		      var_get_intXC("hicu_do_variation_checking"));
  const char *gridType = var_get_stringXC("type");
  if (strcasecmp(gridType, "HICU") == 0) {
    gss.gridType = Dft::GridParams::TYPE_HICU;
  } else if (strcasecmp(gridType, "GC2") == 0) {
    gss.radialGridScheme =  Dft::GridParams::GC2;
  } else if (strcasecmp(gridType, "Turbo") == 0) {
    gss.radialGridScheme =  Dft::GridParams::TURBO;
  } else if (strcasecmp(gridType, "LMG") == 0) {
    gss.radialGridScheme =  Dft::GridParams::LMG;
  } else {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "Unknown radial grid type '%s'", gridType);
    return -1;
  }
    

  /* end of permutation initialization */
  const char *initial_density_fname = var_get_string(NULL, "initial_density");
  if(initial_density_fname && !*initial_density_fname)
    initial_density_fname = NULL;

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "Running %s%s", mode,
            initial_density_fname ? " (restarted)" : "");

  ergo_real threshold_integrals_1el = var_get_realJK("threshold_1el");

  /* Set charges in extraChargesMolecule according to parameters "extra_charges_mol_charge_h", "extra_charges_mol_charge_o" etc. */
  ergo_real extra_charges_atom_charge_h = var_get_real("extra_charges_atom_charge_h");
  ergo_real extra_charges_atom_charge_o = var_get_real("extra_charges_atom_charge_o");
  for(int i = 0; i < ergo.extraChargesMolecule.getNoOfAtoms(); i++) {
    Atom atom = ergo.extraChargesMolecule.getAtom(i);
    if(atom.charge == 1)
      atom.charge = extra_charges_atom_charge_h;
    else if(atom.charge == 8)
      atom.charge = extra_charges_atom_charge_o;
    else {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		"Error processing extraChargesMolecule: only O and H atoms are supported. Found atom with charge %5.2f", atom.charge);
      return -1;
    }
    ergo.extraChargesMolecule.replaceAtom(i, atom);
  }

  BasisInfoStruct basisInfoDensFit(0); /* FIXME */

  int noOfElectrons = ergo.molecule.getNumberOfElectrons();
  try {
    if (noOfElectrons % 2 == 1 || 
	scf.force_unrestricted ||
	alpha_beta_diff != 0) {
      // unrestricted SCF
      SCF_unrestricted SCF(ergo.molecule, 
			   ergo.extraChargesMolecule,
			   *Basis_info, 
			   basisInfoDensFit, /* FIXME */
			   *ergoIntegralInfo,
			   initial_density_fname,
			   ergo.jkOptions,
			   gss,
			   scf,
                           ergo.matOptions,
			   threshold_integrals_1el,
			   alpha_beta_diff);
      SCF.do_SCF_iterations();
      // Optionally use results to perform CI calculation.
      int do_ci_after_scf = var_get_int("do_ci_after_scf");
      if(do_ci_after_scf == 1) {
	// Do CI
	do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Preparing full matrices needed by CI routine..");
	// Get stuff from SCF object.
	symmMatrix S_matrix;
	SCF.get_overlap_matrix(S_matrix);
	symmMatrix H_core;
	SCF.get_H_core_matrix(H_core);
	symmMatrix FockMatrix_a;
	symmMatrix FockMatrix_b;
	SCF.get_Fock_matrices(FockMatrix_a, FockMatrix_b);
	int noOfElectrons_a, noOfElectrons_b;
	SCF.get_no_of_electrons(noOfElectrons_a, noOfElectrons_b);
	ergo_real energy;
	ergo_real nuclearEnergy;
	SCF.get_energy(energy, nuclearEnergy);
	// Create full matrices needed by CI routine.
	int n = Basis_info->noOfBasisFuncs;
	std::vector<ergo_real> S(n*n);
	std::vector<ergo_real> F_a(n*n);
	std::vector<ergo_real> F_b(n*n);
	std::vector<ergo_real> H_1(n*n);
	S_matrix.fullMatrix(S, 
			    ergo.matOptions.inversePermutationHML,
			    ergo.matOptions.inversePermutationHML);
	H_core.fullMatrix(H_1, 
			  ergo.matOptions.inversePermutationHML,
			  ergo.matOptions.inversePermutationHML);
	FockMatrix_a.fullMatrix(F_a, 
				ergo.matOptions.inversePermutationHML,
				ergo.matOptions.inversePermutationHML);
	FockMatrix_b.fullMatrix(F_b, 
				ergo.matOptions.inversePermutationHML,
				ergo.matOptions.inversePermutationHML);
	// Use default CI options
	CI::Options ci_options;
	// Call CI routine.
	if(do_CI(*Basis_info,
		 *ergoIntegralInfo,
		 ci_options,
		 &S[0],
		 &H_1[0],
		 &F_a[0],
		 &F_b[0],
		 noOfElectrons_a,
		 noOfElectrons_b,
		 nuclearEnergy,
		 energy
		 ) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "Error in do_CI.");
	    throw "error in DO_CI";
	  }
	do_output(LOG_CAT_RESULTS, LOG_AREA_SCF, "CI routine finished OK.");
      }
    }
    else {
      // restricted SCF
      SCF_restricted SCF(ergo.molecule, 
			 ergo.extraChargesMolecule,
			 *Basis_info, 
			 basisInfoDensFit, /* FIXME */
			 *ergoIntegralInfo,
			 initial_density_fname,
			 ergo.jkOptions,
			 gss,
			 scf,
                         ergo.matOptions,
			 threshold_integrals_1el);
      SCF.do_SCF_iterations();
      // Optionally do electron dynamics calculation. */
      int do_electron_dynamics_after_scf = var_get_int("do_electron_dynamics_after_scf");
      if(do_electron_dynamics_after_scf == 1) {
	// Get stuff from SCF object.
	symmMatrix S_matrix;
	SCF.get_overlap_matrix(S_matrix);
	symmMatrix H_core;
	SCF.get_H_core_matrix(H_core);
	symmMatrix FockMatrix;
	SCF.get_Fock_matrix(FockMatrix);
	triangMatrix invCholFactor;
	SCF.get_invCholFactor_matrix(invCholFactor);
	symmMatrix densityMatrix;
	SCF.get_density_matrix(densityMatrix);
	JK::ExchWeights CAM_params;
	get_hf_weight_and_cam_params(scf.use_dft, &CAM_params.alpha,
				     &CAM_params.beta, &CAM_params.mu);
	CAM_params.computeRangeSeparatedExchange = CAM_params.beta != ergo_real(0.0);
	do_tdhf_dynamics(*Basis_info,
			 *ergoIntegralInfo,
			 ergo.molecule,
			 ergo.extraChargesMolecule,
			 ergo.matOptions,
			 CAM_params,
			 ergo.jkOptions,
			 FockMatrix,
			 densityMatrix,
			 S_matrix,
			 invCholFactor,
			 ed);
      }

    }
  } 
  catch (const std::bad_alloc & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
	      "\n"
	      "=============================================================\n"
	      "std::bad_alloc caught in es_run: '%s'",  e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
  }
  catch (const std::ios_base::failure & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
	      "\n"
	      "=============================================================\n"
	      "std::ios_base::failure caught in es_run: '%s'\n"
	      "Out of disk space?", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
  }  
  catch (const std::exception& e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Exception caught: %s\n", e.what());
    fprintf(stderr, "Exception caught: %s\n", e.what());
  } catch (const char* s) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Exception caught: %s\n", s);
    fprintf(stderr, "Exception caught: %s\n", s);
  }
  grid_free_files();
  if(var_get_intMA("write_to_file"))
    es_rmdir_with_content(subdir.c_str());
  return 0;
}

#if 0
static void
printmat(int n, const ergo_real *m, const char *name)
{
  printf("Printing matrix %s\n", name);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      printf("%10.5f", m[i + j*n]);
    puts("");
  }
}
#endif

/** ErgoE2Evaluator implements the linear tranformation of the trial
    vector/transition density matrix by the E[2] operator. The
    transition density matrix is supplied in @param dmat. The result
    is returned in @param fmat. */
class ErgoE2Evaluator : public LR::E2Evaluator {
  BasisInfoStruct *bi;
  Molecule *mol;
  bool use_xc;
public:
  ErgoE2Evaluator(BasisInfoStruct *bis, Molecule *m, const char *mode)
    : bi(bis), mol(m) {
    use_xc = mode && strcasecmp(mode, "HF") != 0;
  }
  virtual bool transform(const ergo_real *dmat, ergo_real *fmat) {
    JK::ExchWeights CAM_params;
    int nbast = bi->noOfBasisFuncs;
    get_hf_weight_and_cam_params(use_xc,
                                 &CAM_params.alpha, 
                                 &CAM_params.beta, 
                                 &CAM_params.mu);
    CAM_params.computeRangeSeparatedExchange = CAM_params.beta != ergo_real(0.0);

    //printmat(bi->noOfBasisFuncs, dmat, "transition density");
    memset(fmat, 0, nbast*nbast*sizeof(ergo_real));

    bool res = false;
#if 1
    jkparams_set_from_vars(ergo.jkOptions);
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "calling compute_J_by_boxes_nosymm");
    if(compute_J_by_boxes_nosymm(*bi,
				 *ergoIntegralInfo,
				 ergo.jkOptions,
				 fmat,
				 dmat) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_J_by_boxes_nosymm");
	return false;
      }

    if(CAM_params.alpha != 0.0 || CAM_params.computeRangeSeparatedExchange)
      {
	do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "calling compute_K_by_boxes");
	int n = bi->noOfBasisFuncs;
	ergo_real* K = new ergo_real[n*n];
	memset(K, 0, n*n*sizeof(ergo_real));
	int symmetryFlag = 0;

	if(compute_K_by_boxes(*bi,
			      *ergoIntegralInfo,
			      CAM_params,
			      ergo.jkOptions,
			      K,
			      NULL,
			      dmat,
			      NULL,
			      symmetryFlag) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_K_by_boxes");
	    return false;
	  }
	int i;
	ergo_real CAMhf_weight =
	  CAM_params.computeRangeSeparatedExchange ? 1.0 : CAM_params.alpha;
	for(i = 0; i < n*n; i++)
	  fmat[i] += CAMhf_weight * K[i];
	delete []K;
      }
    res = true;
#else
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "ergo_scripted.cc: calling compute_2e_matrix_simple\n");
    res =
      compute_2e_matrix_simple(bi, ergoIntegralInfo, hf_weight,
				   fmat, dmat) == 0;
#endif
    if(res && use_xc) {
      ergo_real       *dens_matrix = NULL;
      BasisInfoStruct *basis_read  = NULL;
      Dft::GridParams gss(var_get_realXC("radint"),
			  var_get_intXC("angmin"),
			  var_get_intXC("angint"),
			  var_get_realXC("box_size"),
			  var_get_intXC("force_cubic_boxes"),
			  var_get_realXC("hicu_max_error"),
			  var_get_realXC("hicu_box_size"),
			  var_get_realXC("hicu_start_box_size_debug"),
			  var_get_intXC("hicu_use_error_per_volume"),
			  var_get_intXC("hicu_do_double_checking"),
			  var_get_intXC("hicu_compare_to_refined"),
			  var_get_intXC("hicu_use_energy_criterion"),
			  var_get_intXC("hicu_use_energy_criterion_only"),
			  var_get_intXC("hicu_do_variation_checking"));

      if(ddf_load_density("density.bin", 1, *ergoIntegralInfo,
			  &basis_read, &dens_matrix)) {
	do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Cannot load last Fock matrix from potential.bin");
	return false;
      }
      dft_lin_resp_mt(*basis_read, *mol, gss, dens_matrix, dmat, fmat);
      ergo_free(dens_matrix);
      delete basis_read;
    }
    return res;
  }
};


class ErgoOperator : public LR::OneElOperator {
  int px, py, pz;
  public:
  ErgoOperator(int pow_x, int pow_y, int pow_z)
    : px(pow_x), py(pow_y), pz(pow_z){}
  void setDipoleOp(int pow_x, int pow_y, int pow_z) {
    px = pow_x; py = pow_y; pz = pow_z;
  }
  virtual void getOper(ergo_real *res) {
    compute_operator_matrix_full(*Basis_info, *Basis_info, px, py, pz, res);
  }
};

/** Computes the specified number of excited states. @param no_exc
    specifies number of the excited states to be computed, @param mode
    specifies the calculation type (HF, LDA, etc). */
int
es_getexc(const char *mode, int no_exc)
{
  if (no_exc<=0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "Number of excited states must be larger than 0\n");
    return 1;
  }
  if(es_run(mode, 1) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "es_run failed");
    return 2;
  }

  int nocc = ergo.molecule.getNumberOfElectrons();
  if(nocc%2 != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "I work only for the closed shell.\n");
    return 3;
  }
  nocc /= 2;
  int nbast = Basis_info->noOfBasisFuncs;
  
  /** FIXME: consider passing callback functions instead of entire
      matrices. The callback functions fill in specified blocks of
      data with overlap matrix and the Fock matrix.  Current solution
      keeps these two potentially huge data blocks allocated all the
      time in memory. */
  ergo_real       *fock_matrix = NULL;
  BasisInfoStruct *basis_read = NULL;
  if(ddf_load_density("potential.bin", 1, *ergoIntegralInfo,
                      &basis_read, &fock_matrix)) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "Cannot load last Fock matrix from potential.bin");
    return -1;
  } 
  //printmat(Basis_info->noOfBasisFuncs, fock_matrix, "FOCK");

  ergo_real *overlap_matrix = new ergo_real[nbast*nbast];
  if(compute_overlap_matrix(*Basis_info, *basis_read, overlap_matrix) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_overlap_matrix");
    ergo_free(fock_matrix);
    return -2;
  }
  try {
    const char *tmpdir = var_get_string(NULL, "tmpdir");
    ErgoE2Evaluator e2(basis_read, &ergo.molecule, mode);
    LR::EigenSolver solver(nbast, nocc, fock_matrix, overlap_matrix,
			   no_exc);
    solver.convThreshold = var_get_realLR("convergence_threshold");
    solver.increaseSubspaceLimit(var_get_intLR("max_iterations")*2);
    if(!solver.solve(e2, tmpdir && *tmpdir ) ) {
      printf("Not converged!\n");
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Not converged\n");
    }
    ErgoOperator dx(1,0,0);
    ErgoOperator dy(0,1,0);
    ErgoOperator dz(0,0,1);
    solver.computeMoments(dx, dy, dz);
    for(int i=0; i<no_exc; i++) {
      printf("Eigenvalue %2i: %15.9f Tran.Mom.: %15.9g\n", i+1,
             (double)solver.getFreq(i), (double)sqrt(solver.getTransitionMoment2(i)));
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
                "Eigenvalue %2i: %15.9f Tran.Mom.: %15.9g", i+1,
                (double)solver.getFreq(i),  (double)sqrt(solver.getTransitionMoment2(i)));
    }
  } catch(const char*s){
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Error encountered: %s\n", s);
  }
  ergo_free(fock_matrix);
  grid_free_files();
  delete []overlap_matrix;
  delete basis_read;
  return 0; /* success */
}

static const int*
getOperatorParams(int opname) {
  static const int OpX[] = { 1, 0, 0 };
  static const int OpY[] = { 0, 1, 0 };
  static const int OpZ[] = { 0, 0, 1 };
  switch( toupper(opname) ) {
  case 'X': return OpX;
  case 'Y': return OpY;
  case 'Z': return OpZ;
  default: 
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "polarisability: known operators are X, Y and Z.");
    return NULL;
  }
}

static void
solveForRHS(LR::SetOfEqSolver& solver, ErgoE2Evaluator& e2, int opName,
            const char *tmpdir, ergo_real freq)
{
  const int *opL, *opR;

  if(! (opR=getOperatorParams(opName)) ) 
    throw "Unknown operator name";

  ErgoOperator op(opR[0], opR[1], opR[2]);

  solver.setRHS(op);
  if(!solver.solve(e2, tmpdir && *tmpdir )) {
    printf("LR not converged!\n");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "LR not converged\n");
  }

  char opLs[] = "XYZ";
  for(int i = 0; opLs[i]; i++) {

    opL = getOperatorParams(opLs[i]);
    op.setDipoleOp(opL[0], opL[1], opL[2]);
    ergo_real polarisability = -solver.getPolarisability(op);

    printf("Response << %c | %c >> at freq %15.9f: %15.10g\n", opLs[i],
           opName, (double)freq, (double)polarisability);
    do_output(LOG_CAT_RESULTS, LOG_AREA_MAIN,
              "Response << %c | %c >> at freq %15.9f: %15.10g", opLs[i],
              opName, (double)freq, (double)polarisability);
  }
}

/** Computes a dynamical polarizability for an operator specified by
    the @param opName and frequency @param freq - please check what
    does the literature say about computing multiple operators and/or
    frequencies at the same time. Consider using enumerated constants
    for operators instead of arbitrary strings to enforce parameter
    checking. It can be too early in this place for that - the
    operator names should be checked down the execution pipeline.

    @param mode is the type of Hamiltonian (HF, or the xc functional).
    @param freq tells the frequency.
 */
int
es_get_polarisability(const char *mode, const char *opName, double freq) 
{
  int nocc = ergo.molecule.getNumberOfElectrons();

  if(nocc%2 != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "es_get_polarisability works only for a closed shell.\n");
    return 3;
  }
  nocc /= 2;

  if(opName && strlen(opName) != 1 && !getOperatorParams(opName[0]) ) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "polarisability: opname is a 1-character string eg. \"Z\".");
    return 4;
  }

  if(es_run(mode, 1) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "es_run failed");
    return 2;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN,
            "Polarisability calculation with FreQ: %g", freq);

  int nbast = Basis_info->noOfBasisFuncs;

  ergo_real       *fockMatrix = NULL;
  BasisInfoStruct *basis_read = NULL;
  if(ddf_load_density("potential.bin", 1, *ergoIntegralInfo,
                      &basis_read, &fockMatrix)) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "Cannot load last Fock matrix from potential.bin");
    return -1;
  } 
  //printmat(Basis_info->noOfBasisFuncs, fock_matrix, "FOCK");

  std::vector<ergo_real> overlapMatrix(nbast*nbast);
  if(compute_overlap_matrix(*Basis_info, *basis_read, &overlapMatrix[0]) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in compute_overlap_matrix");
    ergo_free(fockMatrix);
    return -2;
  }
  try {
    const char *tmpdir = var_get_string(NULL, "tmpdir");
    ErgoE2Evaluator e2(basis_read, &ergo.molecule, mode);

    LR::SetOfEqSolver solver(nbast, nocc, fockMatrix, &overlapMatrix[0],
                             freq);
    overlapMatrix.clear();
    delete fockMatrix;    fockMatrix = NULL;

    solver.convThreshold = var_get_realLR("convergence_threshold");
    solver.increaseSubspaceLimit(var_get_intLR("max_iterations")*2);
    
    if(opName) {
      for(int i=0; opName[i]; i++)
        solveForRHS(solver, e2, opName[i], tmpdir, freq);
    } else {
      for(const char *r="XYZ"; *r; r++) 
        solveForRHS(solver, e2, *r, tmpdir, freq);
    }

  } catch(const char*s) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
              "SetOfEqSolver encountered error: %s\n", s);
  }
  if(fockMatrix)    ergo_free(fockMatrix);
  delete basis_read;
  return 0; /* success */
}

/** initializes the input module by registering all the recognized
    variables, their types and default values. If configuration
    objects exist for some part of calculations, we make effort to
    take the default values they provide. */
void
Ergo::registerInputVariables()
{
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  const char *tmpdir = getenv("TMPDIR");

  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }

  if( !tmpdir ) {
    tmpdir = "/tmp";
  }
  int defThreadsOMP;
#ifdef _OPENMP
  defThreadsOMP = defThreads;
#else
  defThreadsOMP = 1;
#endif

#define KL ergo.J_K_params
#define KWJK(n,type,d) KW(ergo.J_K_params,n,type,jkOptions.n,(d))
  KWJK(use_fmm,                      VAR_INT,"Use multipole method for Coulomb matrix construction. This also enables linear scaling HF-exchange matrix construction.");
#ifndef SKIP_UNOFFICIAL_INPUT_PARAMS
  KWJK(use_differential_density,     VAR_INT,"Use \"differential density\" procedure to try to speed up Fock matrix construction.");
  KWJK(use_densfit_for_J,            VAR_INT,"Use \"density fitting\" for Coulomb matrix construction.");
#endif
  KW(KL,threshold_2el_J,             VAR_FLOAT,jkOptions.threshold_J,"Threshold value for Coulomb matrix construction.");
  KW(KL,threshold_2el_K,             VAR_FLOAT,jkOptions.threshold_K,"Threshold value for HF exchange matrix construction.");
  KW(KL,threshold_1el,               VAR_FLOAT, 1e-12,"Threshold value for one-electron (core) Hamiltonian matrix construction.");
  KW(KL,threads_K,                   VAR_INT,   defThreads,"Number of threads to use in Coulomb matrix construction.");
  KW(KL,threads_J,                   VAR_INT,   defThreads,"Number of threads to use in HF exchange matrix construction.");
  KWJK(multipole_threshold_factor,   VAR_FLOAT,"Factor to apply to threshold value in multipole part (far field) of Coulomb matrix construction.");
  KWJK(fmm_no_of_branches,           VAR_INT,"Number of branches to use in Coulomb matrix construction. If 0, default branch settings are used based on box size and max extent.");
  KWJK(fmm_branch_splitter_extent_5, VAR_FLOAT,"Extent splitter value 5 for branch division in Coulomb matrix construction.");
  KWJK(fmm_branch_splitter_extent_4, VAR_FLOAT,"Extent splitter value 4 for branch division in Coulomb matrix construction.");
  KWJK(fmm_branch_splitter_extent_3, VAR_FLOAT,"Extent splitter value 3 for branch division in Coulomb matrix construction.");
  KWJK(fmm_branch_splitter_extent_2, VAR_FLOAT,"Extent splitter value 2 for branch division in Coulomb matrix construction.");
  KWJK(fmm_branch_splitter_extent_1, VAR_FLOAT,"Extent splitter value 1 for branch division in Coulomb matrix construction.");
  KWJK(fmm_box_size,                 VAR_FLOAT,"Smallest box size to use in Coulomb matrix construction.");
  KWJK(exchange_box_size,            VAR_FLOAT,"Smallest box size to use in HF exchange matrix construction.");
#undef KL

#define KL ergo.XC_params
  KW(KL,type,                         VAR_STRING, "Turbo","Type of the radial quadrature, one of Turbo and LMG.");
  KW(KL,sparse_mode,                  VAR_INT,    1,"Enable sparse mode for DFT exchange-correlation matrix computation.");
  KW(KL,radint,                       VAR_FLOAT,  5e-9,"Accuracy of the radial integration for atomic grids");
  KW(KL,force_cubic_boxes,            VAR_INT,    0,"Forces cubic boxes in grid generation as opposed to rectangular cuboid");
  KW(KL,box_size,                     VAR_FLOAT,  5.0,"Upper limit on grid box size");
  KW(KL,angmin,                       VAR_INT,    6,"Minimal order of pruned angular atomic grids");
  KW(KL,angint,                       VAR_INT,    29,"Default order of angular atomc grids");
  KW(KL,hicu_max_error,               VAR_FLOAT,  1e-7,"Threshold value (max error per box) for Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_box_size,                VAR_FLOAT,  1.5,"Box size for Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_start_box_size_debug,    VAR_FLOAT,  0.0,"Debug box size param for Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_use_error_per_volume,    VAR_INT,    0,"Use \"error-per-volume\" measure in Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_do_double_checking,      VAR_INT,    1,"Do \"double-checking\" of errors in Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_compare_to_refined,      VAR_INT,    0,"Compare to refined grid in Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_use_energy_criterion,    VAR_INT,    0,"Use energy criterion in Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_use_energy_criterion_only,    VAR_INT,    0,"Use only energy criterion in Hierarchical Cubature (HiCu) DFT grid generation.");
  KW(KL,hicu_do_variation_checking,   VAR_INT,    0,"Do variation checking in Hierarchical Cubature (HiCu) DFT grid generation.");
#undef KL

#define KL ergo.output_params
#undef KL

#define KL ergo.mat_params
#define KWMAT(n,type,d) KW(ergo.mat_params,n,type,matOptions.n,d)
  KW(KL,write_to_file,                VAR_INT, 0,"Write unused matrices to disk to save memory.");
  KWMAT(threshold_inch,               VAR_FLOAT,"Threshold value for removal of small matrix elements within the inverse Cholesky (inch) computation.");
  KW(KL,threads,                      VAR_INT, defThreadsOMP,"Number of threads to use in matrix library.");
  KWMAT(sparse_threshold,             VAR_FLOAT,"Threshold for sparse matrix truncation (removal of small matrix elements).");
  KWMAT(sparse_matrix_block_factor_3, VAR_INT,"Block size factor determining block size at fourth lowest level.");
  KWMAT(sparse_matrix_block_factor_2, VAR_INT,"Block size factor determining block size at third lowest level.");
  KWMAT(sparse_matrix_block_factor_1, VAR_INT,"Block size factor determining block size at second lowest level: sparse_matrix_block_size * sparse_matrix_block_factor_1.");
  KWMAT(sparse_matrix_block_size,     VAR_INT,"Lowest level submatrix block size.");
  KWMAT(parallelLevel,                VAR_INT,"Level in matrix hierarchy for parallelization.");
  KWMAT(no_of_buffers_per_allocator,  VAR_INT,"Number of buffers per allocator to use in matrix memory allocation manager. A large value means that large chunks of memory will be allocated at a time.");
#undef KWMAT
#undef KL

#define KL ergo.lr_params
  KW(KL,max_iterations,              VAR_INT,    10    ,"Maximum number of iterations for linear response calculations.");
  KW(KL,convergence_threshold,       VAR_FLOAT,  9e-4  ,"Convergence threshold for linear response calculations.");
#undef KL

#define KL ergo.ed_params
#define KWED(n,type,d) KW(ergo.ed_params,n,type,edOptions.n,d)
  KW(ergo.ed_params,electron_dynamics_field_type,   VAR_STRING, "none", "Field type to use in electron dynamics simulation: \"none\", \"dc-pulse\", \"ac-pulse\".");
  KWED(max_time, VAR_FLOAT, "How long time to run electron dynamics. Unit: a.u.");
  KWED(timestep, VAR_FLOAT, "Timestep to use in electron dynamics. Unit: a.u.");
  KWED(dc_pulse_strength, VAR_FLOAT, "Strength of DC pulse to use in electron dynamics. Unit: a.u.");
  KWED(dc_pulse_time, VAR_FLOAT, "Duration of DC pulse to use in electron dynamics. Unit: a.u.");
  KWED(ac_pulse_max, VAR_FLOAT, "Maximum field strength E_max for AC pulse to use in electron dynamics. Unit: a.u.");
  KWED(ac_pulse_omega, VAR_FLOAT, "Frequency omega for AC pulse to use in electron dynamics. Unit: a.u.");
#undef KWED
#undef KL

#define KL ergo.scf_params
#define KWSCF(n,type,d) KW(ergo.scf_params,n,type,scfOptions.n,d)
#ifndef SKIP_UNOFFICIAL_INPUT_PARAMS
  KWSCF(do_f_thresh_verification,     VAR_INT,"When truncating Fock matrix, verify that the error matrix norm is below the requested threshold.");
  KWSCF(do_comparison_to_simple_purification,     VAR_INT,"In each density matrix construction, do a comparison to a simple purification scheme.");
  KWSCF(do_puri_mmul_tests,           VAR_INT,"In each density matrix construction, do some extra purification test runs comparing different matrix-matrix multiplication and truncation approaches.");
  KWSCF(output_statistics_mfiles,     VAR_INT,"Output m-file with statistics (timings etc) each SCF iteration.");
  KWSCF(do_sparsity_investigation,    VAR_INT,"Create distance-vs-magnitude plot files.");
  KWSCF(do_sparsity_investigation_reppuri,    VAR_INT,"Perform purification 3 times using different norms, and create corresponding distance-vs-magnitude plot files.");
  KWSCF(sparsity_plots_resolution_r,  VAR_INT,"Distance resolution parameter for sparsity plots.");
  KWSCF(sparsity_plots_resolution_m,  VAR_INT,"Magnitude resolution parameter for sparsity plots.");
  KWSCF(use_stochastic_orbs,              VAR_INT,"Use stochastic orbitals as in the Stochastic DFT paper [PRL 111, 106402 (2013)].");
  KWSCF(stochastic_orbs_no_of_vectors,    VAR_INT,"Number of stochastic orbitals to use, if use_stochastic_orbs is turned on.");
  KWSCF(stochastic_orbs_use_unit_vectors, VAR_INT,"If use_stochastic_orbs is turned on, use unit vectors instead of random vectors.");
  KWSCF(stochastic_orbs_rand_param,       VAR_INT,"If use_stochastic_orbs is turned on, this value affects the way random numbers are generated. Setting this parameter to 2 means that a more elaborate, hopefully better, method is used to generate random numbers.");
  KWSCF(write_guess_density_only,     VAR_INT,"Only generate starting guess density, write it to file and exit.");
  KWSCF(compute_core_density,     VAR_INT,"Each time the density matrix is computed, also compute a 'core density matrix' using only the core electrons.");
  KWSCF(no_of_core_electrons,     VAR_INT,"If compute_core_density is set, use this number of core electrons when coputing the 'core density matrix'.");
  KWSCF(skip_H_core,                  VAR_INT,"Skip computation of 1-electron (core) Hamiltonian matrix. This gives bogus results, only useful if looking at results of initialization part or first cycle where H_core does not matter yet.");
  KWSCF(do_acc_scan_J,                VAR_INT,"Perform \"accuracy scan\" for Coulomb matrix.");
  KWSCF(do_acc_scan_K,                VAR_INT,"Perform \"accuracy scan\" for HF exchange matrix.");
  KWSCF(do_acc_scan_Vxc,              VAR_INT,"Perform \"accuracy scan\" for DFT exchange-correlation matrix.");
  KWSCF(scan_no_of_steps,             VAR_INT,"Number of steps to use in \"accuracy scans\" of J, K, Vxc matrices.");
  KWSCF(scan_start_thresh,            VAR_FLOAT,"Scan start threshold value to use for \"accuracy scans\" of J, K, Vxc matrices.");
  KWSCF(scan_step_factor,             VAR_FLOAT,"Step factor to use for \"accuracy scans\" of J, K, Vxc matrices.");
  KWSCF(purification_create_m_files,  VAR_INT,"Create m files with information about purification."),
  KWSCF(purification_use_rand_perturbation_for_alleigsint, VAR_INT,"Use random perturbation to attempt getting easier Lanczos convergence when determining min/max eigenvalues of Fock matrix before purification."),
#endif
  KWSCF(use_simple_dense_H_core,      VAR_INT,"Use simple dense matrix computation of 1-electron (core) Hamiltonian matrix.");
  KWSCF(use_diis_always,              VAR_INT,"Always use DIIS even if the energy goes up.");
  KWSCF(write_overlap_matrix,         VAR_INT,"Write overlap matrix to file.");
  KWSCF(use_diagonalization,          VAR_INT,"Use diagonalization instead of purification.");
  KWSCF(use_diag_on_error,            VAR_INT,"Use diagonalization if purification fails.");
  KWSCF(use_diag_on_error_guess,      VAR_INT,"Use diagonalization if purification fails during starting guess projection.");
  KWSCF(output_mulliken_pop,          VAR_INT,"Output Mulliken population analysis stuff (atomic charges and atomic spin densities) after SCF finished.");
  KWSCF(output_density_images,        VAR_INT,"Output density image files (including Gabedit gcube density files) after SCF finished.");
  KWSCF(output_density_images_only,   VAR_INT,"Output density image files (including Gabedit gcube density files) directly, without any SCF procedure, using the given starting guess density.");
  KWSCF(output_density_images_boxwidth, VAR_FLOAT, "Box width to use for density image files (including Gabedit gcube density files).");
  KWSCF(compute_gradient_fixeddens,   VAR_INT,"Compute gradient of energy with respect to nuclear positions for fixed electron density (Hellmann-Feynman forces).");
  KWSCF(verify_gradient_fixeddens,    VAR_INT,"Verify gradient of energy with respect to nuclear positions (Hellmann-Feynman forces), using finite differences. Used only when compute_gradient_fixeddens is set.");
  KWSCF(step_length_start,            VAR_FLOAT,"Start value of SCF step length parameter.");
  KWSCF(step_length_giveup,           VAR_FLOAT,"When step length in SCF procedure becomes smaller than this value, give up.");
  KWSCF(starting_guess_disturbance,   VAR_FLOAT,"Magnitude of random disturbance to apply to starting guess density matrix.");
  KWSCF(save_final_potential,         VAR_INT,"Save final Fock/Kohn-Sham matrix to file after SCF finished.");
  KWSCF(output_density_at_every_step, VAR_INT,"Write current density matrix to file in each iteration.");
  KW(KL,no_of_threads_for_V,          VAR_INT, defThreads,"Number of threads to use in V (1el electron-nuclear interaction) matrix construction.");
  KWSCF(no_of_impr_req_for_diis,      VAR_INT,"Number of consecutive energy improvements required before attempting to use DIIS.");
  KWSCF(min_number_of_iterations,     VAR_INT,"Minimum number of SCF iterations.");
  KWSCF(max_restart_count,            VAR_INT,"Max number of times to attempt to restart SCF procedure when it appears stuck.");
  KWSCF(max_number_of_iterations,     VAR_INT,"Maximum number of SCF iterations.");
  KWSCF(no_of_careful_first_scf_steps,VAR_INT,"Number of \"careful\" steps in beginning of SCF procedure. May be useful if starting guess is very bad.");
  KWSCF(do_report_density_diff,       VAR_INT,"In each SCF cycle, output norm of difference between current and previous density matrix.");
  KWSCF(max_no_of_diis_matrices,      VAR_INT,"Max number of DIIS matrices to use. Lower value may be useful to reduce disk/memory usage.");
  KWSCF(force_unrestricted,           VAR_INT,"Use unrestricted SCF even if numbers of alpha- and beta-electrons are equal.");
  KWSCF(force_restricted,             VAR_INT,"\"force_restricted\" parameter used for restricted open-shell calculations.");
  KWSCF(error_maxabs_for_diis,        VAR_FLOAT,"Do not use DIIS if error measure is larger than this value.");
  KWSCF(purification_subspace_err_limit,        VAR_FLOAT,"Requested accuracy in the occupied invariant subspace of the density matrix resulting from purification as measured by the sinus of the largest canonical angle between the exact and approximate subspaces.");
  KWSCF(purification_eigvalue_err_limit,        VAR_FLOAT,"Requested accuracy in eigenvalues of the density matrix resulting from purification.");
  KWSCF(puri_eig_acc_factor_for_guess, VAR_FLOAT,"When doing purification of starting guess density, use the normal eigvalue_err_limit multiplied by this factor.");
  KWSCF(gap_expected_lower_bound,        VAR_FLOAT,"Expected lower bound of HOMO-LUMO gap, used in early SCF iterations when accurate gap info is not yet available.");
  KWSCF(shift_using_prev_density_matrix, VAR_FLOAT, "As input to density matrix construction, Fock/KS matrix modified by subtracting previous density matrix times this factor. This is a way to artificially increase the HOMO-LUMO gap. A.k.a. level shifting.");
  KWSCF(electronic_temperature,          VAR_FLOAT,"Electronic temperature for Fermi-Dirac smearing in density matrix construction. Unit: a.u.");
  KW(ergo.scf_params,purification_truncation_norm,   VAR_STRING, "mixed","Matrix norm to be used for truncation of small matrix elements in purification, one of frob, mixed, and eucl.");
  KW(ergo.scf_params,electric_field_z,VAR_FLOAT, scfOptions.electric_field[2],"External electric field, z component.");
  KW(ergo.scf_params,electric_field_y,VAR_FLOAT, scfOptions.electric_field[1],"External electric field, y component.");
  KW(ergo.scf_params,electric_field_x,VAR_FLOAT, scfOptions.electric_field[0],"External electric field, x component.");
  KWSCF(sparse_threshold_for_S,       VAR_FLOAT,"Threshold value for truncation of overlap matrix S.");
  KWSCF(sparse_threshold_for_Z,       VAR_FLOAT,"Threshold value for truncation of inverse Cholesky factor Z.");
  KWSCF(convergence_threshold,        VAR_FLOAT,"Convergence threshold for SCF procedure; terminate SCF procedure when FDS-SDF error measure is below threshold.");
  KWSCF(break_on_energy_increase,     VAR_INT,"Break SCF procedure if energy increases.");
  KWSCF(create_basis_func_coord_file, VAR_INT,"Write text file with coordinates of basis functions.");
  KWSCF(output_homo_and_lumo_eigenvectors, VAR_INT,"Output HOMO and LUMO molecular orbital coefficients.");
  KWSCF(create_mtx_file_S,            VAR_INT,"Write overlap matrix to file in matrix market format (mtx).");
  KWSCF(create_mtx_file_H_core,       VAR_INT,"Write core Hamiltonian matrix to file in matrix market format (mtx).");
  KWSCF(create_mtx_files_F,           VAR_INT,"Write effective Hamiltonian matrices to file in matrix market format (mtx).");
  KWSCF(create_mtx_files_D,           VAR_INT,"Write density matrices to file in matrix market format (mtx).");
  KWSCF(create_mtx_files_dipole,      VAR_INT,"Write dipole matrices to file in matrix market format (mtx).");

  KW(ergo.scf_params,calculation_identifier,   VAR_STRING, "N/A","String to identify the calculation, used for example in mtx-files.");

  KWSCF(create_2el_integral_m_file,   VAR_INT,"Create m-file with rank-4 tensor containing the values of all 2-electron integrals. Very large; O(N^4); only use for small cases.");

#undef KWSCF
#undef KL

#define KL ergo.var_list
  KW(KL,use_simple_starting_guess,   VAR_INT,    1     ,"Use simple diagonal matrix starting guess.");
  KW(KL,tmpdir,                      VAR_STRING, tmpdir,"Directory for temporary files, usually something like \"/scratch\" or \"/tmp\".");
  KW(KL,spin_polarization,           VAR_INT,    0,"Spin polarization: difference between number of alpha- and beta-spin electrons; 0 for closed shell case.");
  KW(KL,output_basis,                VAR_INT,    0     ,"Write information about basis set to output file.");
  KW(KL,initial_density,             VAR_STRING, "","Filename of binary file from which starting guess density should be read.");
  KW(KL,ghost_basis,                 VAR_STRING, ""    ,"Basis set to use for \"ghost molecule\".");
  KW(KL,enable_memory_usage_output,  VAR_INT,    0,"Write information about memory usage to output file.");
  KW(KL,do_ci_after_scf,             VAR_INT,    0,"Perform a CI (configuration interaction) calculation after the SCF calculation.");
  KW(KL,do_electron_dynamics_after_scf,             VAR_INT,    0,"Perform electron dynamics (e.g. TDHF) calculation after the SCF calculation.");
  KW(KL,use_6_d_functions,           VAR_INT,    0,"If set to 1, use 6 basis functions in d-type shells instead of the usual 5 functions.");
  KW(KL,rand_seed,                   VAR_INT,    0,"Random seed initializing sequence of random numbers.");
  KW(KL,charge,                      VAR_INT,    0,"Net charge of molecule; this number (together with the nuclei) determines the number of electrons in the system.");
  KW(KL,extra_charges_atom_charge_o, VAR_FLOAT, -0.82, "Charge assigned to oxygen atoms in 'extra charges' molecule. The value -0.82 corresponds to the SPC model charge.");
  KW(KL,extra_charges_atom_charge_h, VAR_FLOAT, 0.41, "Charge assigned to hydrogen atoms in 'extra charges' molecule. The value 0.41 corresponds to the SPC model charge.");
  KW(KL,basis,                       VAR_STRING, ""    ,"Basis set file name.");
  KW(KL,scf,                         VAR_LIST,   ergo.scf_params,"List of input variables related to SCF procedure.");
  // ELIAS NOTE 2011-02-22: the "output" variable list is now empty,
  // and then it causes seg fault if included, so I commented it out.
  //KW(KL,output,                      VAR_LIST,   ergo.output_params,"List of input variables related to output.");
  KW(KL,mat,                         VAR_LIST,   ergo.mat_params,"List of input variables related to hierarchic matrix library.");
  KW(KL,lr,                          VAR_LIST,   ergo.lr_params,"List of input variables related to linear response calculations.");
  KW(KL,ed,                          VAR_LIST,   ergo.ed_params,"List of input variables related to electron dynamics calculations.");
  KW(KL,XC,                          VAR_LIST,   ergo.XC_params,"List of input variables related to the computation of the DFT exchange-correlation contribution.");
  KW(KL,J_K,                         VAR_LIST,   ergo.J_K_params,"List of input variables related to the computation of Coulomb and HF exchange matrices.");

#undef KL

}

static void
benchmark_mm()
{
  static const int SZ = 2000;
  static const double ALPHA=1.0, BETA=0.0;
  struct tms t1, t2;
  times(&t1);
  std::vector<double> a(SZ*SZ);
  std::vector<double> b(SZ*SZ);
  std::vector<double> c(SZ*SZ);
  for(int i=SZ*SZ-1; i>=0; i--) { a[i] = 2-i; b[i] = i;}
  mat::gemm("N","N", &SZ, &SZ, &SZ, &ALPHA, &a[0], &SZ,
            &b[0], &SZ, &BETA, &c[0], &SZ);
  times(&t2);
  printf("MM time: %6.3f s\nReference values:\nPhenom 2.4GHz/GOTO: 1.9s\n"
         "Phenom 2.4GHz/generic: 23.6s\n",
         double(t2.tms_utime-t1.tms_utime)/sysconf(_SC_CLK_TCK));
}

extern "C" int yyparse(void);
extern "C" void* yy_scan_string(const char *str);
extern "C" void* yy_create_buffer(FILE *f, int sz);
extern "C" void yy_switch_to_buffer(void *);
extern "C" void yy_delete_buffer(void *);
extern FILE *yyin; /**< file used by the lex-generated parser.*/
static void
ergo_parse_file(FILE *inputFile)
{
  void *buffer = yy_create_buffer(inputFile, 256);
  yy_switch_to_buffer(buffer);
  ergo_scanner_reading_stdin = isatty(fileno(inputFile));
  yyparse();
  yy_delete_buffer(buffer);
}
static void
ergo_parse_string(const char *str)
{
  void *buffer = yy_scan_string(str);
  yy_switch_to_buffer(buffer);
  yyparse();
  yy_delete_buffer(buffer);
}

void
es_warranty(void)
{
  puts(ERGO_LICENSE_TEXT_LONG "\n\n");
}

int
main(int argc, char *argv[])
{
  try {
    enable_output();
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, ERGO_LICENSE_TEXT_LONG);
  bool executed_something = false;
  FILE *inp_file = NULL;
  ergo.registerInputVariables();
  es_set_nthreads_string("detect");
  dft_init();
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "ERGO version %s\n", VERSION);

  for(int i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-e") == 0) {
      if(i+1<argc) {
        do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "parsing string: %s\n", argv[i+1]);
        ergo_parse_string(argv[i+1]);
        executed_something = true;
      } else {
        fprintf(stderr, "option -e encountered without argument\n");
        return 1;
      }
      i++;
    } else if(strcmp(argv[i], "-m") == 0) {
      if(i+1<argc) {
        do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Reading molecule file: %s\n", argv[i+1]);
	if(es_mol_read_molecule(argv[i+1], MOL_MAIN) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "es_mol_read_molecule() failed.");
	  return 1;
	}
      } else {
        fprintf(stderr, "option -m encountered without argument\n");
        return 1;
      }
      i++;
    } else if(strcmp(argv[i], "-c") == 0) {
      if(i+1<argc) {
        do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Reading extra charges from file: %s\n", argv[i+1]);
	if(readMoleculeFileInXyzFormat(ergo.extraChargesMolecule, 
				       argv[i+1], 
				       0, /* This param has no meaning in this case. */
				       false) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "readMoleculeFileInXyzFormat() for extra charges failed.");
	  return 1;
	}
      } else {
        fprintf(stderr, "option -c encountered without argument\n");
        return 1;
      }
      i++;
    } else if(strcmp(argv[i], "-g") == 0) {
      if(i+1<argc) {
        do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Reading ghost molecule from file: %s\n", argv[i+1]);
	if(readMoleculeFileInXyzFormat(ergo.ghostMolecule, 
				       argv[i+1], 
				       0, /* This param has no meaning in this case. */
				       false) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "readMoleculeFileInXyzFormat() for ghost molecule failed.");
	  return 1;
	}
      } else {
        fprintf(stderr, "option -g encountered without argument\n");
        return 1;
      }
      i++;
    } else if(strcmp(argv[i], "-h") == 0) {
      printf("Usage: ergo [args...]\n"
             "args can be: input file name\n"
             "             -e \"input line\"\n"
             "             -d VARIABLE_NAME variable to describe\n"
             "             -m molecule file name\n"
             "             -c extra-charges molecule file name\n"
             "             -h this message\n"
             "Arguments are interpreted in the order of encounter.\n"
             "Predefined variables:\n");
      es_print_help();
      return 0;
    } else if(strcmp(argv[i], "-d") == 0) {
      if(i+1<argc) {
        const struct variable *v = es_find_var(NULL, argv[i+1]);
        if (v) {
          es_print_help_var(v);
          return 0;
        } else {
          fprintf(stderr, "Variable %s not found.\n", argv[i+1]);
          return 1;
        }
      } else {
        fprintf(stderr, "option -d encountered without argument\n");
        return 1;
      }
    } else if(strcmp(argv[i], "-b") == 0) {
      benchmark_mm();
    } else {
      inp_file = fopen(argv[i], "rt");
      if (!inp_file) {
        fprintf(stderr, "Could not open '%s' for reading input.\n", argv[i]);
        return 1;
      }
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "reading input from %s\n", argv[i]);
      ergo_parse_file(inp_file);
      fclose(inp_file);
      executed_something = true;
    }
  }
  if(!executed_something) {
     if(isatty(fileno(stdin))) {
       printf(ERGO_LICENSE_TEXT_BRIEF "\n\n");
       printf("ERGO is ready!\n> ");
    }
    ergo_parse_file(stdin);
  }
  if(ergoIntegralInfo)
    delete ergoIntegralInfo;
  if(Basis_info)
    delete Basis_info;

  }
  catch (mat::AcceptableMaxIter & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "mat::AcceptableMaxIter caught in ergo main: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
  }
  catch (mat::Failure & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "mat::Failure caught in ergo main: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
  }
  catch (char const * e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "char* exception caught in ergo main: '%s'", e);
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
  }
  catch (std::bad_alloc & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "std::bad_alloc caught in ergo main: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
  }
  catch (std::ios_base::failure & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "std::ios_base::failure caught in ergo main: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Out of disk space?");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
  }
  catch (std::exception & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "std::exception caught in ergo main: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_MAIN, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "=============================================================");
  }

  return 0;
}
/* ===================================================================
   Replacement routines for some of the ERGO code.
*/
