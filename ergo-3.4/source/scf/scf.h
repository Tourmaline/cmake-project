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

#ifndef SCFHEADER
#define SCFHEADER

#include <string.h>

#include "molecule.h"
#include "basisinfo.h"
#include "integrals_2el.h"
#include "matrix_typedefs.h"


namespace SCF {

static const int DISTURB_ELEMENT_MAX_COUNT = 60;

struct Options {
  std::string calculation_identifier;
  std::string method_and_basis_set;
  Vector3D electric_field;
  ergo_real electronic_temperature;
  ergo_real sparse_threshold_for_S;
  ergo_real sparse_threshold_for_Z;
  ergo_real convergence_threshold;
  ergo_real step_length_giveup;
  ergo_real step_length_start;
  ergo_real puri_eig_acc_factor_for_guess;
  ergo_real purification_conv_limit;
  ergo_real purification_eigvalue_err_limit;
  ergo_real purification_subspace_err_limit;  
  ergo_real gap_expected_lower_bound;
  mat::normType purification_truncation_norm;
  ergo_real subspace_factor_fock;
  ergo_real subspace_factor_dens;
  int use_artificial_subspace_disturbances;
  int no_of_threads_for_V;
  int purification_maxmul;
  int purification_create_m_files;
  int purification_ignore_failure;
  int purification_use_rand_perturbation_for_alleigsint;
  int use_dft;
  int use_simple_starting_guess;
  int use_diag_guess_from_file;
  int write_diag_dens_to_file;
  ergo_real starting_guess_disturbance;
  int sg_disturb_specific_elements;
  int disturbedElementIndexVector[DISTURB_ELEMENT_MAX_COUNT];
  ergo_real shift_using_prev_density_matrix;
  int skip_H_core;
  int use_simple_dense_H_core;
  int break_on_energy_increase;
  int force_restricted;  /**< use a restricted determinant for open shell. */
  int force_unrestricted; /**< use an unrestricted det. for closed shell. */
  int spin_flip_atom_count;
  int starting_guess_spin_diff;
  int max_no_of_diis_matrices;
  int max_restart_count;
  int no_of_impr_req_for_diis;
  int use_diis_always;
  int do_f_thresh_verification;
  int do_comparison_to_simple_purification;
  int do_puri_mmul_tests;
  int output_statistics_mfiles;
  int do_sparsity_investigation;
  int do_sparsity_investigation_reppuri;
  int sparsity_plots_resolution_r;
  int sparsity_plots_resolution_m;
  int no_of_careful_first_scf_steps;
  int do_report_density_diff;
  ergo_real error_maxabs_for_diis;
  int min_number_of_iterations;
  int max_number_of_iterations;
  int output_density_at_every_step;
  int output_csr_matrices_for_gao;
  int output_density_images;
  int output_density_images_only;
  int write_guess_density_only;
  int compute_core_density;
  int no_of_core_electrons;
  ergo_real output_density_images_boxwidth;
  int image_view_axis;
  int save_final_potential;
  int use_diagonalization;
  int use_diag_on_error;
  int use_diag_on_error_guess;
  int write_overlap_matrix;
  int save_full_matrices_for_matlab;
  int analyze_result_after_scf;
  int do_acc_scan_J;
  int do_acc_scan_K;
  int do_acc_scan_Vxc;
  int scan_do_invcholfactor_transf;
  int scan_no_of_steps;
  ergo_real scan_start_thresh;
  ergo_real scan_step_factor;
  int create_mtx_file_S;  
  int create_mtx_file_H_core;
  int create_mtx_files_F;  
  int create_mtx_files_D;  
  int create_mtx_files_dipole;
  int create_2el_integral_m_file;
  int create_basis_func_coord_file;
  int output_homo_and_lumo_eigenvectors;
  int output_mulliken_pop;
  int compute_gradient_fixeddens;
  int verify_gradient_fixeddens;

  int use_stochastic_orbs;
  int stochastic_orbs_no_of_vectors;
  int stochastic_orbs_use_unit_vectors;
  int stochastic_orbs_rand_param;

  /** Initializes all the fields to sane values. */
  Options() : calculation_identifier("N/A"),
       method_and_basis_set("N/A"),
       electric_field(0,0,0),
       electronic_temperature(0),
       sparse_threshold_for_S(1e-9),
       sparse_threshold_for_Z(1e-8),
       convergence_threshold(2e-7),
       step_length_giveup(0.00005),
       step_length_start(0.4),
       puri_eig_acc_factor_for_guess(1e-2),
       purification_conv_limit(0.1),
       purification_eigvalue_err_limit(1e-8),
       purification_subspace_err_limit(1e-6),
       gap_expected_lower_bound(0.05),
       purification_truncation_norm(mat::euclNorm),
       subspace_factor_fock(0.1),
       subspace_factor_dens(0.1),
       use_artificial_subspace_disturbances(0),
       no_of_threads_for_V(1),
       purification_maxmul(100),
       purification_create_m_files(0),
       purification_ignore_failure(0),
       purification_use_rand_perturbation_for_alleigsint(0),
       use_dft(0),
       use_simple_starting_guess(0),
       use_diag_guess_from_file(0),
       write_diag_dens_to_file(0),
       starting_guess_disturbance(0.0),
       sg_disturb_specific_elements(0),
       shift_using_prev_density_matrix(0.0),
       skip_H_core(0),
       use_simple_dense_H_core(0),
       break_on_energy_increase(0),
       force_restricted(0),
       force_unrestricted(0),
       spin_flip_atom_count(0),
       starting_guess_spin_diff(0),
       max_no_of_diis_matrices(10),
       max_restart_count(2),
       no_of_impr_req_for_diis(4),
       use_diis_always(0),
       do_f_thresh_verification(0),
       do_comparison_to_simple_purification(0),
       do_puri_mmul_tests(0),
       output_statistics_mfiles(0),
       do_sparsity_investigation(0),
       do_sparsity_investigation_reppuri(0),
       sparsity_plots_resolution_r(100),
       sparsity_plots_resolution_m(100),
       no_of_careful_first_scf_steps(0),
       do_report_density_diff(1),
       error_maxabs_for_diis(0.5),
       min_number_of_iterations(),
       max_number_of_iterations(),
       output_density_at_every_step(1),
       output_csr_matrices_for_gao(0),
       output_density_images(0),
       output_density_images_only(0),
       write_guess_density_only(0),
       compute_core_density(0),
       no_of_core_electrons(0),
       output_density_images_boxwidth(0.5),
       image_view_axis(),
       save_final_potential(0),
       use_diagonalization(0),
       use_diag_on_error(1),
       use_diag_on_error_guess(1),
       write_overlap_matrix(0),
       save_full_matrices_for_matlab(0),
       analyze_result_after_scf(0),
       do_acc_scan_J(0),
       do_acc_scan_K(0),
       do_acc_scan_Vxc(0),
       scan_do_invcholfactor_transf(1),
       scan_no_of_steps(16),
       scan_start_thresh(1e-9),
       scan_step_factor(sqrt((ergo_real)10)),
       create_mtx_file_S(0),
       create_mtx_file_H_core(0),
       create_mtx_files_F(0),
       create_mtx_files_D(0),
       create_mtx_files_dipole(0),
       create_2el_integral_m_file(0),
       create_basis_func_coord_file(0),
       output_homo_and_lumo_eigenvectors(0),
       output_mulliken_pop(0),
       compute_gradient_fixeddens(0),
       verify_gradient_fixeddens(0),
       use_stochastic_orbs(0),
       stochastic_orbs_no_of_vectors(100),
       stochastic_orbs_use_unit_vectors(0),
       stochastic_orbs_rand_param(1)
  { 
    memset(disturbedElementIndexVector, 0,
           sizeof(disturbedElementIndexVector));
  }
};

/** An object respresenting the configuration of the matrix
    library. All the thresholds and relevant parameters are collected
    in one object for the purposes of the input processing. */
struct MatOptions {
  mat::SizesAndBlocks size_block_info;
  std::vector<int> permutationHML;
  std::vector<int> inversePermutationHML;
  ergo_real sparse_threshold; /**< threshold value for sparse matrix
				 truncation. */
  ergo_real threshold_inch; /**< Truncation threshold in INCH function. */
  int sparse_matrix_block_size;
  int sparse_matrix_block_factor_3;
  int sparse_matrix_block_factor_2;
  int sparse_matrix_block_factor_1;
  int threads;
  int parallelLevel;
  int no_of_buffers_per_allocator;

  MatOptions() :
    sparse_threshold(1e-8),
    threshold_inch(0),
    sparse_matrix_block_size(32),
    sparse_matrix_block_factor_3(8),
    sparse_matrix_block_factor_2(8),
    sparse_matrix_block_factor_1(32),
    threads(1),
    parallelLevel(1),
    /* FIXME: there should be a param to set no_of_buffers_per_allocator, for large calculations it needs to be larger, e.g. 10 x larger seems to give much better performance of matrix operations for large cases. 
       This is also connected to blocksize, maybe the best solution would be to have a param determining the number of MegaBytes per allocator or something like that.  */
    no_of_buffers_per_allocator(20000) 
  {};
  ~MatOptions() {
  }
  /** after the parameters are called, this routine is to be called
      to figure out the basis set permutation. */
  void prepare(const BasisInfoStruct& basisInfo);
};

struct OutputOptions {
  OutputOptions() 
  {}
    
};

} /* end of SCF name space */



#endif
