#!/bin/sh

    top_builddir=.
    top_srcdir=.



# Prefer gawk - we know exactly what it can do.
# awk on Sun does not support functions, need to use nawk for this
if gawk '{print 1}'</dev/null > /dev/null 2>&1; then
   AWK=gawk
elif nawk '{print 1}'</dev/null > /dev/null 2>&1; then
   AWK=nawk
else
   AWK=awk
fi

. "$top_srcdir"/test/functions

errorfilename=ergoscf.out.error.lowacc

echo

echo Testing H UHF STO-3G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
H    0      0     0
EOF
basis = "STO-3G"
charge = 0
spin_polarization = 1
use_simple_starting_guess = 0
scf.use_diag_on_error = 0
scf.min_number_of_iterations = 2
J_K.use_fmm = 0
run "HF"
EOINPUT
if
check_final_energy -0.4665819 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


echo Testing H UHF 6-31G**
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
H    0      0     0
EOF
basis = "6-31Gss"
charge = 0
spin_polarization = 1
use_simple_starting_guess = 0
scf.use_diag_on_error = 0
scf.use_diag_on_error_guess = 0
scf.purification_eigvalue_err_limit = 1e-5
scf.purification_subspace_err_limit = 1e-3
J_K.use_fmm = 0
run "HF"
EOINPUT
if
check_final_energy -0.4982329 1e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


echo Testing H2[+] UHF STO-3G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
H    0      0     0
H    0    1.4     0
EOF
basis = "STO-3G"
charge = 1
spin_polarization = 1
use_simple_starting_guess = 0
scf.use_diag_on_error = 0
scf.purification_eigvalue_err_limit = 1e-6
scf.purification_subspace_err_limit = 1e-4
J_K.use_fmm = 0
run "HF"
EOINPUT
if
check_final_energy -0.5385113 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin



echo Testing h2o HF 6-31G**
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
scf.convergence_threshold = 2e-3
scf.purification_eigvalue_err_limit = 1e-4
scf.purification_subspace_err_limit = 1e-2
scf.use_diag_on_error = 0
run "HF"
EOINPUT
if
check_final_energy -76.0226431 1e-4 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin



echo Testing Be[-] UHF STO-2G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
Be   0    0   0
EOF
basis = "STO-2G"
charge = -1
spin_polarization = 1
scf.convergence_threshold = 1e-5
scf.purification_eigvalue_err_limit = 1e-5
scf.purification_subspace_err_limit = 1e-3
scf.use_diag_on_error = 1
scf.min_number_of_iterations = 2
run "HF"
EOINPUT
if
check_final_energy -13.6625792 1e-4 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


echo Testing Be[-] UHF 6-31G*
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
Be   0    0   0
EOF
basis = "6-31Gs"
charge = -1
scf.convergence_threshold = 1e-4
scf.purification_eigvalue_err_limit = 1e-4
scf.purification_subspace_err_limit = 1e-2
scf.use_diag_on_error = 1
spin_polarization = 1
run "HF"
EOINPUT
if
check_final_energy -14.4939039 1e-4 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin



echo
echo Low accuracy tests completed successfully!
echo
