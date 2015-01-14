#!/bin/sh

    top_builddir=.
    top_srcdir=.


if test `"$top_builddir"/source/ergo -e precision` = single; then
    echo SKIPPED
    exit 0
fi

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

errorfilename=ergoscf.out.error.h2o_125_blyp

echo

guess_input='
basis = "STO-2G";
use_simple_starting_guess=1;
J_K.threshold_2el_J = 1e-6;
J_K.threshold_2el_K = 1e-6;
J_K.fmm_box_size = 10.0;
scf.convergence_threshold = 1e-2;
run "HF";
'

echo Testing h2o_125_1 BLYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input| "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_1.xyz > /dev/null
rm -f ergoscf.out
echo running BLYP/6-31Gss...
"$top_builddir"/source/ergo -m $top_srcdir/mol/h2o_125_1.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "BLYP"
EOINPUT
if 
check_final_energy -9549.3769105 2e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.1.blyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo Testing h2o_125_2 BLYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_2.xyz  > /dev/null
rm -f ergoscf.out
echo running BLYP/6-31Gss...
"$top_builddir"/source/ergo -m ../mol/h2o_125_2.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "BLYP"
EOINPUT
if 
check_final_energy -9549.3829921 2e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.2.blyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo Testing h2o_125_3 BLYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_3.xyz > /dev/null
rm -f ergoscf.out
echo running BLYP/6-31Gss...
"$top_builddir"/source/ergo -m ../mol/h2o_125_3.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "BLYP"
EOINPUT
if 
check_final_energy -9549.3737683 2e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.3.blyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo Testing h2o_125_4 BLYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_4.xyz > /dev/null
rm -f ergoscf.out
echo running BLYP/6-31Gss...
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_4.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "BLYP"
EOINPUT
if 
check_final_energy -9549.3769639 2e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.4.blyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



rm density.bin

echo
echo h2o_125 blyp tests completed successfully!
echo
