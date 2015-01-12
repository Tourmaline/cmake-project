#!/bin/sh

if test "$top_builddir" = ""; then
    top_builddir=..
fi
if test "$top_srcdir" = ""; then
    top_srcdir=..
fi

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

errorfilename=ergoscf.out.error.h2o_125_b3lyp

guess_input='
basis = "STO-2G";
use_simple_starting_guess=1;
J_K.threshold_2el_J = 1e-6;
J_K.threshold_2el_K = 1e-6;
J_K.fmm_box_size = 10.0;
scf.convergence_threshold = 1e-2;
run "HF";
'

echo

echo Testing h2o_125_1 B3LYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_1.xyz > /dev/null
rm -f ergoscf.out
echo running B3LYP/6-31Gss...
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_1.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "B3LYP-G"
EOINPUT
if 
check_final_energy -9552.1494446 1e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.1.b3lyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo Testing h2o_125_2 B3LYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_2.xyz > /dev/null
rm -f ergoscf.out
echo running B3LYP/6-31Gss...
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_2.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "B3LYP-G"
EOINPUT
if 
check_final_energy -9552.1558339 1e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.2.b3lyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo Testing h2o_125_3 B3LYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_3.xyz > /dev/null
rm -f ergoscf.out
echo running B3LYP/6-31Gss...
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_3.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
scf.convergence_threshold = 1e-5
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "B3LYP-G"
EOINPUT
if 
check_final_energy -9552.1454731 1e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.3.b3lyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo Testing h2o_125_4 B3LYP/6-31Gss
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_125_4.xyz > /dev/null
rm -f ergoscf.out
echo running B3LYP/6-31Gss...
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
run "B3LYP-G"
EOINPUT
if 
check_final_energy -9552.1496239 1e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o125.4.b3lyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



rm density.bin

echo
echo h2o_125 b3lyp tests completed successfully!
echo
