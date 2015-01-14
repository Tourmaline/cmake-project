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

errorfilename=ergoscf.out.error.h2o_1000

echo

echo Testing h2o_1000 HF/6-31G**
rm -f ergoscf.out
../source/ergo -m ../mol/h2o_1000.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
use_simple_starting_guess=1
J_K.fmm_box_size = 44
J_K.exchange_box_size = 44
scf.convergence_threshold = 1e-5
run "HF"
EOINPUT
if 
check_final_energy -76022.6431000 1e-4 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.h2o_1000
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm density.bin

echo
echo h2o_1000 Hartree-Fock tests completed successfully!
echo
