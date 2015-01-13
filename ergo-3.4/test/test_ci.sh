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

errorfilename=ergoscf.out.error.ci

echo

echo Testing h2 FCI 6-31G**
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
H     0.0   0.0   0.0
H     0.0   0.0   1.4
EOF
basis = "6-31Gss"
scf.force_unrestricted = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if 
check_final_ci_corr_energy -0.03387 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing h2o FCI STO-3G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "STO-3G"
scf.force_unrestricted = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if 
check_final_ci_energy -75.0124258 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing h2o[+] FCI STO-3G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "STO-3G"
scf.force_unrestricted = 1
charge = 1
spin_polarization = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if 
check_final_ci_energy -74.6947713 1e-7 ; 
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
echo Configuration Interaction [CI] tests completed successfully!
echo
