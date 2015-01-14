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

errorfilename=ergoscf.out.error.dftsparse

echo

echo Testing c2h8 LDA/6-31G*
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
C     0.0       0.0       0.0
C     0.0       0.0      10.0
H      0.000000     0.000000     1.084800
H      1.022759     0.000000    -0.361600
H     -0.511380     0.885735    -0.361600
H     -0.511380    -0.885735    -0.361600
H      0.000000     0.000000    11.084800
H      1.022759     0.000000     9.638400
H     -0.511380     0.885735     9.638400
H     -0.511380    -0.885735     9.638400
EOF
basis = "4-31G"
XC.sparse_mode = 1
run "LDA"
EOINPUT
if 
check_final_energy -80.097373423 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing h2o BLYP/6-31G**
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
XC.sparse_mode = 1
XC.radint=1e-13
XC.type="LMG"
run "BLYP"
EOINPUT
if 
check_final_energy -76.396247 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo
echo Sparse DFT tests completed successfully!
echo
