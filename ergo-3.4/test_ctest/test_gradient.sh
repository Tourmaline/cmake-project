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

errorfilename1=ergoscf.out.error.gradient.1
errorfilename2=ergoscf.out.error.gradient.2
errorfilename3=ergoscf.out.error.gradient.3

epsilon=0.0001
tolerance=0.00001

echo

echo Testing gradient computation for h2o HF/6-31G** using finite differences and fixed basis set, achieved via ghost basis.


coord1=`echo 0.3 + $epsilon | bc` ; coord1=`printf "%9.6f" $coord1`
coord2=`echo 0.3 - $epsilon | bc` ; coord2=`printf "%9.6f" $coord2`
echo coord1 is $coord1
echo coord2 is $coord2

rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.gradient

"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     $coord1  0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav1

"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     $coord2  0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav2

if 
check_gradient_component_by_finite_diff 0 1 $epsilon $tolerance
then
echo OK
rm ergoscf.out.sav1 ergoscf.out.sav2 ergoscf.out.gradient
else
echo ERROR
mv ergoscf.out.sav1 $errorfilename1
mv ergoscf.out.sav2 $errorfilename2
mv ergoscf.out.gradient $errorfilename3
echo output files saved as $errorfilename1 $errorfilename2 $errorfilename3
exit 1
fi


coord1=`echo -1.8 + $epsilon | bc` ; coord1=`printf "%9.6f" $coord1`
coord2=`echo -1.8 - $epsilon | bc` ; coord2=`printf "%9.6f" $coord2`
echo coord1 is $coord1
echo coord2 is $coord2

rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.gradient

"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H     $coord1 0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav1

"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H     $coord2 0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav2

if 
check_gradient_component_by_finite_diff 1 0 $epsilon $tolerance
then
echo OK
rm ergoscf.out.sav1 ergoscf.out.sav2 ergoscf.out.gradient
else
echo ERROR
mv ergoscf.out.sav1 $errorfilename1
mv ergoscf.out.sav2 $errorfilename2
mv ergoscf.out.gradient $errorfilename3
echo output files saved as $errorfilename1 $errorfilename2 $errorfilename3
exit 1
fi


coord1=`echo 0.2 + $epsilon | bc` ; coord1=`printf "%9.6f" $coord1`
coord2=`echo 0.2 - $epsilon | bc` ; coord2=`printf "%9.6f" $coord2`
echo coord1 is $coord1
echo coord2 is $coord2

rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.gradient

"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      $coord1
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav1

"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      $coord2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav2

if 
check_gradient_component_by_finite_diff 2 2 $epsilon $tolerance
then
echo OK
rm ergoscf.out.sav1 ergoscf.out.sav2 ergoscf.out.gradient
else
echo ERROR
mv ergoscf.out.sav1 $errorfilename1
mv ergoscf.out.sav2 $errorfilename2
mv ergoscf.out.gradient $errorfilename3
echo output files saved as $errorfilename1 $errorfilename2 $errorfilename3
exit 1
fi


rm -f ergoscf.out
rm -f density.bin

echo
echo Gradient tests completed successfully!
echo
