#!/bin/sh


    top_builddir=.
    top_srcdir=.


if test `"$top_builddir"/source/ergo -e precision` = 'single'; then
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

errorfilename=ergoscf.out.error.lr_pol

echo

echo Testing CO SVWN5/4-31G polarizability
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT 
basis = "4-31G"
molecule_inline
N      0.00    -0.0   0.0
H      0.00    -1.9   0.1
H      1.64     0.9   0.0
H     -1.64     0.9   0.0
EOF
scf.convergence_threshold = 1e-6
J_K.use_fmm = 0
lr.convergence_threshold = 1e-5
XC.type="LMG"
XC.radint=1e-9
get_polarisability "SVWN5" all 0.0
initial_density="density.bin"
get_polarisability "SVWN5" all 0.2
EOINPUT
echo  # this is to get an extra newline
if 
check_final_energy -55.985633009 1e-5 ;
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

check_pol () {
if 
check_polarisability $1 $2 $3 $4 $5;
then
echo Polarisability $1 $2 at frequency $3 OK
else
echo ERROR Polarisability $1 $2 at frequency $3
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
}

check_pol X X   0.0   -7.877230629 1e-5
check_pol Y X   0.0    0.0 1e-5
check_pol Z X   0.0    0.0 1e-5
check_pol X Y   0.0    0.0 1e-5
check_pol Y Y   0.0   -7.936291744 1e-5
check_pol Z Y   0.0   0.2016103778 1e-5
check_pol X Z   0.0   0.0 1e-5
check_pol Y Z   0.0   0.2016104728 1e-5
check_pol Z Z   0.0    -2.45887058 1e-5
check_pol X X   0.2   -8.774076869 1e-5
check_pol Y X   0.2    0.0 1e-5
check_pol Z X   0.2    0.0 1e-5
check_pol X Y   0.2    0.0 1e-5
check_pol Y Y   0.2   -8.880466698  1e-5
check_pol Z Y   0.2    0.2234719254 1e-5
check_pol X Z   0.2    0.0 1e-5
check_pol Y Z   0.2    0.2234719254 1e-5
check_pol Z Z   0.2   -2.849382309 1e-5

rm density.bin potential.bin

echo Testing HF CAM-B3LYP/aug-cc-pVDZ polarizability
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT 
basis = "aug-cc-pVDZ"
molecule_inline Angstrom
F   0.0 0.0 0.0
H   0   0   0.924719
EOF
scf.convergence_threshold = 1e-6
J_K.use_fmm = 1
lr.convergence_threshold = 1e-5
XC.type="GC2"
XC.radint=1e-9
# FIXME: XC.sparse_mode= not implemented in response calculations.
XC.sparse_mode = 0
get_polarisability "camb3lyp" all 0.0
initial_density="density.bin"
get_polarisability "camb3lyp" all 0.1
EOINPUT
echo  # this is to get an extra newline


if 
check_final_energy -100.43540064008 1e-5 ;
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
check_pol X X   0.000000000    -4.458509918 1e-5
check_pol Y X   0.000000000  0.0  1e-5
check_pol Z X   0.000000000  0.0  1e-5
check_pol X Y   0.000000000  0.0  1e-5
check_pol Y Y   0.000000000    -4.458509918 1e-5
check_pol Z Y   0.000000000  0.0  1e-5
check_pol X Z   0.000000000  0.0  1e-5
check_pol Y Z   0.000000000  0.0  1e-5
check_pol Z Z   0.000000000    -6.444241063 4e-5
check_pol X X   0.100000000    -4.583001179 3e-5
check_pol Y X   0.100000000  0.0  1e-5
check_pol Z X   0.100000000  0.0  1e-5
check_pol X Y   0.100000000  0.0  1e-5
check_pol Y Y   0.100000000    -4.583001187 3e-5
check_pol Z Y   0.100000000  0.0  1e-5
check_pol X Z   0.100000000  0.0  1e-5
check_pol Y Z   0.100000000  0.0  1e-5
check_pol Z Z   0.100000000    -6.584672755 1e-5

rm density.bin potential.bin

echo
echo LR polarizability tests completed successfully!
echo

exit 0
