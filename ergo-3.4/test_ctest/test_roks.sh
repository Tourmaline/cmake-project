#!/bin/sh


    top_builddir=.
    top_srcdir=.


if test `"$top_builddir"/source/ergo -e precision` = 'single'; then
    echo SKIPPED
    exit 0
fi

. "$top_srcdir"/test/functions

errorfilename=ergoscf.out.error.roks

echo

echo Testing Li ROHF 3-21G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
Li    0      0     0
EOF
basis = "3-21G"
spin_polarization = 1
scf.min_number_of_iterations = 2
scf.force_restricted = 1
run "HF"
EOINPUT
if
check_final_energy -7.381510983 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out density.bin

echo Testing O2 SVWN5/6-31Gss
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O   0   0   0
O   0   0   2.2810
EOF
basis = "6-31Gss"
spin_polarization=2
scf.force_restricted = 1
XC.radint = 1e-9
XC.angint = 30
XC.angmin = 6
XC.type = "LMG"
run "LDA"
EOINPUT
if
check_final_energy -149.250773843 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out density.bin


echo Testing H2O[+] B3LYP/6-31Gss
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
O        -0.034119    0.000000   -0.287789
H         0.798604    0.000000    0.168673
H        -0.764486    0.000000    0.319117
EOF
basis = "6-31Gss"
charge = 1
spin_polarization = 1
scf.force_restricted = 1
#use_simple_starting_guess = 1
XC.type="GC2"
XC.radint = 1e-9
XC.angint = 30
XC.angmin = 6
scf.convergence_threshold = 1e-5
use_simple_starting_guess = 1
run "B3LYP"
EOINPUT
if
check_final_energy -75.860561961 4e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out density.bin




echo
echo Restricted Open-shell tests completed successfully!
echo
