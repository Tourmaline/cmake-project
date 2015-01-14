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

errorfilename=ergoscf.out.error.uhf

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
scf.min_number_of_iterations = 2
run "HF"
EOINPUT
if
check_final_energy -0.4665819 1e-7 ;
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
run "HF"
EOINPUT
if
check_final_energy -0.4982329 1e-7 ;
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
run "HF"
EOINPUT
if
check_final_energy -0.5385113 1e-7 ;
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


echo Testing H2 with bond distance 5 au, UHF/6-31G**
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
H   -2.5  0  0
H    2.5  0  0
EOF
basis = "6-31Gss"
scf.force_unrestricted = 1
 scf.starting_guess_disturbance = 0.01
run "HF"
EOINPUT
if
check_final_energy -0.9970912 1e-7 ;
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
use_simple_starting_guess = 0
scf.use_diag_on_error = 0
scf.min_number_of_iterations = 2
run "HF"
EOINPUT
if
check_final_energy -13.6625792 1e-5 ;
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
scf.use_diag_on_error = 0
scf.use_diag_on_error_guess = 1
spin_polarization = 1
use_simple_starting_guess = 0
run "HF"
EOINPUT
if
check_final_energy -14.4939039 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.752685 1e-6 ;
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
echo Unrestricted Hartree-Fock tests completed successfully!
echo
