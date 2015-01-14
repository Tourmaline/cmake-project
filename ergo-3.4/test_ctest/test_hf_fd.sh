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

errorfilename=ergoscf.out.error.hffd

echo


echo Testing h2o HF/Ahlrichs-VTZ with Fermi-Dirac smearing at T=20000K
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
O   0.000000000000000   0.000000000000000   0.000000000000000
H  -0.957281569814700   0.000000000000000   0.000000000000000
H   0.240007793647257   0.926706239896334   0.000000000000000
EOF
basis = "Ahlrichs-VTZ"
# Using conversion factor 3.166815e-06 from Kelvin to a.u.
scf.electronic_temperature = 0.0633363
scf.use_diagonalization = 1
scf.use_diis_always = 1
run "HF"
EOINPUT
if 
check_final_energy -76.00838544 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing cnof HF/6-31G with Fermi-Dirac smearing at T=20000K
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
C   0.000000000000000   0.000000000000000   0.000000000000000
N   0.687930370790000   0.105835441660000   0.211670883320000
O   0.793765812450000  -0.211670883320000   1.058354416600000
F  -0.476259487470000   1.058354416600000   0.264588604150000
EOF
basis = "6-31G"
# Using conversion factor 3.166815e-06 from Kelvin to a.u.
scf.electronic_temperature = 0.0633363
scf.use_diagonalization = 1
scf.use_diis_always = 1
run "HF"
EOINPUT
if 
check_final_energy -264.0295861799 1e-6 ; 
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
echo Hartree-Fock Fermi-Dirac T tests completed successfully!
echo
