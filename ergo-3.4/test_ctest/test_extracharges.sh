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

errorfilename=ergoscf.out.error.extracharges


echo

echo '1'                         > temp_molecule.xyz
echo ''                         >> temp_molecule.xyz
echo 'H        0.0   0.0   0.0' >> temp_molecule.xyz

echo '1'                         > temp_extra_charges.xyz
echo ''                         >> temp_extra_charges.xyz
echo 'H        0.0   0.0   0.5' >> temp_extra_charges.xyz

echo Testing H2 HF/cc-pVTZ with one of the H atoms given as "extra charge"
rm -f ergoscf.out
"$top_builddir"/source/ergo -m temp_molecule.xyz -c temp_extra_charges.xyz -g temp_extra_charges.xyz <<EOINPUT > /dev/null
basis = "cc-pVTZ"
ghost_basis = "cc-pVTZ"
charge = -1
extra_charges_atom_charge_h = 1
run "HF"
EOINPUT
if 
# Use reference energy corrected for missing nuclear_repulsion_energy
check_final_energy -2.12166118866 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo

echo '3'                                                     > temp_molecule.xyz
echo ''                                                     >> temp_molecule.xyz
echo 'O        0.457627840    -0.311951930     0.082447166' >> temp_molecule.xyz
echo 'H        0.457627840     0.594655950    -0.201387561' >> temp_molecule.xyz
echo 'H       -0.302917789    -0.472814751     0.628518998' >> temp_molecule.xyz

echo '3'                                                     > temp_extra_charges.xyz
echo ''                                                     >> temp_extra_charges.xyz
echo 'O       -0.463040796     0.126316403     2.021879464' >> temp_extra_charges.xyz
echo 'H       -0.463040796    -0.809137147     2.187489411' >> temp_extra_charges.xyz
echo 'H        0.381674058     0.386221704     1.673444396' >> temp_extra_charges.xyz

echo Testing two-h2o PBE0/3-21G with one of the molecules given as "extra charges"
rm -f ergoscf.out
"$top_builddir"/source/ergo -m temp_molecule.xyz -c temp_extra_charges.xyz -g temp_extra_charges.xyz <<EOINPUT > /dev/null
basis = "3-21G"
ghost_basis = "3-21G"
charge = -10
extra_charges_atom_charge_h = 1
extra_charges_atom_charge_o = 8
# make sure to use a grid that gives grid points where the density is; not only on atoms
XC.type="HICU"
XC.hicu_max_error = 1e-8
run "PBE0"
EOINPUT
if 
check_final_energy -186.35911209 5e-5 ; 
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
rm temp_molecule.xyz
rm temp_extra_charges.xyz

echo
echo extra-charges tests completed successfully!
echo
