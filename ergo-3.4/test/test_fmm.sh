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

errorfilename=ergoscf.out.error.fmm

echo

echo Testing four_h2o HF/6-31G using FMM
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
H                  8.675000    0.000000  -16.640000
H                  9.579936    0.000000  -17.920455
H                  7.203213    0.037757  -17.673336
H                  8.108149    0.037757  -18.953791
H                  9.241330    0.571501   -2.493601
H                 10.146266    0.571501   -3.774055
H                  7.134485    0.659897   -2.982496
H                  8.039421    0.659897   -4.262951
O                  8.675000    0.000000  -17.600000
O                  7.203213    0.037757  -18.633336
O                  9.241330    0.571501   -3.453601
O                  7.134485    0.659897   -3.942496
EOF
basis = "6-31G"
run "HF"
EOINPUT
if 
check_final_energy -303.7606521 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out


echo Testing four_h2o HF/6-31G* using FMM
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
H                  8.675000    0.000000  -16.640000
H                  9.579936    0.000000  -17.920455
H                  7.203213    0.037757  -17.673336
H                  8.108149    0.037757  -18.953791
H                  9.241330    0.571501   -2.493601
H                 10.146266    0.571501   -3.774055
H                  7.134485    0.659897   -2.982496
H                  8.039421    0.659897   -4.262951
O                  8.675000    0.000000  -17.600000
O                  7.203213    0.037757  -18.633336
O                  9.241330    0.571501   -3.453601
O                  7.134485    0.659897   -3.942496
EOF
basis = "6-31Gs"
initial_density = "density.bin"
run "HF"
EOINPUT
if 
check_final_energy -303.8352704 1e-7 ; 
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
echo FMM tests completed successfully!
echo
