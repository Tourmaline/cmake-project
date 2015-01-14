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

errorfilename=ergoscf.out.error.mixedbasis

echo

echo Testing three h2o with HF and mixed basis 6-31G 6-31G* 6-31G**
rm -f ergoscf.out

cat  <<EOM | "$top_builddir"/source/ergo
molecule_inline Angstrom
O        2.634829230    -0.987901199     1.894282198
H        2.634829230    -0.081293319     1.610447471
H        1.874283601    -1.148764019     2.440354030
O       -2.398712033    -0.928257960    79.638548870
H       -2.398712033    -1.863711511    79.804158818
H       -1.553997179    -0.668352659    79.290113803
O        0.126860479     0.972719769   154.235672931
H        0.126860479     1.052575258   155.182310719
H       -0.728382869     0.680994393   153.942490279
EOF
basis = "6-31Gss"
J_K.threshold_1el = 1e-13
J_K.threshold_2el_J = 1e-11
J_K.threshold_2el_K = 1e-11
use_simple_starting_guess = 0
range 1 = 0 3 "6-31G"
range 2 = 3 3 "6-31Gs"
run "HF"
EOM


if 
check_final_energy -228.016771 1e-6 ; 
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
echo Mixed basis set tests completed successfully!
echo
