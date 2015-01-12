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

errorfilename=ergoscf.out.error.ghost

echo


echo Testing h2o HF/6-31G** with extra basis funcs from ghost molecule
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
O        0.457627840    -0.311951930     0.082447166
H        0.457627840     0.594655950    -0.201387561
H       -0.302917789    -0.472814751     0.628518998
EOF
basis = "6-31Gss"
ghost_inline Angstrom
O       -0.463040796     0.126316403     2.021879464
H       -0.463040796    -0.809137147     2.187489411
H        0.381674058     0.386221704     1.673444396
EOF
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
if 
check_final_energy -76.0259735062 1e-7 ; 
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
echo Ghost-molecule tests completed successfully!
echo
