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

errorfilename=ergoscf.out.error.lrhf

echo

echo Testing CNOF HF/STO-2G
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT 
basis = "STO-2G"
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
J_K.use_fmm = 0
scf.convergence_threshold = 1e-6
lr.convergence_threshold = 1e-5
get_excited_state "HF" 4
EOINPUT
if 
check_final_energy -252.028041454588873 1e-5 ; 
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 1 0.11364134 1e-5 ;
then
echo Eigenvalue 1 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 2 0.23956153 1e-5 ;
then
echo Eigenvalue 2 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 3 0.31319248 1e-5 ;
then
echo Eigenvalue 3 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if [ -n "$KEEP" ]; then
    mv ergoscf.out ergoscf.out_lr_exc_hf_ok
else
 rm ergoscf.out
fi

rm density.bin potential.bin

echo
echo TD-HF tests completed successfully!
echo
