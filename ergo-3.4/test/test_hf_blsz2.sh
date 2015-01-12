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

errorfilename=ergoscf.out.error.hfblsz1

echo

echo Testing h2o HF/6-31G** with matrix blocksize = 2
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
mat.sparse_matrix_block_size = 2
mat.sparse_matrix_block_factor_1 = 2
mat.sparse_matrix_block_factor_2 = 2
mat.sparse_matrix_block_factor_3 = 2
set_nthreads(1)
run "HF"
EOINPUT
if 
check_final_energy -76.0226431 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing cnof HF/6-31G with matrix blocksize = 2
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "6-31G"
mat.sparse_matrix_block_size = 2
mat.sparse_matrix_block_factor_1 = 2
mat.sparse_matrix_block_factor_2 = 2
mat.sparse_matrix_block_factor_3 = 2
set_nthreads(1)
run "HF"
EOINPUT
if 
check_final_energy -264.0542682 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing nh3[+] UHF/6-31G** with matrix blocksize = 2
rm -f ergoscf.out
"$top_builddir"/source/ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
N      0.000000     0.000000     0.000000
H      0.000000     0.000000     1.012316
H      0.969771     0.000000    -0.290392
H     -0.390071     0.887881    -0.290336
EOF
basis = "6-31Gss"
mat.sparse_matrix_block_size = 2
mat.sparse_matrix_block_factor_1 = 2
mat.sparse_matrix_block_factor_2 = 2
mat.sparse_matrix_block_factor_3 = 2
charge = 1
spin_polarization = 1
set_nthreads(1)
run "HF"
EOINPUT
if
check_final_energy -55.8546235 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.757013 1e-6 ;
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
echo Block size 2 tests completed successfully!
echo
