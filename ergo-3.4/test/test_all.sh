#! /bin/sh

test_list="
test_hf.sh
test_uhf.sh
test_hf_blsz1.sh 
test_dft_pure.sh 
test_dft_hybrid.sh
test_lr_exc.sh 
test_lr_exc_hf.sh 
test_fmm.sh 
test_fmm_b3lyp.sh
test_ext_elec_field.sh
test_mixedbasis.sh
test_ghost.sh
test_ci.sh 
"

s=0
for i in $test_list; do
    if ./$i; then 
        s=`echo $s + 1| bc -l`
    else
       echo "Test no $s -> $i failed."
       break
    fi
done

echo $s tests performed.
