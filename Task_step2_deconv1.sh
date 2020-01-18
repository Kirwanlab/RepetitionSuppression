#!/bin/bash


# First-level regression analysis script for "Repetition of Computer Security Warnings Results in Differential Repetition Suppression Effects as Revealed with Functional MRI"
# Data are available at https://openneuro.org/datasets/ds002363
# This script is written to run locally (i.e., not on a compute cluster).
#
# Based on scripts written by Nathan Muncy (https://github.com/nmuncy)


parDir=/Volumes/Yorick/RSE1_BIDS
derDir=${parDir}/derivatives
phase=rse1

cd $derDir

for i in sub*; do 
    cd $i

    runLength=()
    runLength=(${runLength[@]} `cat dfile.run-1_${phase}.1D | wc -l`)
    runLength=(${runLength[@]} `cat dfile.run-2_${phase}.1D | wc -l`)

    cat dfile.run-* > dfile_rall_rse1.1D

    1d_tool.py -infile dfile_rall_${phase}.1D -set_run_lengths ${runLength[@]} -demean -write motion_demean_${phase}.1D
    1d_tool.py -infile dfile_rall_${phase}.1D -set_run_lengths ${runLength[@]}  -derivative -demean -write motion_deriv_${phase}.1D
    1d_tool.py -infile motion_demean_${phase}.1D -set_run_lengths ${runLength[@]}  -split_into_pad_runs mot_demean_${phase}
    1d_tool.py -infile dfile_rall_${phase}.1D -set_run_lengths ${runLength[@]}  -show_censor_count -censor_prev_TR -censor_motion 0.3 motion_${phase}

    cat out.cen.run-*${phase}.1D > outcount_censor_${phase}.1D
    1deval -a motion_${phase}_censor.1D -b outcount_censor_${phase}.1D -expr "a*b" > censor_${phase}_combined.1D

#deconvolution 3, with (hopefully) correct timing files
# controls for stimulus size



#run the regression analysis
if [ ! -f ${phase}+tlrc.HEAD ]; then 
3dDeconvolve -input run-1_${phase}_scale+tlrc run-2_${phase}_scale+tlrc \
-polort A -float \
-num_stimts 25 \
-stim_file  1  mot_demean_${phase}.r01.1D'[0]' -stim_label 1  "Roll_1"  -stim_base 1 \
-stim_file  2  mot_demean_${phase}.r01.1D'[1]' -stim_label 2  "Pitch_1" -stim_base 2 \
-stim_file  3  mot_demean_${phase}.r01.1D'[2]' -stim_label 3  "Yaw_1"   -stim_base 3 \
-stim_file  4  mot_demean_${phase}.r01.1D'[3]' -stim_label 4  "dS_1"    -stim_base 4 \
-stim_file  5  mot_demean_${phase}.r01.1D'[4]' -stim_label 5  "dL_1"    -stim_base 5 \
-stim_file  6  mot_demean_${phase}.r01.1D'[5]' -stim_label 6  "dP_1"    -stim_base 6 \
-stim_file  7  mot_demean_${phase}.r02.1D'[0]' -stim_label 7  "Roll_2"  -stim_base 7 \
-stim_file  8  mot_demean_${phase}.r02.1D'[1]' -stim_label 8  "Pitch_2" -stim_base 8 \
-stim_file  9  mot_demean_${phase}.r02.1D'[2]' -stim_label 9  "Yaw_2"   -stim_base 9 \
-stim_file  10 mot_demean_${phase}.r02.1D'[3]' -stim_label 10 "dS_2"    -stim_base 10 \
-stim_file  11 mot_demean_${phase}.r02.1D'[4]' -stim_label 11 "dL_2"    -stim_base 11 \
-stim_file  12 mot_demean_${phase}.r02.1D'[5]' -stim_label 12 "dP_2"    -stim_base 12 \
-stim_file  13 timing_files/stimSizes.txt -stim_label 13 "stimSize" \
-stim_times 14 timing_files/bus1.txt 'BLOCK(3,1)' -stim_label 14   "bus1" \
-stim_times 15 timing_files/bus2.txt 'BLOCK(3,1)' -stim_label 15  "bus2" \
-stim_times 16 timing_files/bus3.txt 'BLOCK(3,1)' -stim_label 16  "bus3" \
-stim_times 17 timing_files/bus4.txt 'BLOCK(3,1)' -stim_label  17 "bus4" \
-stim_times 18 timing_files/bus5.txt 'BLOCK(3,1)' -stim_label  18 "bus5" \
-stim_times 19 timing_files/bus6.txt 'BLOCK(3,1)' -stim_label  19 "bus6" \
-stim_times 20 timing_files/mal1.txt 'BLOCK(3,1)' -stim_label  20 "mal1" \
-stim_times 21 timing_files/mal2.txt 'BLOCK(3,1)' -stim_label  21 "mal2" \
-stim_times 22 timing_files/mal3.txt 'BLOCK(3,1)' -stim_label  22 "mal3" \
-stim_times 23 timing_files/mal4.txt 'BLOCK(3,1)' -stim_label  23 "mal4" \
-stim_times 24 timing_files/mal5.txt 'BLOCK(3,1)' -stim_label  24 "mal5" \
-stim_times 25 timing_files/mal6.txt 'BLOCK(3,1)' -stim_label  25 "mal6" \
-noFDR -nofullf_atall \
-x1D X.${phase}.xmat.1D \
-xjpeg X.${phase}.jpg \
-x1D_uncensored X.${phase}.nocensor.xmat.1D \
-errts ${phase}_errts \
-bucket ${phase} \
-jobs 12 \
-GOFORIT 12
fi


cd $derDir

done

