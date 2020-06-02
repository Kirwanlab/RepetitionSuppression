#!/bin/bash

# Script to perform group-level analysis for "Repetition of Computer Security Warnings Results in Differential Repetition Suppression Effects as Revealed with Functional MRI"
# Data are available at https://openneuro.org/datasets/ds002363
# This performs a repeated measures multi-variate model (MVM) with factors for stimulus type (2 levels) and repetition (6 levels). 
#
# Based on scripts written by Nathan Muncy (https://github.com/nmuncy)


###--- Notes



### --- Set up --- ###										###??? update variables/arrays
#
# This is where the script will orient itself.
# Notes are supplied, and is the only section
# that really needs to be changed for each
# experiment.


# General variables
#parDir=~/fsl_groups/fslg_KirwanLab/compute/SittingBIDS
parDir=/Volumes/Yorick/RSE1_BIDS
workDir=${parDir}/derivatives								# par dir of data
outDir=${workDir}/grp-analysis-2020-01-14						# where output will be written (should match step3)
refFile=${workDir}/sub-02/run-1_rse1_scale+tlrc			# reference file, for finding dimensions etc

tempDir=/Volumes/Yorick/Templates/vold2_mni 							# desired template
priorDir=${tempDir}/priors_ACT								# location of atropos priors
mask=Intersection_GM_mask+tlrc								# this will be made, just specify name for the interesection gray matter mask


thr=0.3														# thresh value for Group_EPI_mask, ref Group_EPI_mean


### --- Create Masks --- ###
#
# This section will create a group mean intersection mask
# then threshold it at $thr to create a binary intersection mask.
# A gray matter mask will be constructed, and then the GM mask
# will be multiplied with the intersection mask to create a
# single GM intersection mask


cd $outDir


# intersection mask
if [ ! -f Group_epi_mask.nii.gz ]; then

    #this assumes that everybody in the derivatives folder is going into the mask
    list=`ls ${workDir}/sub*/mask_epi_anat+tlrc.HEAD`
    
	3dMean -prefix ${outDir}/Group_epi_mean.nii.gz $list
	3dmask_tool -input $list -frac $thr -prefix ${outDir}/Group_epi_mask.nii.gz
fi


# make $mask
if [ ! -f ${mask}.HEAD ]; then

	# GM mask
	c3d ${priorDir}/Prior2.nii.gz ${priorDir}/Prior4.nii.gz -add -o tmp_Prior_GM.nii.gz
	3dresample -master $refFile -rmode NN -input tmp_Prior_GM.nii.gz -prefix tmp_Template_GM_mask.nii.gz

	# combine GM and intersection mask
    c3d tmp_Template_GM_mask.nii.gz Group_epi_mask.nii.gz -multiply -o tmp_Intersection_GM_prob_mask.nii.gz
    c3d tmp_Intersection_GM_prob_mask.nii.gz -thresh 0.1 1 1 0 -o tmp_Intersection_GM_mask.nii.gz
    3dcopy tmp_Intersection_GM_mask.nii.gz $mask
	rm tmp*
fi

### --- MVM --- ###
#
# This will blur both stats and errts files according to $blurM, and
# then the blurred errts files will be used to model noise with
# an auto-correlation function. MVM scripts will be written and run,
# and none of this will happen on participants who move too much.
# A variable number of bx/wi subj variables is accepted, but this
# will not run a t-test.
#
# Currently, MVM post-hoc comparisons are permutations of bx/wi-subj
# variables. I.e. Behaviors A B C for groups Aut Con will yield
# comparisons of Aut-Con A-B, Aut-Con A-C, Aut-Con B-C. I could build
# more comparisons in the future.

blurInt=5
pref=rse1
print=ACF_raw_${pref}.txt
# outPre=${pref}_MVM_1
# outPre=${pref}_MVM_2
outPre=${pref}_MVM_4

if [ ! -s $print ]; then

    for k in sub-02 sub-03 sub-04 sub-05 sub-06 sub-07 sub-08 sub-09 sub-10 sub-11 sub-12 sub-13 sub-14 sub-15 sub-17 sub-18 sub-19 sub-20 sub-21 sub-22 sub-23 sub-24; do
        for m in rse1 rse1_errts; do
            hold=${workDir}/${k}/${m}
            file=${workDir}/${k}/rse1_errts_blur${blurInt}+tlrc
            # blur
            if [ ! -f ${hold}_blur${blurInt}+tlrc.HEAD ]; then
                3dmerge -prefix ${hold}_blur${blurInt} -1blur_fwhm $blurInt -doall ${hold}+tlrc
            fi
        done
        # parameter estimate
        3dFWHMx -mask $mask -input $file -acf >> $print
    done
fi

# simulate noise, determine thresholds
if [ ! -s ACF_MC_${pref}.txt ]; then

	sed '/ 0  0  0    0/d' $print > tmp

	xA=`awk '{ total += $1 } END { print total/NR }' tmp`
	xB=`awk '{ total += $2 } END { print total/NR }' tmp`
	xC=`awk '{ total += $3 } END { print total/NR }' tmp`

	3dClustSim -mask $mask -LOTS -iter 10000 -acf $xA $xB $xC > ACF_MC_${pref}.txt
	rm tmp
fi
		

# MVM Script

3dMVM  -prefix $outPre -jobs 12 -mask $mask \
    -wsVars "stimtype*repetition" \
    -num_glt 2 \
    -gltLabel 1 first_vs_second -gltCode 1 'repetition : 1*rep1 -1*rep2' \
    -gltLabel 2 linear_trend    -gltCode 2 'repetition : -5*rep1 -3*rep2 -1*rep3 1*rep4 3*rep5 5*rep6' \
    -dataTable  \
    Subj    stimtype    repetition    InputFile \
	sub-02 bus rep1 ../sub-02/rse1_4_blur5+tlrc[1] \
	sub-02 bus rep2 ../sub-02/rse1_4_blur5+tlrc[2] \
	sub-02 bus rep3 ../sub-02/rse1_4_blur5+tlrc[3] \
	sub-02 bus rep4 ../sub-02/rse1_4_blur5+tlrc[4] \
	sub-02 bus rep5 ../sub-02/rse1_4_blur5+tlrc[5] \
	sub-02 bus rep6 ../sub-02/rse1_4_blur5+tlrc[6] \
	sub-02 mal rep1 ../sub-02/rse1_4_blur5+tlrc[7] \
	sub-02 mal rep2 ../sub-02/rse1_4_blur5+tlrc[8] \
	sub-02 mal rep3 ../sub-02/rse1_4_blur5+tlrc[9] \
	sub-02 mal rep4 ../sub-02/rse1_4_blur5+tlrc[10] \
	sub-02 mal rep5 ../sub-02/rse1_4_blur5+tlrc[11] \
	sub-02 mal rep6 ../sub-02/rse1_4_blur5+tlrc[12] \
	sub-03 bus rep1 ../sub-03/rse1_4_blur5+tlrc[1] \
	sub-03 bus rep2 ../sub-03/rse1_4_blur5+tlrc[2] \
	sub-03 bus rep3 ../sub-03/rse1_4_blur5+tlrc[3] \
	sub-03 bus rep4 ../sub-03/rse1_4_blur5+tlrc[4] \
	sub-03 bus rep5 ../sub-03/rse1_4_blur5+tlrc[5] \
	sub-03 bus rep6 ../sub-03/rse1_4_blur5+tlrc[6] \
	sub-03 mal rep1 ../sub-03/rse1_4_blur5+tlrc[7] \
	sub-03 mal rep2 ../sub-03/rse1_4_blur5+tlrc[8] \
	sub-03 mal rep3 ../sub-03/rse1_4_blur5+tlrc[9] \
	sub-03 mal rep4 ../sub-03/rse1_4_blur5+tlrc[10] \
	sub-03 mal rep5 ../sub-03/rse1_4_blur5+tlrc[11] \
	sub-03 mal rep6 ../sub-03/rse1_4_blur5+tlrc[12] \
	sub-04 bus rep1 ../sub-04/rse1_4_blur5+tlrc[1] \
	sub-04 bus rep2 ../sub-04/rse1_4_blur5+tlrc[2] \
	sub-04 bus rep3 ../sub-04/rse1_4_blur5+tlrc[3] \
	sub-04 bus rep4 ../sub-04/rse1_4_blur5+tlrc[4] \
	sub-04 bus rep5 ../sub-04/rse1_4_blur5+tlrc[5] \
	sub-04 bus rep6 ../sub-04/rse1_4_blur5+tlrc[6] \
	sub-04 mal rep1 ../sub-04/rse1_4_blur5+tlrc[7] \
	sub-04 mal rep2 ../sub-04/rse1_4_blur5+tlrc[8] \
	sub-04 mal rep3 ../sub-04/rse1_4_blur5+tlrc[9] \
	sub-04 mal rep4 ../sub-04/rse1_4_blur5+tlrc[10] \
	sub-04 mal rep5 ../sub-04/rse1_4_blur5+tlrc[11] \
	sub-04 mal rep6 ../sub-04/rse1_4_blur5+tlrc[12] \
	sub-05 bus rep1 ../sub-05/rse1_4_blur5+tlrc[1] \
	sub-05 bus rep2 ../sub-05/rse1_4_blur5+tlrc[2] \
	sub-05 bus rep3 ../sub-05/rse1_4_blur5+tlrc[3] \
	sub-05 bus rep4 ../sub-05/rse1_4_blur5+tlrc[4] \
	sub-05 bus rep5 ../sub-05/rse1_4_blur5+tlrc[5] \
	sub-05 bus rep6 ../sub-05/rse1_4_blur5+tlrc[6] \
	sub-05 mal rep1 ../sub-05/rse1_4_blur5+tlrc[7] \
	sub-05 mal rep2 ../sub-05/rse1_4_blur5+tlrc[8] \
	sub-05 mal rep3 ../sub-05/rse1_4_blur5+tlrc[9] \
	sub-05 mal rep4 ../sub-05/rse1_4_blur5+tlrc[10] \
	sub-05 mal rep5 ../sub-05/rse1_4_blur5+tlrc[11] \
	sub-05 mal rep6 ../sub-05/rse1_4_blur5+tlrc[12] \
	sub-06 bus rep1 ../sub-06/rse1_4_blur5+tlrc[1] \
	sub-06 bus rep2 ../sub-06/rse1_4_blur5+tlrc[2] \
	sub-06 bus rep3 ../sub-06/rse1_4_blur5+tlrc[3] \
	sub-06 bus rep4 ../sub-06/rse1_4_blur5+tlrc[4] \
	sub-06 bus rep5 ../sub-06/rse1_4_blur5+tlrc[5] \
	sub-06 bus rep6 ../sub-06/rse1_4_blur5+tlrc[6] \
	sub-06 mal rep1 ../sub-06/rse1_4_blur5+tlrc[7] \
	sub-06 mal rep2 ../sub-06/rse1_4_blur5+tlrc[8] \
	sub-06 mal rep3 ../sub-06/rse1_4_blur5+tlrc[9] \
	sub-06 mal rep4 ../sub-06/rse1_4_blur5+tlrc[10] \
	sub-06 mal rep5 ../sub-06/rse1_4_blur5+tlrc[11] \
	sub-06 mal rep6 ../sub-06/rse1_4_blur5+tlrc[12] \
	sub-07 bus rep1 ../sub-07/rse1_4_blur5+tlrc[1] \
	sub-07 bus rep2 ../sub-07/rse1_4_blur5+tlrc[2] \
	sub-07 bus rep3 ../sub-07/rse1_4_blur5+tlrc[3] \
	sub-07 bus rep4 ../sub-07/rse1_4_blur5+tlrc[4] \
	sub-07 bus rep5 ../sub-07/rse1_4_blur5+tlrc[5] \
	sub-07 bus rep6 ../sub-07/rse1_4_blur5+tlrc[6] \
	sub-07 mal rep1 ../sub-07/rse1_4_blur5+tlrc[7] \
	sub-07 mal rep2 ../sub-07/rse1_4_blur5+tlrc[8] \
	sub-07 mal rep3 ../sub-07/rse1_4_blur5+tlrc[9] \
	sub-07 mal rep4 ../sub-07/rse1_4_blur5+tlrc[10] \
	sub-07 mal rep5 ../sub-07/rse1_4_blur5+tlrc[11] \
	sub-07 mal rep6 ../sub-07/rse1_4_blur5+tlrc[12] \
	sub-08 bus rep1 ../sub-08/rse1_4_blur5+tlrc[1] \
	sub-08 bus rep2 ../sub-08/rse1_4_blur5+tlrc[2] \
	sub-08 bus rep3 ../sub-08/rse1_4_blur5+tlrc[3] \
	sub-08 bus rep4 ../sub-08/rse1_4_blur5+tlrc[4] \
	sub-08 bus rep5 ../sub-08/rse1_4_blur5+tlrc[5] \
	sub-08 bus rep6 ../sub-08/rse1_4_blur5+tlrc[6] \
	sub-08 mal rep1 ../sub-08/rse1_4_blur5+tlrc[7] \
	sub-08 mal rep2 ../sub-08/rse1_4_blur5+tlrc[8] \
	sub-08 mal rep3 ../sub-08/rse1_4_blur5+tlrc[9] \
	sub-08 mal rep4 ../sub-08/rse1_4_blur5+tlrc[10] \
	sub-08 mal rep5 ../sub-08/rse1_4_blur5+tlrc[11] \
	sub-08 mal rep6 ../sub-08/rse1_4_blur5+tlrc[12] \
	sub-09 bus rep1 ../sub-09/rse1_4_blur5+tlrc[1] \
	sub-09 bus rep2 ../sub-09/rse1_4_blur5+tlrc[2] \
	sub-09 bus rep3 ../sub-09/rse1_4_blur5+tlrc[3] \
	sub-09 bus rep4 ../sub-09/rse1_4_blur5+tlrc[4] \
	sub-09 bus rep5 ../sub-09/rse1_4_blur5+tlrc[5] \
	sub-09 bus rep6 ../sub-09/rse1_4_blur5+tlrc[6] \
	sub-09 mal rep1 ../sub-09/rse1_4_blur5+tlrc[7] \
	sub-09 mal rep2 ../sub-09/rse1_4_blur5+tlrc[8] \
	sub-09 mal rep3 ../sub-09/rse1_4_blur5+tlrc[9] \
	sub-09 mal rep4 ../sub-09/rse1_4_blur5+tlrc[10] \
	sub-09 mal rep5 ../sub-09/rse1_4_blur5+tlrc[11] \
	sub-09 mal rep6 ../sub-09/rse1_4_blur5+tlrc[12] \
	sub-10 bus rep1 ../sub-10/rse1_4_blur5+tlrc[1] \
	sub-10 bus rep2 ../sub-10/rse1_4_blur5+tlrc[2] \
	sub-10 bus rep3 ../sub-10/rse1_4_blur5+tlrc[3] \
	sub-10 bus rep4 ../sub-10/rse1_4_blur5+tlrc[4] \
	sub-10 bus rep5 ../sub-10/rse1_4_blur5+tlrc[5] \
	sub-10 bus rep6 ../sub-10/rse1_4_blur5+tlrc[6] \
	sub-10 mal rep1 ../sub-10/rse1_4_blur5+tlrc[7] \
	sub-10 mal rep2 ../sub-10/rse1_4_blur5+tlrc[8] \
	sub-10 mal rep3 ../sub-10/rse1_4_blur5+tlrc[9] \
	sub-10 mal rep4 ../sub-10/rse1_4_blur5+tlrc[10] \
	sub-10 mal rep5 ../sub-10/rse1_4_blur5+tlrc[11] \
	sub-10 mal rep6 ../sub-10/rse1_4_blur5+tlrc[12] \
	sub-11 bus rep1 ../sub-11/rse1_4_blur5+tlrc[1] \
	sub-11 bus rep2 ../sub-11/rse1_4_blur5+tlrc[2] \
	sub-11 bus rep3 ../sub-11/rse1_4_blur5+tlrc[3] \
	sub-11 bus rep4 ../sub-11/rse1_4_blur5+tlrc[4] \
	sub-11 bus rep5 ../sub-11/rse1_4_blur5+tlrc[5] \
	sub-11 bus rep6 ../sub-11/rse1_4_blur5+tlrc[6] \
	sub-11 mal rep1 ../sub-11/rse1_4_blur5+tlrc[7] \
	sub-11 mal rep2 ../sub-11/rse1_4_blur5+tlrc[8] \
	sub-11 mal rep3 ../sub-11/rse1_4_blur5+tlrc[9] \
	sub-11 mal rep4 ../sub-11/rse1_4_blur5+tlrc[10] \
	sub-11 mal rep5 ../sub-11/rse1_4_blur5+tlrc[11] \
	sub-11 mal rep6 ../sub-11/rse1_4_blur5+tlrc[12] \
	sub-12 bus rep1 ../sub-12/rse1_4_blur5+tlrc[1] \
	sub-12 bus rep2 ../sub-12/rse1_4_blur5+tlrc[2] \
	sub-12 bus rep3 ../sub-12/rse1_4_blur5+tlrc[3] \
	sub-12 bus rep4 ../sub-12/rse1_4_blur5+tlrc[4] \
	sub-12 bus rep5 ../sub-12/rse1_4_blur5+tlrc[5] \
	sub-12 bus rep6 ../sub-12/rse1_4_blur5+tlrc[6] \
	sub-12 mal rep1 ../sub-12/rse1_4_blur5+tlrc[7] \
	sub-12 mal rep2 ../sub-12/rse1_4_blur5+tlrc[8] \
	sub-12 mal rep3 ../sub-12/rse1_4_blur5+tlrc[9] \
	sub-12 mal rep4 ../sub-12/rse1_4_blur5+tlrc[10] \
	sub-12 mal rep5 ../sub-12/rse1_4_blur5+tlrc[11] \
	sub-12 mal rep6 ../sub-12/rse1_4_blur5+tlrc[12] \
	sub-13 bus rep1 ../sub-13/rse1_4_blur5+tlrc[1] \
	sub-13 bus rep2 ../sub-13/rse1_4_blur5+tlrc[2] \
	sub-13 bus rep3 ../sub-13/rse1_4_blur5+tlrc[3] \
	sub-13 bus rep4 ../sub-13/rse1_4_blur5+tlrc[4] \
	sub-13 bus rep5 ../sub-13/rse1_4_blur5+tlrc[5] \
	sub-13 bus rep6 ../sub-13/rse1_4_blur5+tlrc[6] \
	sub-13 mal rep1 ../sub-13/rse1_4_blur5+tlrc[7] \
	sub-13 mal rep2 ../sub-13/rse1_4_blur5+tlrc[8] \
	sub-13 mal rep3 ../sub-13/rse1_4_blur5+tlrc[9] \
	sub-13 mal rep4 ../sub-13/rse1_4_blur5+tlrc[10] \
	sub-13 mal rep5 ../sub-13/rse1_4_blur5+tlrc[11] \
	sub-13 mal rep6 ../sub-13/rse1_4_blur5+tlrc[12] \
	sub-14 bus rep1 ../sub-14/rse1_4_blur5+tlrc[1] \
	sub-14 bus rep2 ../sub-14/rse1_4_blur5+tlrc[2] \
	sub-14 bus rep3 ../sub-14/rse1_4_blur5+tlrc[3] \
	sub-14 bus rep4 ../sub-14/rse1_4_blur5+tlrc[4] \
	sub-14 bus rep5 ../sub-14/rse1_4_blur5+tlrc[5] \
	sub-14 bus rep6 ../sub-14/rse1_4_blur5+tlrc[6] \
	sub-14 mal rep1 ../sub-14/rse1_4_blur5+tlrc[7] \
	sub-14 mal rep2 ../sub-14/rse1_4_blur5+tlrc[8] \
	sub-14 mal rep3 ../sub-14/rse1_4_blur5+tlrc[9] \
	sub-14 mal rep4 ../sub-14/rse1_4_blur5+tlrc[10] \
	sub-14 mal rep5 ../sub-14/rse1_4_blur5+tlrc[11] \
	sub-14 mal rep6 ../sub-14/rse1_4_blur5+tlrc[12] \
	sub-15 bus rep1 ../sub-15/rse1_4_blur5+tlrc[1] \
	sub-15 bus rep2 ../sub-15/rse1_4_blur5+tlrc[2] \
	sub-15 bus rep3 ../sub-15/rse1_4_blur5+tlrc[3] \
	sub-15 bus rep4 ../sub-15/rse1_4_blur5+tlrc[4] \
	sub-15 bus rep5 ../sub-15/rse1_4_blur5+tlrc[5] \
	sub-15 bus rep6 ../sub-15/rse1_4_blur5+tlrc[6] \
	sub-15 mal rep1 ../sub-15/rse1_4_blur5+tlrc[7] \
	sub-15 mal rep2 ../sub-15/rse1_4_blur5+tlrc[8] \
	sub-15 mal rep3 ../sub-15/rse1_4_blur5+tlrc[9] \
	sub-15 mal rep4 ../sub-15/rse1_4_blur5+tlrc[10] \
	sub-15 mal rep5 ../sub-15/rse1_4_blur5+tlrc[11] \
	sub-15 mal rep6 ../sub-15/rse1_4_blur5+tlrc[12] \
	sub-17 bus rep1 ../sub-17/rse1_4_blur5+tlrc[1] \
	sub-17 bus rep2 ../sub-17/rse1_4_blur5+tlrc[2] \
	sub-17 bus rep3 ../sub-17/rse1_4_blur5+tlrc[3] \
	sub-17 bus rep4 ../sub-17/rse1_4_blur5+tlrc[4] \
	sub-17 bus rep5 ../sub-17/rse1_4_blur5+tlrc[5] \
	sub-17 bus rep6 ../sub-17/rse1_4_blur5+tlrc[6] \
	sub-17 mal rep1 ../sub-17/rse1_4_blur5+tlrc[7] \
	sub-17 mal rep2 ../sub-17/rse1_4_blur5+tlrc[8] \
	sub-17 mal rep3 ../sub-17/rse1_4_blur5+tlrc[9] \
	sub-17 mal rep4 ../sub-17/rse1_4_blur5+tlrc[10] \
	sub-17 mal rep5 ../sub-17/rse1_4_blur5+tlrc[11] \
	sub-17 mal rep6 ../sub-17/rse1_4_blur5+tlrc[12] \
	sub-18 bus rep1 ../sub-18/rse1_4_blur5+tlrc[1] \
	sub-18 bus rep2 ../sub-18/rse1_4_blur5+tlrc[2] \
	sub-18 bus rep3 ../sub-18/rse1_4_blur5+tlrc[3] \
	sub-18 bus rep4 ../sub-18/rse1_4_blur5+tlrc[4] \
	sub-18 bus rep5 ../sub-18/rse1_4_blur5+tlrc[5] \
	sub-18 bus rep6 ../sub-18/rse1_4_blur5+tlrc[6] \
	sub-18 mal rep1 ../sub-18/rse1_4_blur5+tlrc[7] \
	sub-18 mal rep2 ../sub-18/rse1_4_blur5+tlrc[8] \
	sub-18 mal rep3 ../sub-18/rse1_4_blur5+tlrc[9] \
	sub-18 mal rep4 ../sub-18/rse1_4_blur5+tlrc[10] \
	sub-18 mal rep5 ../sub-18/rse1_4_blur5+tlrc[11] \
	sub-18 mal rep6 ../sub-18/rse1_4_blur5+tlrc[12] \
	sub-19 bus rep1 ../sub-19/rse1_4_blur5+tlrc[1] \
	sub-19 bus rep2 ../sub-19/rse1_4_blur5+tlrc[2] \
	sub-19 bus rep3 ../sub-19/rse1_4_blur5+tlrc[3] \
	sub-19 bus rep4 ../sub-19/rse1_4_blur5+tlrc[4] \
	sub-19 bus rep5 ../sub-19/rse1_4_blur5+tlrc[5] \
	sub-19 bus rep6 ../sub-19/rse1_4_blur5+tlrc[6] \
	sub-19 mal rep1 ../sub-19/rse1_4_blur5+tlrc[7] \
	sub-19 mal rep2 ../sub-19/rse1_4_blur5+tlrc[8] \
	sub-19 mal rep3 ../sub-19/rse1_4_blur5+tlrc[9] \
	sub-19 mal rep4 ../sub-19/rse1_4_blur5+tlrc[10] \
	sub-19 mal rep5 ../sub-19/rse1_4_blur5+tlrc[11] \
	sub-19 mal rep6 ../sub-19/rse1_4_blur5+tlrc[12] \
	sub-20 bus rep1 ../sub-20/rse1_4_blur5+tlrc[1] \
	sub-20 bus rep2 ../sub-20/rse1_4_blur5+tlrc[2] \
	sub-20 bus rep3 ../sub-20/rse1_4_blur5+tlrc[3] \
	sub-20 bus rep4 ../sub-20/rse1_4_blur5+tlrc[4] \
	sub-20 bus rep5 ../sub-20/rse1_4_blur5+tlrc[5] \
	sub-20 bus rep6 ../sub-20/rse1_4_blur5+tlrc[6] \
	sub-20 mal rep1 ../sub-20/rse1_4_blur5+tlrc[7] \
	sub-20 mal rep2 ../sub-20/rse1_4_blur5+tlrc[8] \
	sub-20 mal rep3 ../sub-20/rse1_4_blur5+tlrc[9] \
	sub-20 mal rep4 ../sub-20/rse1_4_blur5+tlrc[10] \
	sub-20 mal rep5 ../sub-20/rse1_4_blur5+tlrc[11] \
	sub-20 mal rep6 ../sub-20/rse1_4_blur5+tlrc[12] \
	sub-21 bus rep1 ../sub-21/rse1_4_blur5+tlrc[1] \
	sub-21 bus rep2 ../sub-21/rse1_4_blur5+tlrc[2] \
	sub-21 bus rep3 ../sub-21/rse1_4_blur5+tlrc[3] \
	sub-21 bus rep4 ../sub-21/rse1_4_blur5+tlrc[4] \
	sub-21 bus rep5 ../sub-21/rse1_4_blur5+tlrc[5] \
	sub-21 bus rep6 ../sub-21/rse1_4_blur5+tlrc[6] \
	sub-21 mal rep1 ../sub-21/rse1_4_blur5+tlrc[7] \
	sub-21 mal rep2 ../sub-21/rse1_4_blur5+tlrc[8] \
	sub-21 mal rep3 ../sub-21/rse1_4_blur5+tlrc[9] \
	sub-21 mal rep4 ../sub-21/rse1_4_blur5+tlrc[10] \
	sub-21 mal rep5 ../sub-21/rse1_4_blur5+tlrc[11] \
	sub-21 mal rep6 ../sub-21/rse1_4_blur5+tlrc[12] \
	sub-22 bus rep1 ../sub-22/rse1_4_blur5+tlrc[1] \
	sub-22 bus rep2 ../sub-22/rse1_4_blur5+tlrc[2] \
	sub-22 bus rep3 ../sub-22/rse1_4_blur5+tlrc[3] \
	sub-22 bus rep4 ../sub-22/rse1_4_blur5+tlrc[4] \
	sub-22 bus rep5 ../sub-22/rse1_4_blur5+tlrc[5] \
	sub-22 bus rep6 ../sub-22/rse1_4_blur5+tlrc[6] \
	sub-22 mal rep1 ../sub-22/rse1_4_blur5+tlrc[7] \
	sub-22 mal rep2 ../sub-22/rse1_4_blur5+tlrc[8] \
	sub-22 mal rep3 ../sub-22/rse1_4_blur5+tlrc[9] \
	sub-22 mal rep4 ../sub-22/rse1_4_blur5+tlrc[10] \
	sub-22 mal rep5 ../sub-22/rse1_4_blur5+tlrc[11] \
	sub-22 mal rep6 ../sub-22/rse1_4_blur5+tlrc[12] \
	sub-23 bus rep1 ../sub-23/rse1_4_blur5+tlrc[1] \
	sub-23 bus rep2 ../sub-23/rse1_4_blur5+tlrc[2] \
	sub-23 bus rep3 ../sub-23/rse1_4_blur5+tlrc[3] \
	sub-23 bus rep4 ../sub-23/rse1_4_blur5+tlrc[4] \
	sub-23 bus rep5 ../sub-23/rse1_4_blur5+tlrc[5] \
	sub-23 bus rep6 ../sub-23/rse1_4_blur5+tlrc[6] \
	sub-23 mal rep1 ../sub-23/rse1_4_blur5+tlrc[7] \
	sub-23 mal rep2 ../sub-23/rse1_4_blur5+tlrc[8] \
	sub-23 mal rep3 ../sub-23/rse1_4_blur5+tlrc[9] \
	sub-23 mal rep4 ../sub-23/rse1_4_blur5+tlrc[10] \
	sub-23 mal rep5 ../sub-23/rse1_4_blur5+tlrc[11] \
	sub-23 mal rep6 ../sub-23/rse1_4_blur5+tlrc[12] \
	sub-24 bus rep1 ../sub-24/rse1_4_blur5+tlrc[1] \
	sub-24 bus rep2 ../sub-24/rse1_4_blur5+tlrc[2] \
	sub-24 bus rep3 ../sub-24/rse1_4_blur5+tlrc[3] \
	sub-24 bus rep4 ../sub-24/rse1_4_blur5+tlrc[4] \
	sub-24 bus rep5 ../sub-24/rse1_4_blur5+tlrc[5] \
	sub-24 bus rep6 ../sub-24/rse1_4_blur5+tlrc[6] \
	sub-24 mal rep1 ../sub-24/rse1_4_blur5+tlrc[7] \
	sub-24 mal rep2 ../sub-24/rse1_4_blur5+tlrc[8] \
	sub-24 mal rep3 ../sub-24/rse1_4_blur5+tlrc[9] \
	sub-24 mal rep4 ../sub-24/rse1_4_blur5+tlrc[10] \
	sub-24 mal rep5 ../sub-24/rse1_4_blur5+tlrc[11] \
	sub-24 mal rep6 ../sub-24/rse1_4_blur5+tlrc[12]


