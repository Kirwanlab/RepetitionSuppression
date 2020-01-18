#!/bin/bash

# Wrapper script for pre-processing script for "Repetition of Computer Security Warnings Results in Differential Repetition Suppression Effects as Revealed with Functional MRI"
# Data are available at https://openneuro.org/datasets/ds002363
#
# Based on scripts written by Nathan Muncy (https://github.com/nmuncy)


# stderr and stdout are written to ${outDir}/error_* and ${outDir}/output_* for troubleshooting.
# job submission output are time stamped for troubleshooting


workDir=~/compute/RSE1_BIDS   ###??? update this

scriptDir=${workDir}/code
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/TS1_${time}


mkdir -p $outDir

cd ${workDir}/derivatives
for i in sub*; do

    sbatch \
    -o ${outDir}/output_TS1_${i}.txt \
    -e ${outDir}/error_TS1_${i}.txt \
    ${scriptDir}/Task_step1_sbatch_preproc.sh $i

    sleep 1
done
