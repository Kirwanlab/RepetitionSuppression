#!/bin/bash

# Script to blur the output of the first-level regression analysis for "Repetition of Computer Security Warnings Results in Differential Repetition Suppression Effects as Revealed with Functional MRI"
# Data are available at https://openneuro.org/datasets/ds002363
# Written to run locally.
#
# Based on scripts written by Nathan Muncy (https://github.com/nmuncy)

parDir=/Volumes/Yorick/RSE1_BIDS
derDir=${parDir}/derivatives
phase=rse1

cd $derDir

for i in sub*; do 
    cd $i

    #blur the functional dataset by 5mm
    3dmerge -prefix rse1_4_blur5 -1blur_fwhm 5 -doall rse1_4+tlrc

    cd $derDir
done

