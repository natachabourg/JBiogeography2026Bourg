#!/bin/bash

#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=72:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
#conda activate opendrift
source activate opendrift
python OpenDrift_bluebottles.py
