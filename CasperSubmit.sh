#!/bin/bash -l
### Job Name
#PBS -N Physics-Impacts
### Charging account
#PBS -A P66770001
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=20GB
### Allow job to run up to 40 minutes
#PBS -l walltime=4:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe


# to check process:  qstat -u $USER

# module load python/2.7.14
ncar_pylib
./Physics-Impacts_CAPE-CIN-Clouds.py RUN


