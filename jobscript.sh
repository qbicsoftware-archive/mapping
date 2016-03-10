#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -l nodes=1:ppn=10:cfc
#PBS -l walltime=200:00:00
#PBS -l mem=50g
# properties = {properties}

set -e

echo Running on machine $(hostname)

module load qbic/anaconda
module load bio/samtools/1.2
module load qbic/bwa


{exec_job}
exit 0
