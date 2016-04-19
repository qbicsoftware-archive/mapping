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
module load qbic/samtools/1.3
module load qbic/bwa
module load qbic/fastqc/0.11.4
module load qbic/picard/git

{exec_job}
exit 0
