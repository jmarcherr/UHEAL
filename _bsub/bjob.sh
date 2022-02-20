#!/bin/sh
#BSUB -J run_preproc
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o bjob_log/Output_%J.txt
#BSUB -e bjob_log/Error_%J.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd /work1/jonmarc/UHEAL_master/UHEAL/
matlab -nodisplay -r main -logfile log_preproc

