#!/bin/sh
#BSUB -J run_analysis_midlate
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o bjob_log/Output.txt
#BSUB -e bjob_logError.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd /work1/jonmarc/UHEAL_master/UHEAL/_EEG/_analysis/par_analysis
matlab -nodisplay -r run_midlate_analysis_par -logfile log_analysis_midlate

