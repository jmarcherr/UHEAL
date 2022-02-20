#!/bin/sh
#BSUB -J run_aep_analysis
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o bjob_log/aep_analysis_Output_%J.txt
#BSUB -e bjob_log/aep_analysis_Error_%J.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd /work1/jonmarc/UHEAL_master/UHEAL/_EEG/_analysis/par_analysis
matlab -nodisplay -r run_AEP_analysis_par -logfile log_preproc_aep_analysis

