#!/bin/sh
#BSUB -J run_preproc_AEP
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o bjob_log/AEP_Output_%J.txt
#BSUB -e bjob_lob/AEP_Error_%J.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd /work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocessing/_parprocessing 
matlab -nodisplay -r run_preproc_eeg_AEP -logfile log_preproc_aep

