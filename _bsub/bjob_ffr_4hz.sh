#!/bin/sh
#BSUB -J run_preproc_ffr_4Hz
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o bjob_log/Output_4_hz_%J.txt
#BSUB -e bjob_lob/Error_4Hz_%J.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd /work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocessing/_parprocessing 
matlab -nodisplay -r run_preproc_eeg_FFR_4Hz -logfile log_preproc_ffr_4hz

