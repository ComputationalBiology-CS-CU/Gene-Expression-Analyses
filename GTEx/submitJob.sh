# Qsub command to submit the job for linearModelGtex.R
qsub -cwd -N linearModel -l h_data=200G,h_rt=48:00:00 -b y /nfs/apps/R/3.1.2/bin/Rscript linearModelGtex.R 
