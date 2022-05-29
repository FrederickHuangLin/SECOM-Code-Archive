#!/bin/bash

#SBATCH --partition=largemem # Job to run on partition 'partname'. (default: 'norm')
#SBATCH --ntasks=1 # Number of task to be run
#SBATCH --cpus-per-task=20 # Number of CPUs required for each task
#SBATCH --mem=500g # Memory required for the job. 
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END,FAIL  # Mail events
#SBATCH --mail-user=huang.lin@nih.gov   # Where to send mail 
#SBATCH --time=2-00:00:00 # Runtime in day-hour:min:sec
#SBATCH --requeue # If an allocated node hangs, the job should be requeued
#SBATCH --output=se_mb.out # Location of stdout file
#SBATCH --error=se_mb.err # Location of std class="softBottom"err file

pwd; hostname; date

echo "start my job"
module load R/4.1.0
echo "finish loading R"

Rscript se_mb.R
echo "done"
date