#!/bin/bash

#SBATCH --partition=norm # Job to run on partition 'partname'. (default: 'norm')
#SBATCH --ntasks=1 # Number of task to be run
#SBATCH --cpus-per-task=5 # Number of CPUs required for each task
#SBATCH --mem=200g # Memory required for the job. 
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END,FAIL  # Mail events
#SBATCH --mail-user=huang.lin@nih.gov   # Where to send mail 
#SBATCH --time=0-10:00:00 # Runtime in day-hour:min:sec
#SBATCH --requeue # If an allocated node hangs, the job should be requeued
#SBATCH --output=forehead_tax.out # Location of stdout file
#SBATCH --error=forehead_tax.err # Location of std class="softBottom"err file

pwd; hostname; date

echo "start my job"
module load QIIME/2-2021.4
echo "finish loading QIIME2"

conda activate QIIME/2-2021.4

qiime tools import \
  --input-path forehead_seqs.fa \
  --output-path forehead_seqs.qza \
  --type 'FeatureData[Sequence]'

qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads forehead_seqs.qza \
  --o-classification forehead_taxonomy.qza \
  --p-n-jobs -1 

echo "done"
date
