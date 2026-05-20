#!/bin/bash
#SBATCH --job-name=cq_fit
#SBATCH --array=0-39             # 0-19 for minus, 20-39 for plus
#SBATCH --output=logs/fit_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G                 # Adjust based on your needs
#SBATCH --time=00:20:00          # Adjust based on your run time

# 1. Logic for Charge (-c)
# Tasks 0-19 will be minus (-1), Tasks 20-39 will be plus (1)
if [ $SLURM_ARRAY_TASK_ID -lt 20 ]; then
    INDEX=$SLURM_ARRAY_TASK_ID
    SUFFIX="minus"
else
    INDEX=$((SLURM_ARRAY_TASK_ID - 20))
    SUFFIX="plus"
fi

# 2. Logic for Start Event (-s)
# Each index processes 30M events
START_HIST=${INDEX}

# 3. Path Setup
IN_DIR="/project/femtoscopy/students/toronyibandi/HEP/data/Q_full"
OUT_DIR="/project/femtoscopy/students/toronyibandi/HEP/data/fit_params"
mkdir -p $OUT_DIR

# 4. Execution
# We use the unique index and suffix for the filename
./../build/fitCQ \
    -i ${IN_DIR}/Q_${SUFFIX}.root\
    -s $START_HIST \
    -n 1 \
    -o ${OUT_DIR}/fitpars${INDEX}_${SUFFIX}.txt