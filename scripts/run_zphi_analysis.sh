#!/bin/bash
#SBATCH --job-name=zphi_ana
#SBATCH --array=0-19             # 0-9 for minus, 10-19 for plus
#SBATCH --output=logs/ana_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G                 # Adjust based on your needs
#SBATCH --time=00:45:00          # Adjust based on your run time

# 1. Logic for Charge (-c)
# Tasks 0-9 will be minus (-1), Tasks 10-19 will be plus (1)
if [ $SLURM_ARRAY_TASK_ID -lt 10 ]; then
    CHARGE=-1
    INDEX=$SLURM_ARRAY_TASK_ID
    SUFFIX="minus"
else
    CHARGE=1
    INDEX=$((SLURM_ARRAY_TASK_ID - 10))
    SUFFIX="plus"
fi

# 2. Logic for Start Event (-s)
# Each index processes 30M events
START_EVENT=$((INDEX * 30000000))

# 3. Path Setup
OUT_DIR="/project/femtoscopy/students/toronyibandi/HEP/data/zphi_full"
mkdir -p $OUT_DIR

# 4. Execution
# We use the unique index and suffix for the filename
./../build/analyzeCzphi \
    -s $START_EVENT \
    -n 30000000 \
    -w -1 \
    -c $CHARGE \
    -o ${OUT_DIR}/zphi${INDEX}_${SUFFIX}.root