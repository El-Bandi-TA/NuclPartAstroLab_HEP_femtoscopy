#!/bin/bash
#SBATCH --job-name=zphi_ana
#SBATCH --array=0-59             # 0-29 for minus, 30-59 for plus
#SBATCH --output=logs/ana_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G                 # Adjust based on your needs
#SBATCH --time=00:10:00          # Adjust based on your run time

# 1. Logic for Charge (-c)
# Tasks 0-29 will be minus (-1), Tasks 30-59 will be plus (1)
if [ $SLURM_ARRAY_TASK_ID -lt 30 ]; then
    CHARGE=-1
    INDEX=$SLURM_ARRAY_TASK_ID
    SUFFIX="minus"
else
    CHARGE=1
    INDEX=$((SLURM_ARRAY_TASK_ID - 30))
    SUFFIX="plus"
fi

# 2. Logic for Start Event (-s)
# Each index processes 30M events
START_EVENT=$((INDEX * 10000000))

# 3. Path Setup
OUT_DIR="/project/femtoscopy/students/toronyibandi/HEP/data/zphi_cut"
mkdir -p $OUT_DIR

# 4. Execution
# We use the unique index and suffix for the filename
./../build/analyzeCzphi \
    -s $START_EVENT \
    -n 10000000 \
    -c $CHARGE \
    --cut \
    -w -1 \
    -o ${OUT_DIR}/zphi_cut${INDEX}_${SUFFIX}.root