#!/bin/bash
#SBATCH --job-name=h1_kerr_test
#SBATCH --partition=oldgpu
#SBATCH --nodelist=pinatubo
#SBATCH --array=0-0           # Creates 1 tasks (one per config file)
#SBATCH --mem=200G                       # Request 128GB memory per task
#SBATCH --time=7-00:00:00               # Request 1 week of wall time
#SBATCH --cpus-per-task=16         # Dynamically set to number of MathKernels
#SBATCH --output=mmodes_%A_%a.out       # STDOUT log
#SBATCH --error=mmodes_%A_%a.err        # STDERR log
#SBATCH --mail-type=END,FAIL            # Email on job completion or failure
#SBATCH --mail-user=benjamin.leather@aei.mpg.de

# Load any necessary modules (e.g., module load gcc)

# Define an array of r0 values corresponding to each configuration file.
r0_values=(6.0)

# Determine the config file name for this array task.
CONFIG_FILE="generated_configs/ConfigFile_${r0_values[$SLURM_ARRAY_TASK_ID]}.txt"

# Verify that the config file exists before running.
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

# Run your script with the config file as an argument.
/home/MATHEMATICA/12.3.1/bin/math -script ./h1-mmodes-script.wls "$CONFIG_FILE"
