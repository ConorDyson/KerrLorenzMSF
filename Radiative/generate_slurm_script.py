#!/usr/bin/env python3
import os
import subprocess

# Base configuration template (from your ConfigFile.txt)
# The {r0} placeholder will be replaced with each specific r0 value.
base_config_template = """a       0.6
r0      {r0}
mmax    1
tmax    250.0
n       10
angres  4
twr     16
twq     8
nu      0.800000
puncord 4
gv_opt  2
glg_opt 0
ic_opt  0
ic_ampl 0.0
lmax    30
lplot   0
nterms  8
inford  7
horord  6
rinf    4000.0
xhor    0.0001
rgrid   1
rmax    30.0
accgoal 16
kapord  8
prec    50
rstmin  2.0
rstmax  30.0
dformat Real64
"""

# # Slurm batch script template using a job array.
# slurm_script_template = """#!/bin/bash
# #SBATCH --job-name=h1_kerr
# #SBATCH --partition=oldgpu
# #SBATCH --nodelist=krakatoa
# #SBATCH --array=0-{max_index}           # Creates {total_jobs} tasks (one per config file)
# #SBATCH --mem=99G                       # Request 99GB memory per task
# #SBATCH --time=7-00:00:00               # Request 1 week of wall time
# #SBATCH --cpus-per-task={cpus}         # Dynamically set to number of MathKernels
# #SBATCH --output=mmodes_%A_%a.out       # STDOUT log
# #SBATCH --error=mmodes_%A_%a.err        # STDERR log
# #SBATCH --mail-type=END,FAIL            # Email on job completion or failure
# #SBATCH --mail-user=benjamin.leather@aei.mpg.de

# # Load any necessary modules (e.g., module load gcc)

# # Define an array of r0 values corresponding to each configuration file.
# r0_values=({r0_array})

# # Determine the config file name for this array task.
# CONFIG_FILE="generated_configs/ConfigFile_${{r0_values[$SLURM_ARRAY_TASK_ID]}}.txt"

# # Verify that the config file exists before running.
# if [ ! -f "$CONFIG_FILE" ]; then
#     echo "Error: Config file $CONFIG_FILE not found!"
#     exit 1
# fi

# # Run your script with the config file as an argument.
# /home/MATHEMATICA/12.3.1/bin/math -script ./h1-mmodes-script.wls "$CONFIG_FILE"
# """

# def main():
#     # Specify the set of r0 values for each run.
#     # Change or extend this list to match your desired r0 values.
#     # r0_values = [8.100000000000001, 8.2, 8.3, 8.4, 8.5]
#     r0_values = [8., 8.100000000000001, 8.2, 8.3, 8.4]
#     # r0_values = [8.5, 8.600000000000001, 8.7, 8.8, 8.9]

#     # Create an output directory for generated files (or use current directory)
#     output_dir = "generated_configs"
#     os.makedirs(output_dir, exist_ok=True)

#     # Generate individual configuration files with r0 value in the filename.
#     for r0 in r0_values:
#         config_content = base_config_template.format(r0=r0)
#         # File name includes the r0 value (e.g., "ConfigFile_8.0.txt")
#         config_filename = os.path.join(output_dir, f"ConfigFile_{r0}.txt")
#         with open(config_filename, "w") as f:
#             f.write(config_content)
#         print(f"Generated {config_filename}")

#     # Prepare the r0 values bash array as a space-separated string.
#     # This preserves the decimal values.
#     r0_array_str = " ".join(str(r0) for r0 in r0_values)

#     # Generate the Slurm batch script.
#     total_jobs = len(r0_values)
#     max_index = total_jobs - 1
#     slurm_script_content = slurm_script_template.format(
#         max_index=max_index,
#         total_jobs=total_jobs,
#         r0_array=r0_array_str
#     )
#     slurm_script_filename = os.path.join(output_dir, "run_mmodes.sh")
#     with open(slurm_script_filename, "w") as f:
#         f.write(slurm_script_content)
#     print(f"Generated {slurm_script_filename}")

# if __name__ == "__main__":
#     main()

# Function to query Mathematica for the number of available MathKernels
# Uses $ProcessorCount as a proxy for license-allowed kernel count
# Adjust the Wolfram Language expression if you have a licensing function available

def get_mathkernel_count():
    wl_expr = 'Print[$ProcessorCount]; Exit[];'
    try:
        output = subprocess.check_output([
            '/home/MATHEMATICA/12.3.1/bin/math', '-noprompt', '-run', wl_expr
        ], stderr=subprocess.DEVNULL, text=True)
        return int(output.strip())
    except Exception as e:
        print(f"Warning: could not determine MathKernel count ({e}); defaulting to 4")
        return 4


def main():
    # Specify the set of r0 values for each run.
    r0_values = [6.0]

    # Detect number of MathKernels and set CPUs accordingly
    cpus = get_mathkernel_count()
    print(f"Using {cpus} CPUs per task (MathKernels available)")

    # Create an output directory for generated files
    output_dir = "generated_configs"
    os.makedirs(output_dir, exist_ok=True)

    # Generate individual configuration files with r0 value in the filename.
    for r0 in r0_values:
        config_content = base_config_template.format(r0=r0)
        config_filename = os.path.join(output_dir, f"ConfigFile_{r0}.txt")
        with open(config_filename, "w") as f:
            f.write(config_content)
        print(f"Generated {config_filename}")

    # Prepare the r0 values bash array as a space-separated string.
    r0_array_str = " ".join(str(r0) for r0 in r0_values)

    # Determine job array indices
    total_jobs = len(r0_values)
    max_index = total_jobs - 1

    # Build the Slurm batch script with an f-string to embed variables directly
    slurm_script_content = f"""#!/bin/bash
#SBATCH --job-name=h1_kerr_test
#SBATCH --partition=oldgpu
#SBATCH --nodelist=pinatubo
#SBATCH --array=0-{max_index}           # Creates {total_jobs} tasks (one per config file)
#SBATCH --mem=200G                       # Request 128GB memory per task
#SBATCH --time=7-00:00:00               # Request 1 week of wall time
#SBATCH --cpus-per-task={cpus}         # Dynamically set to number of MathKernels
#SBATCH --output=mmodes_%A_%a.out       # STDOUT log
#SBATCH --error=mmodes_%A_%a.err        # STDERR log
#SBATCH --mail-type=END,FAIL            # Email on job completion or failure
#SBATCH --mail-user=benjamin.leather@aei.mpg.de

# Load any necessary modules (e.g., module load gcc)

# Define an array of r0 values corresponding to each configuration file.
r0_values=({r0_array_str})

# Determine the config file name for this array task.
CONFIG_FILE="generated_configs/ConfigFile_${{r0_values[$SLURM_ARRAY_TASK_ID]}}.txt"

# Verify that the config file exists before running.
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file $CONFIG_FILE not found!"
    exit 1
fi

# Run your script with the config file as an argument.
/home/MATHEMATICA/12.3.1/bin/math -script ./h1-mmodes-script.wls "$CONFIG_FILE"
"""

    # Write the Slurm script to file
    slurm_script_filename = os.path.join(output_dir, "run_mmodes.sh")
    with open(slurm_script_filename, "w") as f:
        f.write(slurm_script_content)
    print(f"Generated {slurm_script_filename}")


if __name__ == "__main__":
    main()