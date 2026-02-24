#!/usr/bin/env python3
import os
import subprocess

def modify_config_file(original_file, new_file, new_r0):
    """
    Reads the original configuration file, replaces the value of r0,
    and writes the result to new_file.
    """
    with open(original_file, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        # Replace the r0 line with the new value
        if line.strip().startswith('r0'):
            new_line = f"r0\t{new_r0}\n"
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    with open(new_file, 'w') as f:
        f.writelines(new_lines)

def get_n_cores(default=10):
    """
    Automatically set the number of cores from an environment variable.
    Here we use SLURM_CPUS_PER_TASK, which is common on SLURM clusters.
    If not found, fallback to a default value (e.g., 4).
    """
    return int(os.environ.get("SLURM_CPUS_PER_TASK", default))

def run_simulation(script_file, n_cores):
    """
    Runs the Mathematica Wolfram Language script with a limited number of cores.
    The environment is updated to restrict the number of threads.
    """
    # Copy the current environment and set thread limits
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(n_cores)
    env["MKL_NUM_THREADS"] = str(n_cores)
    
    # Use wolframscript to execute the Wolfram Language script
    cmd = f"/home/MATHEMATICA/12.3.1/bin/math -script {script_file}"
    print(f"Executing: {cmd} with OMP_NUM_THREADS={n_cores}")
    result = subprocess.run(cmd, shell=True, env=env)
    return result.returncode

def main():
    # File names for the config file and the Mathematica script
    original_config = "ConfigFile.txt"
    simulation_script = "h1-mmodes-script.wls"
    
    # Automatically determine the number of cores allocated on the cluster
    n_cores = get_n_cores(default=4)
    print(f"Using {n_cores} cores for this simulation.")

    # Define the range of r0 values you wish to run
    r0_values = [7.1, 7.2, 7.3, 7.4, 7.6, 7.7, 7.8, 7.9, 8.1, 8.2, 8.3, 8.4, 8.6, 8.7, 8.8, 8.9]

    for r0 in r0_values:
        print(f"\n=== Running simulation with r0 = {r0} ===")
        # Create a temporary config file for this r0 value
        temp_config = f"ConfigFile_r0_{r0}.txt"
        modify_config_file(original_config, temp_config, r0)
        
        # Replace the original config with the updated one
        os.replace(temp_config, original_config)
        
        # Run the simulation (Mathematica script) with the current config and n_cores
        ret = run_simulation(simulation_script, n_cores)
        if ret != 0:
            print(f"Simulation with r0 = {r0} failed with return code {ret}")
        else:
            print(f"Simulation with r0 = {r0} completed successfully")

if __name__ == "__main__":
    main()