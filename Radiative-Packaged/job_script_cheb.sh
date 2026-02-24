#!/bin/bash
#SBATCH --job-name=h1_Lorentz # shows up in the output of squeue
#SBATCH --time=00-11:59:00       # specify the requested wall-time
#SBATCH --mem=99G
#SBATCH --partition=astro3_short #sba specify the partition to run on
#SBATCH --ntasks=1                  # number of tasks
#SBATCH --nodes=1               #request full node
#SBATCH --cpus-per-task=20      # number of OpenMP threads per MPI rank
#SBATCH --threads-per-core=1    # number of threads active per CPU core (=1: Hyperthreading off)
#SBATCH -o std_output.txt
#SBATCH -e std_error.txt

# Load Mathematica and setup port forwarding
ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf

# Execute the code
wolframscript  h1-mmodes-script.wls
  
pkill -fu $USER astro01
