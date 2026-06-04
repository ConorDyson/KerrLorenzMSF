#!/bin/bash
#SBATCH --job-name=h1_radial
#SBATCH --time=04-23:59:00
#SBATCH --mem=30G
#SBATCH --partition=astro2_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --array=1-15%15
#SBATCH --output=/groups/astro/dyson/Open/KerrLorenzMSF/Ret-Numerics/h1-radial-gen/logsTube/task_%a.txt
#SBATCH --error=/groups/astro/dyson/Open/KerrLorenzMSF/Ret-Numerics/h1-radial-gen/logsTube/task_%a.txt

PAIRS_FILE="/groups/astro/dyson/Open/KerrLorenzMSF/Ret-Numerics/h1-radial-gen/logsTube/pairs.txt"
SCRIPT_DIR="/groups/astro/dyson/Open/KerrLorenzMSF/Ret-Numerics/h1-radial-gen"
SAVEPATH="/lustre/astro/dyson/DataFiles/KerrCircSourceFiles/h1RadialTube"

# Read this task's (a, r0) pair
read -r a r0 <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PAIRS_FILE")

# Redirect stdout+stderr to a descriptive log file now that we know a and r0
LOG="${SCRIPT_DIR}/logsTube/h1_a${a}_r${r0}.log"
exec > "$LOG" 2>&1

echo "Task $SLURM_ARRAY_TASK_ID on $(hostname): a=$a r0=$r0"

# --- Shared SSH tunnel (one per node, shared across co-located tasks) ---
# Lock is node-scoped so tasks on different nodes don't share a lock file
TUNNEL_LOCK="/tmp/mathtunnel_${USER}_$(hostname -s).lock"

# Jitter so tasks landing on the same node don't all race to open the tunnel
sleep $(( RANDOM % 30 ))

(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"$TUNNEL_LOCK"

# Run the Mathematica job (quote args to handle scientific-notation values safely)
wolframscript "${SCRIPT_DIR}/RunFile.wls" "$a" "$r0" "$SAVEPATH"
MATH_EXIT=$?

# --- Tear down tunnel only when no WolframKernel remains on this node ---
(
  flock -x 200
  sleep 5  # grace period: let any sibling kernel finish its own shutdown
  if ! pgrep -u "$USER" -x WolframKernel > /dev/null 2>&1; then
    pkill -fu "$USER" astro01 || true
  fi
) 200>"$TUNNEL_LOCK"

exit $MATH_EXIT
