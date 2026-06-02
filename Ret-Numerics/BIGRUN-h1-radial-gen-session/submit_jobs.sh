#!/bin/bash

# Edit these lists to control which (a, r0) pairs get submitted
# A_VALS=(0.001 0.1 0.2 0.3 0.6 0.7 0.8 0.9 0.999)
# R0_VALS=(6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5 15.0 15.5 16.0 16.5 17.0 17.5 18.0 18.5 19.0 19.5 20.0 20.5 21.0 21.5 22.0 22.5 23.0 23.5 24.0 24.5 25.0 25.5 26.0 26.5 27.0 27.5 28.0 28.5 29.0 29.5 30.0 30.5 31.0 31.5 32.0 32.5 33.0 33.5 34.0 34.5 35.0 35.5 36.0 36.5 37.0 37.5 38.0 38.5 39.0 39.5 40.0 40.5 41.0 41.5 42.0 42.5 43.0 43.5 44.0 44.5 45.0 45.5 46.0 46.5 47.0 47.5 48.0 48.5 49.0 49.5 50.0)

A_VALS=(0.00000000000001 0.6 0.99)
R0_VALS=(7.8 7.9 8.0 8.1 8.2)


# Maximum number of array tasks running concurrently (SLURM distributes across nodes)
MAX_CONCURRENT=50

SCRIPT_DIR="/groups/astro/dyson/Open/KerrLorenzMSF/GenerationCodes/Numerics-WorldTube/h1-radial-gen"
SAVEPATH="/lustre/astro/dyson/DataFiles/KerrCircDiskData/h1_radial"


mkdir -p "${SCRIPT_DIR}/logs"

# Build the pairs list file (one "a r0" per line, 1-indexed by array task ID)
PAIRS_FILE="${SCRIPT_DIR}/logs/pairs.txt"
> "$PAIRS_FILE"
for a in "${A_VALS[@]}"; do
  for r0 in "${R0_VALS[@]}"; do
    echo "$a $r0" >> "$PAIRS_FILE"
  done
done

NTOTAL=$(wc -l < "$PAIRS_FILE")
echo "Total (a, r0) pairs: ${NTOTAL}"

# Write the array job script
JOBFILE="${SCRIPT_DIR}/logs/array_job.sh"
cat > "$JOBFILE" << EOF
#!/bin/bash
#SBATCH --job-name=h1_radial
#SBATCH --time=00-11:59:00
#SBATCH --mem=3G
#SBATCH --partition=astro2_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --array=1-${NTOTAL}%$((NTOTAL < MAX_CONCURRENT ? NTOTAL : MAX_CONCURRENT))
#SBATCH --output=${SCRIPT_DIR}/logs/task_%a.txt
#SBATCH --error=${SCRIPT_DIR}/logs/task_%a.txt

PAIRS_FILE="${PAIRS_FILE}"
SCRIPT_DIR="${SCRIPT_DIR}"
SAVEPATH="${SAVEPATH}"

# Read this task's (a, r0) pair
read -r a r0 <<< \$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\$PAIRS_FILE")

# Redirect stdout+stderr to a descriptive log file now that we know a and r0
LOG="\${SCRIPT_DIR}/logs/h1_a\${a}_r\${r0}.log"
AKEY="\$(printf '%.0f' \$(echo "\$a * 1000" | bc))"
RKEY="\$(printf '%.0f' \$(echo "\$r0 * 10" | bc))"
SAVELOG="\${SAVEPATH}/data\${AKEY}/data\${AKEY}-\${RKEY}/h1_a\${a}_r\${r0}.log"
mkdir -p "\$(dirname "\$SAVELOG")"
exec > >(tee "\$LOG" "\$SAVELOG") 2>&1

echo "Task \$SLURM_ARRAY_TASK_ID on \$(hostname): a=\$a r0=\$r0"

# --- Shared SSH tunnel (one per node, shared across co-located tasks) ---
# Lock is node-scoped so tasks on different nodes don't share a lock file
TUNNEL_LOCK="/tmp/mathtunnel_\${USER}_\$(hostname -s).lock"

# Jitter so tasks landing on the same node don't all race to open the tunnel
sleep \$(( RANDOM % 30 ))

(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"\$TUNNEL_LOCK"

# Run the Mathematica job (quote args to handle scientific-notation values safely)
wolframscript "\${SCRIPT_DIR}/RunFile.wls" "\$a" "\$r0" "\$SAVEPATH"
MATH_EXIT=\$?

# --- Tear down tunnel only when no WolframKernel remains on this node ---
(
  flock -x 200
  sleep 5  # grace period: let any sibling kernel finish its own shutdown
  if ! pgrep -u "\$USER" -x WolframKernel > /dev/null 2>&1; then
    pkill -fu "\$USER" astro01 || true
  fi
) 200>"\$TUNNEL_LOCK"

exit \$MATH_EXIT
EOF

sbatch "$JOBFILE"
echo "Submitted array job (${NTOTAL} tasks, ${MAX_CONCURRENT} running at a time)"
