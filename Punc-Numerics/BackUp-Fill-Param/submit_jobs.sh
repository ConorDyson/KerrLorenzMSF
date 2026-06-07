#!/bin/bash

# Format: "decimal:fraction" — decimal for log filenames, fraction passed to Mathematica
A_VALS=(
  "0.00000000000001:1/10^14"
  "0.6:3/5"
  "0.99:99/100"
)
R0_VALS=(
  "7.8:78/10"
  "7.9:79/10"
  "8.0:8"
  "8.1:81/10"
  "8.2:82/10"
)

# r grid parameters passed to Mathematica (exact fractions)
RMIN="79/10"
RMAX="81/10"
RSTEP="1/10"

LMAX=10
MMAX=2

JOBS_PER_NODE=1

SCRIPT_DIR="/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics"
# SAVEPATH="/lustre/astro/dyson/DataFiles/h1Puncs"
SAVEPATH="/lustre/astro/dyson/DataFiles/KerrCircSourceFiles/h1Puncs"

mkdir -p "${SCRIPT_DIR}/logs"

PAIRS_FILE="${SCRIPT_DIR}/logs/pairs.txt"
> "$PAIRS_FILE"
for a in "${A_VALS[@]}"; do
  for r0 in "${R0_VALS[@]}"; do
    echo "$a $r0" >> "$PAIRS_FILE"
  done
done

NTOTAL=$(wc -l < "$PAIRS_FILE")
echo "Total (a, r0) pairs: ${NTOTAL}"

JOBFILE="${SCRIPT_DIR}/logs/array_job.sh"
cat > "$JOBFILE" << EOF
#!/bin/bash
#SBATCH --job-name=h1_puncs
#SBATCH --time=00-11:59:00
#SBATCH --mem=500G
#SBATCH --partition=tycho1_short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --threads-per-core=1
#SBATCH --array=1-${NTOTAL}%$((NTOTAL < JOBS_PER_NODE ? NTOTAL : JOBS_PER_NODE))
#SBATCH --output=${SCRIPT_DIR}/logs/task_%a.txt
#SBATCH --error=${SCRIPT_DIR}/logs/task_%a.txt

PAIRS_FILE="${PAIRS_FILE}"
SCRIPT_DIR="${SCRIPT_DIR}"
SAVEPATH="${SAVEPATH}"

read -r a r0 <<< \$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\$PAIRS_FILE")
a_dec="\${a%%:*}"; a_frac="\${a##*:}"
r0_dec="\${r0%%:*}"; r0_frac="\${r0##*:}"

echo "Task \$SLURM_ARRAY_TASK_ID on \$(hostname): a=\$a_dec r0=\$r0_dec"

TUNNEL_LOCK="/tmp/mathtunnel_\${USER}_\$(hostname -s).lock"
sleep \$(( RANDOM % 30 ))
(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"\$TUNNEL_LOCK"

wolframscript "\${SCRIPT_DIR}/RunMathematicaPunctes.wls" "\$a_frac" "\$r0_frac" "${RMIN}" "${RMAX}" "${RSTEP}" "${LMAX}" "${MMAX}" "\$SAVEPATH"
MATH_EXIT=\$?

(
  flock -x 200
  sleep 5
  if ! pgrep -u "\$USER" -x WolframKernel > /dev/null 2>&1; then
    pkill -fu "\$USER" astro01 || true
  fi
) 200>"\$TUNNEL_LOCK"

exit \$MATH_EXIT
EOF

sbatch "$JOBFILE"
echo "Submitted array job (${NTOTAL} tasks, ${JOBS_PER_NODE} running at a time)"
