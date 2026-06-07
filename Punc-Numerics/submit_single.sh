#!/bin/bash

# ── Parameters — edit these ──
A="3/5"           # spin (exact fraction)
A_DEC="0.6"       # spin (decimal, for filenames/logs)
R0="82/10"            # orbital radius (exact fraction)
R0_DEC="8.2"      # orbital radius (decimal, for filenames/logs)
RMIN="30/10"      # r grid min
RMAX="200/10"      # r grid max
RSTEP="1/10"      # r grid step
LMAX=45
MMAX=20
SAVEPATH="/lustre/astro/dyson/DataFiles/KerrCircSourceFiles/h1Puncs"

# ── Paths ──
SCRIPT_DIR="/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics"
mkdir -p "${SCRIPT_DIR}/logs"

JOBFILE="${SCRIPT_DIR}/logs/single_job.sh"
cat > "$JOBFILE" << EOF
#!/bin/bash
#SBATCH --job-name=h1_punc_single
#SBATCH --time=02-23:59:00
#SBATCH --mem=100G
#SBATCH --partition=tycho1_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --threads-per-core=1
#SBATCH --output=${SCRIPT_DIR}/logs/single_a${A_DEC}_r${R0_DEC}_l${LMAX}_m${MMAX}.txt
#SBATCH --error=${SCRIPT_DIR}/logs/single_a${A_DEC}_r${R0_DEC}_l${LMAX}_m${MMAX}.txt

echo "Job started on \$(hostname): a=${A_DEC} r0=${R0_DEC} lmax=${LMAX} mmax=${MMAX}"

TUNNEL_LOCK="/tmp/mathtunnel_\${USER}_\$(hostname -s).lock"
(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"\$TUNNEL_LOCK"

wolframscript "${SCRIPT_DIR}/RunMathematicaPunctesBackUpWorking.wls" \
  "${A}" "${R0}" "${RMIN}" "${RMAX}" "${RSTEP}" "${LMAX}" "${MMAX}" "${SAVEPATH}"
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
echo "Submitted single job: a=${A_DEC} r0=${R0_DEC} lmax=${LMAX} mmax=${MMAX}"
