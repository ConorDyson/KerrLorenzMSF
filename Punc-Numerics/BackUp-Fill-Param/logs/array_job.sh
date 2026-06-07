#!/bin/bash
#SBATCH --job-name=h1_puncs
#SBATCH --time=00-11:59:00
#SBATCH --mem=500G
#SBATCH --partition=tycho1_short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --threads-per-core=1
#SBATCH --array=1-15%1
#SBATCH --output=/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/task_%a.txt
#SBATCH --error=/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/task_%a.txt

PAIRS_FILE="/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/pairs.txt"
SCRIPT_DIR="/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics"
SAVEPATH="/lustre/astro/dyson/DataFiles/KerrCircSourceFiles/h1Puncs"

read -r a r0 <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PAIRS_FILE")
a_dec="${a%%:*}"; a_frac="${a##*:}"
r0_dec="${r0%%:*}"; r0_frac="${r0##*:}"

echo "Task $SLURM_ARRAY_TASK_ID on $(hostname): a=$a_dec r0=$r0_dec"

TUNNEL_LOCK="/tmp/mathtunnel_${USER}_$(hostname -s).lock"
sleep $(( RANDOM % 30 ))
(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"$TUNNEL_LOCK"

wolframscript "${SCRIPT_DIR}/RunMathematicaPunctes.wls" "$a_frac" "$r0_frac" "79/10" "81/10" "1/10" "10" "2" "$SAVEPATH"
MATH_EXIT=$?

(
  flock -x 200
  sleep 5
  if ! pgrep -u "$USER" -x WolframKernel > /dev/null 2>&1; then
    pkill -fu "$USER" astro01 || true
  fi
) 200>"$TUNNEL_LOCK"

exit $MATH_EXIT
