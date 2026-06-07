#!/bin/bash
#SBATCH --job-name=h1_punc_single
#SBATCH --time=02-23:59:00
#SBATCH --mem=100G
#SBATCH --partition=tycho1_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --threads-per-core=1
#SBATCH --output=/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/single_a0.6_r8.2_l45_m20.txt
#SBATCH --error=/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/single_a0.6_r8.2_l45_m20.txt

echo "Job started on $(hostname): a=0.6 r0=8.2 lmax=45 mmax=20"

TUNNEL_LOCK="/tmp/mathtunnel_${USER}_$(hostname -s).lock"
(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"$TUNNEL_LOCK"

wolframscript "/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/RunMathematicaPunctesBackUpWorking.wls"   "3/5" "82/10" "30/10" "200/10" "1/10" "45" "20" "/lustre/astro/dyson/DataFiles/KerrCircSourceFiles/h1Puncs"
MATH_EXIT=$?

(
  flock -x 200
  sleep 5
  if ! pgrep -u "$USER" -x WolframKernel > /dev/null 2>&1; then
    pkill -fu "$USER" astro01 || true
  fi
) 200>"$TUNNEL_LOCK"

exit $MATH_EXIT
