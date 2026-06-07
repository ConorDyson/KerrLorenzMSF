#!/bin/bash
#SBATCH --job-name=punc_test
#SBATCH --time=01:00:00
#SBATCH --mem=30G
#SBATCH --partition=tycho1_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/test_%j.txt
#SBATCH --error=/groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/logs/test_%j.txt

TUNNEL_LOCK="/tmp/mathtunnel_${USER}_$(hostname -s).lock"
(
  flock -x 200
  if ! ss -tlnp 2>/dev/null | grep -q ':16286'; then
    ssh astro01.hpc.ku.dk -L16286:mathlm.nbi.dk:16286 -L16287:mathlm.nbi.dk:16287 -Nf
    sleep 3
  fi
) 200>"$TUNNEL_LOCK"

wolframscript /groups/astro/dyson/Open/KerrLorenzMSF/Punc-Numerics/RunMathematicaPunctes.wls
