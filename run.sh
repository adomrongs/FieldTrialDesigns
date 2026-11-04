#!/bin/bash
# =============================================================================
# run.sh: launches the simulation on Slurm, one job per replicate (job array).
# Run from the root of the repository:
#   bash run.sh burnin          burn-in
#   bash run.sh all             burn-in + every scenario (each replicate of a
#                               scenario starts when its burn-in has finished)
#   bash run.sh GS_pREP_odw     one scenario (needs the burn-in RData)
# =============================================================================

reps="1-30"
resources="-c 1 --mem=32G -t 24:00:00"
scenarios="PS_RCBD GS_RCBD GS_Augmented GS_pREP_odw"

# Inside the Slurm job: one replicate of one script
if [ "$1" == "--job" ]; then
  module load r-bio
  export OMP_NUM_THREADS=1 # sommer/odw: one thread per job (more threads were much slower)
  Rscript "code/scenarios/$2.R" "$SLURM_ARRAY_TASK_ID"
  exit $?
fi

mkdir -p logs
options="--parsable --array=$reps $resources -o logs/%x-%A_%a.out"

case $1 in
  burnin)
    sbatch $options -J BurnIn "$0" --job BurnIn
    ;;
  all)
    burnin_id=$(sbatch $options -J BurnIn "$0" --job BurnIn)
    echo "BurnIn: $burnin_id"
    for scenario in $scenarios; do
      echo "$scenario: $(sbatch $options --dependency=aftercorr:$burnin_id -J $scenario "$0" --job $scenario)"
    done
    ;;
  *)
    if [ ! -f "code/scenarios/$1.R" ]; then
      echo "Usage: bash run.sh burnin | all | <scenario> ($scenarios)"
      exit 1
    fi
    sbatch $options -J "$1" "$0" --job "$1"
    ;;
esac
