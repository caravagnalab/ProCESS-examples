#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=validation.log
#SBATCH --job-name=validation

module load R/4.4.1

# to modify
SPN="SPN01"
DIR="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation"
# combinations you want to validate
COVERAGE=(50 100)
PURITY=(0.3 0.6 0.9)

# to keep as it is
QUEUE=$(sinfo -h -o "%P %a %D %t" | grep -w 'EPYC\|GENOA\|THIN' |awk '$2 == "up" && $4 ~ /idle|mix/ {print $1, $3}' | awk '{sum[$1] += $2} END {for (p in sum) print p, sum[p]}' | sort -k2 -nr | head -n1 | cut -f1 -d " ")
ACCOUNT="cdslab"

# validate normal
echo "Submitting job for normal"
JOBID=$(sbatch --parsable \
               --job-name=validate_comb_${COV}_${PUR} \
               --partition=$QUEUE \
               --nodes=1 \
               --cpus-per-task=1 \
               --mem=24GB \
               --time=4:00:00 \
               --output=out/val_normal.out \
               --error=out/val_normal.err \
               --wrap="$DIR/validate_normal.py $SPN 50 0.3 $DIR -P $QUEUE -A $ACCOUNT")


if [ $? -eq 0 ]; then
    echo "Job submitted successfully: $JOBID"
else
    echo "Failed to submit job for COV=$COV, PUR=$PUR"
    echo "Error: $JOBID"
fi

job_ids_single_combination=()
# validate combination
for COV in "${COVERAGE[@]}"; do
  for PUR in "${PURITY[@]}"; do
    echo "Submitting job for COV=$COV, PUR=$PUR"

    JOBID=$(sbatch --parsable \
                   --job-name=validate_comb_${COV}_${PUR} \
                   --partition=$QUEUE \
                   --nodes=1 \
                   --cpus-per-task=1 \
                   --mem=24GB \
                   --time=4:00:00 \
                   --output=out/val_${COV}_${PUR}.out \
                   --error=out/val_${COV}_${PUR}.err \
                   --wrap="$DIR/validate_combination.py $SPN $COV $PUR $DIR -P $QUEUE -A $ACCOUNT")
    job_ids_single_combination+=("$JOBID") 
    if [ $? -eq 0 ]; then
        echo "Job submitted successfully: $JOBID"
    else
        echo "Failed to submit job for COV=$COV, PUR=$PUR"
        echo "Error: $JOBID"
    fi
  done
done

job_ids_all_combination=$(IFS=, ; echo "${job_ids_single_combination[*]}")
echo "Submitting job for $SPN"

COVERAGE_STR=$(IFS=, ; echo "${COVERAGE[*]}" | sed 's/,/,/g')
PURITY_STR=$(IFS=, ; echo "${PURITY[*]}" | sed 's/,/,/g')

CMD="$DIR/validate_SPN.py $SPN ${COVERAGE_STR} ${PURITY_STR} $DIR -P $QUEUE -A $ACCOUNT"
echo $CMD
JOBID=$(sbatch --parsable \
               --dependency=afterok:$job_ids_all_combination \
               --job-name=validate_${SPN} \
               --partition=$QUEUE \
               --nodes=1 \
               --cpus-per-task=1 \
               --mem=24GB \
               --time=4:00:00 \
               --output=out/val_${SPN}.out \
               --error=out/val_${SPN}.err \
               --wrap="$CMD")


if [ $? -eq 0 ]; then
    echo "Job submitted successfully: $JOBID"
else
    echo "Failed to submit job for SPN=$SPN"
    echo "Error: $JOBID"
fi
