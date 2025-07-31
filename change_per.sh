#!/bin/bash

spn=$1
scout="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
tumourevodir="${scout}/${spn}/sequencing/tumourevo"
sarekdir="${scout}/${spn}/sequencing/sarek"

for file in "${tumourevodir}"/tumourevo_*.sh; do
  # Skip if no matching files
  [ -e "$file" ] || continue

  echo "Processing: $file"

  # Find the line number where 'nextflow' appears
  l=$(grep -n -w nextflow "$file" | cut -f1 -d ":")

  # Only proceed if the line was found
  if [ -n "$l" ]; then
    # Insert 2 lines before
    sed -i "${l}i sg cdslab << EOF\numask 007" "$file"

    # Adjust line number after insertion: 'nextflow' is now 2 lines below
    l=$((l + 2))

    # Insert 1 line after
    sed -i "${l}a EOF" "$file"
  else
    echo "Warning: 'nextflow' not found in $file"
  fi
done


for file in "${sarekdir}"/sequenza_*.sh; do
  # Skip if no matching files
  [ -e "$file" ] || continue

  echo "Processing: $file"

  # Find the line number where 'nextflow' appears
  l=$(grep -n -w nextflow "$file" | cut -f1 -d ":")

  # Only proceed if the line was found
  if [ -n "$l" ]; then
    # Insert 2 lines before
    sed -i "${l}i sg cdslab << EOF\numask 007" "$file"

    # Adjust line number after insertion: 'nextflow' is now 2 lines below
    l=$((l + 2))

    # Insert 1 line after
    sed -i "${l}a EOF" "$file"
  else
    echo "Warning: 'nextflow' not found in $file"
  fi
done


for file in "${sarekdir}"/sarek_*.sh; do
  # Skip if no matching files
  [ -e "$file" ] || continue

  echo "Processing: $file"

  # Find the line number where 'nextflow' appears
  l=$(grep -n -w nextflow "$file" | cut -f1 -d ":")

  # Only proceed if the line was found
  if [ -n "$l" ]; then
    # Insert 2 lines before
    sed -i "${l}i sg cdslab << EOF\numask 007" "$file"

    # Adjust line number after insertion: 'nextflow' is now 2 lines below
    l=$((l + 2))

    # Insert 1 line after
    sed -i "${l}a EOF" "$file"
  else
    echo "Warning: 'nextflow' not found in $file"
  fi
done
