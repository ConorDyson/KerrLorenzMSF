#!/bin/bash
# For each (a, r0) data directory: merge per-m .h5 files into a single file,
# delete the per-m .h5 files, and remove old per-m config files.

SAVEPATH="/lustre/astro/dyson/DataFiles/KerrCircDiskData/h1_radial"

for datadir in "$SAVEPATH"/data*/data*-*/; do
  perM_files=( "$datadir"data/h1_a*_rp*_l*_m*.h5 )

  # Skip if no per-m files exist
  [[ -e "${perM_files[0]}" ]] || continue

  # Extract a and r0 from the first filename
  first="${perM_files[0]}"
  basename=$(basename "$first")
  # filename: h1_a{a}_rp{r0}_l{lmax}_m{m}.h5
  a=$(echo "$basename"   | sed 's/h1_a\(.*\)_rp.*/\1/')
  r0=$(echo "$basename"  | sed 's/h1_a.*_rp\(.*\)_l.*/\1/')

  merged="$datadir/h1_a${a}_r${r0}.h5"
  echo "Processing: $datadir  ->  $(basename $merged)"

  # Merge each per-m group into the combined file
  for f in "${perM_files[@]}"; do
    m=$(echo "$(basename $f)" | sed 's/.*_m\(.*\)\.h5/\1/')
    group="/m_${m}"
    h5copy -i "$f" -o "$merged" -s "$group" -d "$group"
  done

  # Verify the merged file was created before deleting
  if [[ -f "$merged" ]]; then
    rm -f "${perM_files[@]}"
    echo "  Merged ${#perM_files[@]} files, deleted originals"
  else
    echo "  ERROR: merged file not created, skipping deletion for $datadir"
    continue
  fi

  # Remove old per-m config files (keep only the single meta config)
  find "$datadir/config" -name "config*-0*.txt" -delete
  echo "  Removed per-m config files"

done

echo "Done."
