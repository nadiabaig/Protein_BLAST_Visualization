#!/usr/bin/env bash
# Wrapper for Protein BLAST Visualization

BLAST_FILE="blast_results.tsv"
OUTDIR="out"

mkdir -p "$OUTDIR"

echo "Running BLAST visualization on $BLAST_FILE..."
python3 prot_blast.py \
  --in "$BLAST_FILE" \
  --outdir "$OUTDIR" \
  --topn 10 \
  --pdf \
  --xlsx

echo "Done. Outputs are in $OUTDIR/"
