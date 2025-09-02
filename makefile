# Makefile for Protein BLAST Visualization

# Input BLAST file and output directory
BLAST_FILE = blast_results.tsv
OUTDIR     = out

# Default target
all: run

# Run the visualization script (basic: PNG + CSV)
run:
        python3 prot_blast.py --in $(BLAST_FILE) --outdir $(OUTDIR)

# Run with PDF plots too
pdf:
        python3 prot_blast.py --in $(BLAST_FILE) --outdir $(OUTDIR) --pdf

# Run with Excel workbook too
xlsx:
        python3 prot_blast.py --in $(BLAST_FILE) --outdir $(OUTDIR) --xlsx

# Run with everything enabled
full:
        python3 prot_blast.py --in $(BLAST_FILE) --outdir $(OUTDIR) --pdf --xlsx --topn 10

# Clean outputs
clean:
        rm -rf $(OUTDIR)/*
