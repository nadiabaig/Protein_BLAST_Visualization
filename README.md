# Protein BLAST Visualization Toolkit

`Script.py` is a Python script that turns **NCBI BLASTP results** into publication-ready **plots** and **summary tables**. 
It is designed for BLAST results saved in **tabular format (`-outfmt 6`)**, with optional extra fields for query length and taxonomy.

---

## ✨ Features

- **Plots (PNG, optional PDF):**
  - Histogram of % identity 
  - Histogram of query coverage (requires `qlen`) 
  - Scatterplot: % identity vs alignment length 
  - Taxonomic bar chart (requires `sscinames`) 
  - Coverage map for one query (requires `qlen qstart qend`) 

- **Tables (CSV, optional Excel):**
  - All hits (`all_hits_clean.csv`) 
  - Top N hits per query (`top_hits_per_query.csv`) 
  - Per-query summary (`query_summary.csv`) 
  - Per-taxon summary (`taxon_summary.csv`, if taxonomy included) 
  - Optional Excel workbook (`results.xlsx`) with **automatic chunking** if datasets exceed Excel’s row limits 

---

## 📥 Input

A BLAST tabular file (`-outfmt 6`). 
For best results, include these fields:

```bash
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids sscinames salltitles"
