# Protein BLAST Visualization Toolkit     <img width="100" height="100" alt="image" src="https://github.com/user-attachments/assets/5bac13ee-b5c2-4f39-9ec2-affb5b543ba3" />


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
```

##  How to Run the script
- **Install python packages:**
```bash
python3 -m pip install -U pandas numpy matplotlib xlsxwriter openpyxl
```
 - **Generate BLAST output in tabular format:**
   
```bash
blastp -query proteins.faa -db nr -out blast_results.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids sscinames salltitles" 
```
 - **Basic Run:**

```bash
python3 Script.py --in blast_results.tsc --outdir out --pdf   
```

 - **Write an excel sheet:**

```bash
python3 Script.py --in blast_results.tsc --outdir out --xlsx  
```

  - **Show top 10 blast hits per query and force coverage map for a specific query:**

```bash
python3 Script.py --in blast_results.tsc --outdir out --topn 10 --query My_query_ID --pdf --xlsx  
```
##  How to Run the script
Optionally you can use make file

 **Usage:**
 ```bash
make run
make pdf
make xlsx
make full
make clean
```




