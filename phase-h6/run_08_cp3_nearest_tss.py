#!/usr/bin/env python3
# run_08_cp3_nearest_tss.py
# CP3: Compute the nearest TSS (Ensembl rn7) for each QTL interval using Pandas.
import os
import sys
import datetime
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP2_DIR = os.path.join(PHASE_DIR, "CP2_QTL_Universe")
CP3_DIR = os.path.join(PHASE_DIR, "CP3_Positional_Baseline")
REF_DIR = os.path.join(CP3_DIR, "reference")
LOG_PATH = os.path.join(CP3_DIR, "logs", "run_08_cp3_nearest_tss.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def find_nearest_tss(qtl_row, df_tss):
    chrom = qtl_row['chrom']
    q_start = qtl_row['start']
    q_end = qtl_row['end']
    
    # Subset TSS to the same chromosome
    chrom_tss = df_tss[df_tss['chrom'] == chrom]
    if chrom_tss.empty:
        return pd.Series({'nearest_gene': 'None', 'tss_dist': np.nan})
        
    # TSS position (since it's 1bp, start == end-1, we can just use start)
    t_pos = chrom_tss['start'].values
    
    # Calculate distances
    # Distance to interval [q_start, q_end]
    # If t_pos < q_start: q_start - t_pos
    # If t_pos > q_end: t_pos - q_end
    # Else (inside interval): 0
    dists = np.maximum(0, np.maximum(q_start - t_pos, t_pos - q_end))
    
    min_idx = np.argmin(dists)
    nearest_dist = dists[min_idx]
    nearest_gene = chrom_tss.iloc[min_idx]['name']
    
    return pd.Series({'nearest_gene': nearest_gene, 'tss_dist': nearest_dist})

def process_qtl_bed(bed_path, df_tss, out_path):
    if not os.path.exists(bed_path):
        log(f"Warning: {bed_path} not found.")
        return
        
    df_qtl = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    log(f"Loaded {len(df_qtl)} QTLs from {os.path.basename(bed_path)}")
    
    # Apply nearest TSS calculation
    res = df_qtl.apply(lambda row: find_nearest_tss(row, df_tss), axis=1)
    df_result = pd.concat([df_qtl, res], axis=1)
    
    df_result.to_csv(out_path, sep='\t', index=False)
    log(f"Saved {out_path}")

def main():
    log("Started Phase H6 - Task 8 (CP3): Nearest TSS Mapping")
    
    tss_path = os.path.join(REF_DIR, "rn7_protein_coding_tss.bed")
    if not os.path.exists(tss_path):
        log(f"Error: TSS bed not found at {tss_path}. Run Task 7 first.")
        sys.exit(1)
        
    df_tss = pd.read_csv(tss_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    log(f"Loaded {len(df_tss)} TSS records.")
    
    primary_bed = os.path.join(CP2_DIR, "qtl_filtered_primary.bed")
    primary_out = os.path.join(CP3_DIR, "nearest_tss_primary.tsv")
    process_qtl_bed(primary_bed, df_tss, primary_out)
    
    secondary_bed = os.path.join(CP2_DIR, "qtl_broad_secondary.bed")
    secondary_out = os.path.join(CP3_DIR, "nearest_tss_secondary.tsv")
    process_qtl_bed(secondary_bed, df_tss, secondary_out)
    
    log("Task 8 completed successfully.")

if __name__ == "__main__":
    main()
