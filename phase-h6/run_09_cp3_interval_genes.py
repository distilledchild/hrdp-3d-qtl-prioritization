#!/usr/bin/env python3
# run_09_cp3_interval_genes.py
# CP3: Find all genes intersecting QTL intervals, and genes within +/- 1Mb distance.
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
LOG_PATH = os.path.join(CP3_DIR, "logs", "run_09_cp3_interval_genes.log")

WINDOW_SIZE = 1_000_000 # 1 Mb

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def find_overlapping_genes(qtl_row, df_genes, expand_bp=0):
    chrom = qtl_row['chrom']
    q_start = qtl_row['start'] - expand_bp
    q_end = qtl_row['end'] + expand_bp
    
    chrom_genes = df_genes[df_genes['chrom'] == chrom]
    if chrom_genes.empty:
        return []
        
    # Overlap condition: max(start1, start2) <= min(end1, end2)
    # gene_start <= q_end AND gene_end >= q_start
    mask = (chrom_genes['start'] <= q_end) & (chrom_genes['end'] >= q_start)
    overlaps = chrom_genes[mask]
    
    results = []
    for _, g in overlaps.iterrows():
        # calculate exact distance to original unexpanded interval
        orig_start = qtl_row['start']
        orig_end = qtl_row['end']
        
        # If overlaps original interval, distance is 0
        if g['start'] <= orig_end and g['end'] >= orig_start:
            dist = 0
        else:
            # Shortest distance between [g_start, g_end] and [orig_start, orig_end]
            dist = max(0, max(orig_start - g['end'], g['start'] - orig_end))
            
        results.append({
            'qtl_name': qtl_row['name'],
            'qtl_chrom': chrom,
            'qtl_start': orig_start,
            'qtl_end': orig_end,
            'gene_name': g['name'],
            'gene_start': g['start'],
            'gene_end': g['end'],
            'dist_to_interval': dist,
            'is_within_interval': dist == 0
        })
    return results

def process_qtl_bed(bed_path, df_genes, out_prefix):
    if not os.path.exists(bed_path):
        log(f"Warning: {bed_path} not found.")
        return
        
    df_qtl = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    log(f"Loaded {len(df_qtl)} QTLs from {os.path.basename(bed_path)}")
    
    all_results = []
    for _, row in df_qtl.iterrows():
        # find within +/- 1Mb
        hits = find_overlapping_genes(row, df_genes, expand_bp=WINDOW_SIZE)
        all_results.extend(hits)
        
    df_res = pd.DataFrame(all_results)
    if df_res.empty:
        log("No overlaps found at all.")
        return
        
    # Split into strict interval vs distance-weighted
    df_interval = df_res[df_res['is_within_interval'] == True].copy()
    
    interval_out = os.path.join(CP3_DIR, f"interval_genes_{out_prefix}.tsv")
    dist_out = os.path.join(CP3_DIR, f"distance_weighted_genes_{out_prefix}.tsv")
    
    df_interval.to_csv(interval_out, sep='\t', index=False)
    df_res.to_csv(dist_out, sep='\t', index=False)
    
    log(f"Found {len(df_interval)} gene-QTL pairs strictly within interval.")
    log(f"Found {len(df_res)} gene-QTL pairs within +/- 1Mb window.")
    log(f"Saved {interval_out}")
    log(f"Saved {dist_out}")

def main():
    log("Started Phase H6 - Task 9 (CP3): Interval & Distance-Weighted Gene Mapping")
    
    gene_path = os.path.join(REF_DIR, "rn7_protein_coding_genes.bed")
    if not os.path.exists(gene_path):
        log(f"Error: Gene bed not found at {gene_path}. Run Task 7 first.")
        sys.exit(1)
        
    df_genes = pd.read_csv(gene_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    log(f"Loaded {len(df_genes)} protein-coding genes.")
    
    primary_bed = os.path.join(CP2_DIR, "qtl_filtered_primary.bed")
    process_qtl_bed(primary_bed, df_genes, "primary")
    
    secondary_bed = os.path.join(CP2_DIR, "qtl_broad_secondary.bed")
    process_qtl_bed(secondary_bed, df_genes, "secondary")
    
    log("Task 9 completed successfully.")

if __name__ == "__main__":
    main()
