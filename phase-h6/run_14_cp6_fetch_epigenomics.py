#!/usr/bin/env python3
# run_14_cp6_fetch_epigenomics.py
# CP6-1: Compare 3D-assigned target loops against ATAC-seq peaks to find active regulatory regions.

import os
import sys
import datetime
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP4_DIR = os.path.join(PHASE_DIR, "CP4_3D_Mapping")
CP6_DIR = os.path.join(PHASE_DIR, "CP6_Support_Evidence")
RES_DIR = os.path.join(CP6_DIR, "resources")
LOG_PATH = os.path.join(CP6_DIR, "logs", "run_14_cp6_epigenomics.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def get_overlaps(df_a, df_b, a_prefix, b_prefix):
    overlaps = []
    chroms = set(df_a[f'{a_prefix}chrom']).intersection(set(df_b[f'{b_prefix}chrom']))
    for chrom in chroms:
        a_sub = df_a[df_a[f'{a_prefix}chrom'] == chrom]
        b_sub = df_b[df_b[f'{b_prefix}chrom'] == chrom]
        
        a_starts = a_sub[f'{a_prefix}start'].values[:, None]
        a_ends = a_sub[f'{a_prefix}end'].values[:, None]
        
        b_starts = b_sub[f'{b_prefix}start'].values[None, :]
        b_ends = b_sub[f'{b_prefix}end'].values[None, :]
        
        mask = (a_starts <= b_ends) & (a_ends >= b_starts)
        a_idx_rel, b_idx_rel = np.where(mask)
        
        a_idx_abs = a_sub.index[a_idx_rel].values
        b_idx_abs = b_sub.index[b_idx_rel].values
        
        for i, j in zip(a_idx_abs, b_idx_abs):
            overlaps.append((i, j))
    return overlaps

def main():
    os.makedirs(CP6_DIR, exist_ok=True)
    log("Started Phase H6 - Task 14 (CP6): Epigenomic Evidence")
    
    # Load Primary 3D Targets
    targets_path = os.path.join(CP4_DIR, "qtl_loop_target_genes_primary.tsv")
    df_targets = pd.read_csv(targets_path, sep='\t')
    log(f"Loaded {len(df_targets)} primary 3D interactions.")
    
    # Load ATAC-seq peaks (rn7)
    df_atac_list = []
    
    atac_path1 = os.path.join(RES_DIR, "rn7_GSE134935_ATAC.bed")
    if os.path.exists(atac_path1):
        df_atac_list.append(pd.read_csv(atac_path1, sep='\t', header=None, usecols=[0,1,2], names=['a_chrom', 'a_start', 'a_end']))
        
    atac_path2 = os.path.join(RES_DIR, "rn7_Yuan2021_ATAC.bed")
    if os.path.exists(atac_path2):
        df_atac_list.append(pd.read_csv(atac_path2, sep='\t', header=None, usecols=[0,1,2], names=['a_chrom', 'a_start', 'a_end']))
        
    if not df_atac_list:
        log("Warning: No ATAC-seq files found. Will proceed with 0 ATAC support.")
        df_atac = pd.DataFrame(columns=['a_chrom', 'a_start', 'a_end'])
    else:
        df_atac = pd.concat(df_atac_list, ignore_index=True)
        log(f"Loaded {len(df_atac)} combined ATAC-seq peaks from GSE134935 and Yuan2021.")
    
    # We want to check if the QTL-hit Distal Anchor overlaps with an ATAC peak
    # Extract the anchor coordinates that the QTL hit
    def get_hit_chrom(row):
        return row['anchor1_chrom'] if row['qtl_hit_anchor'] == 'Anchor1' else row['anchor2_chrom']
        
    def get_hit_start(row):
        return row['anchor1_start'] if row['qtl_hit_anchor'] == 'Anchor1' else row['anchor2_start']
        
    def get_hit_end(row):
        return row['anchor1_end'] if row['qtl_hit_anchor'] == 'Anchor1' else row['anchor2_end']
        
    df_targets['t_chrom'] = df_targets.apply(get_hit_chrom, axis=1)
    df_targets['t_start'] = df_targets.apply(get_hit_start, axis=1)
    df_targets['t_end'] = df_targets.apply(get_hit_end, axis=1)
    
    log("Computing overlaps between QTL-hit Distal Anchors and ATAC peaks...")
    overlaps = get_overlaps(df_targets, df_atac, 't_', 'a_')
    
    # Mark targets with ATAC support
    supported_idx = set([i for i, j in overlaps])
    df_targets['atac_support'] = [1 if i in supported_idx else 0 for i in df_targets.index]
    
    atac_count = df_targets['atac_support'].sum()
    log(f"Found {atac_count} / {len(df_targets)} interactions with ATAC-seq support at the Distal Anchor.")
    
    out_path = os.path.join(CP6_DIR, "candidate_epigenomic_support.tsv")
    df_targets.to_csv(out_path, sep='\t', index=False)
    log(f"Saved to {out_path}")
    log("Task 14 completed successfully.")

if __name__ == "__main__":
    main()
