#!/usr/bin/env python3
# run_15_cp6_fetch_conservation.py
# CP6-2: Check evolutionary conservation (phastCons100way) for target promoters

import os
import sys
import datetime
import pandas as pd
import pyBigWig
import time

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP6_DIR = os.path.join(PHASE_DIR, "CP6_Support_Evidence")
LOG_PATH = os.path.join(CP6_DIR, "logs", "run_15_cp6_conservation.log")

BW_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/rn7/phastCons100way/rn7.phastCons100way.bw"

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def get_mean_score(bw, chrom, start, end):
    try:
        # UCSC bigwig uses standard chr names (e.g., 'chr1')
        if chrom not in bw.chroms():
            return np.nan
        # Avoid out-of-bounds
        chr_len = bw.chroms(chrom)
        if start >= chr_len or end <= 0: return np.nan
        
        s = max(0, start)
        e = min(chr_len, end)
        
        if s >= e: return np.nan
        
        vals = bw.values(chrom, s, e)
        # Filter None
        valid_vals = [v for v in vals if v is not None and not np.isnan(v)]
        if len(valid_vals) == 0:
            return 0.0
        return sum(valid_vals) / len(valid_vals)
    except Exception as e:
        return np.nan

def main():
    os.makedirs(CP6_DIR, exist_ok=True)
    log("Started Phase H6 - Task 15 (CP6): Evolutionary Conservation (phastCons100way)")
    
    # Load Epigenomics results (cascading support)
    targets_path = os.path.join(CP6_DIR, "candidate_epigenomic_support.tsv")
    if not os.path.exists(targets_path):
        log(f"Error: {targets_path} not found. Run run_14 first.")
        sys.exit(1)
        
    df_targets = pd.read_csv(targets_path, sep='\t')
    log(f"Loaded {len(df_targets)} interactions.")
    
    log(f"Opening remote BigWig file: {BW_URL} (this may take a moment to initialize)")
    
    try:
        bw = pyBigWig.open(BW_URL)
    except Exception as e:
        log(f"Failed to open remote BigWig: {e}")
        log("Trying to download first... Not supported in this simplified script.")
        sys.exit(1)
        
    log("Successfully connected to remote BigWig.")
    
    scores = []
    total = len(df_targets)
    log("Querying mean conservation scores for target promoters (+/- 500bp)...")
    
    # Process unique promoters to save remote query time
    unique_promoters = df_targets[['t_chrom', 't_start', 't_end']].drop_duplicates()
    log(f"Found {len(unique_promoters)} unique promoters to query.")
    
    score_map = {}
    
    for i, row in unique_promoters.iterrows():
        c, s, e = row['t_chrom'], int(row['t_start']), int(row['t_end'])
        key = f"{c}:{s}-{e}"
        score_map[key] = get_mean_score(bw, c, s, e)
        
        if (i+1) % 500 == 0:
            log(f"  Queried {i+1} unique regions...")
            
    bw.close()
    
    # Map back
    df_targets['phastCons_mean'] = df_targets.apply(
        lambda r: score_map.get(f"{r['t_chrom']}:{int(r['t_start'])}-{int(r['t_end'])}", np.nan), axis=1
    )
    
    # Threshold for "conserved" (e.g. > 0.2 average conservation in 1kb is quite high for mammalian promoters)
    df_targets['is_conserved'] = df_targets['phastCons_mean'] > 0.1
    
    cons_count = df_targets['is_conserved'].sum()
    log(f"Found {cons_count} / {len(df_targets)} interactions with conserved promoters (mean > 0.1).")
    
    out_path = os.path.join(CP6_DIR, "candidate_conservation_support.tsv")
    df_targets.to_csv(out_path, sep='\t', index=False)
    log(f"Saved to {out_path}")
    log("Task 15 completed successfully.")

if __name__ == "__main__":
    import numpy as np
    main()
