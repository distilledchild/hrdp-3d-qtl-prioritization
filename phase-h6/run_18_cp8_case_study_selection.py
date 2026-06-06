#!/usr/bin/env python3
# run_18_cp8_case_study_selection.py
# CP8: Select Case Study Loci and generate Virtual 4C bedGraph tracks

import os
import sys
import datetime
import pandas as pd
import numpy as np
import cooler

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP7_DIR = os.path.join(PHASE_DIR, "CP7_Scoring")
CP8_DIR = os.path.join(PHASE_DIR, "CP8_Case_Studies")
LOG_PATH = os.path.join(CP8_DIR, "logs", "run_18_cp8_case_studies.log")
TRACKS_DIR = os.path.join(CP8_DIR, "virtual4c_tracks")

MCOOL_PATH = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic2/phase-h3/607.mcool"
RESOLUTION = 10000

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def extract_virtual_4c(clr, chrom, start, end, tss, target_gene, out_dir):
    # Viewpoint is TSS +/- 10kb
    vp_start = max(0, tss - 10000)
    vp_end = tss + 10000
    
    # We want +/- 1Mb around TSS for the plot
    plot_start = max(0, tss - 1000000)
    plot_end = tss + 1000000
    
    region_str = f"{chrom}:{plot_start}-{plot_end}"
    vp_str = f"{chrom}:{vp_start}-{vp_end}"
    
    log(f"  Extracting Virtual 4C for {target_gene} (Viewpoint: {vp_str})")
    
    # Extract submatrix for plot region
    try:
        mat = clr.matrix(balance=True).fetch(region_str)
        bins = clr.bins().fetch(region_str)
    except Exception as e:
        log(f"  Failed to fetch region: {e}")
        return
    
    # Find bins corresponding to viewpoint
    vp_mask = (bins['start'] < vp_end) & (bins['end'] > vp_start)
    vp_indices = np.where(vp_mask)[0]
    
    if len(vp_indices) == 0:
        log(f"  No bins found for viewpoint {vp_str}.")
        return
        
    # Sum contacts over the viewpoint bins
    v4c_profile = np.nansum(mat[vp_indices, :], axis=0)
    
    # Write bedGraph
    bg_path = os.path.join(out_dir, f"{target_gene}_v4c_10kb.bedGraph")
    with open(bg_path, "w") as f:
        f.write(f"track type=bedGraph name=\"{target_gene} Virtual 4C\" description=\"Virtual 4C from {vp_str} at 10kb\"\n")
        for i, val in enumerate(v4c_profile):
            if np.isnan(val) or val == 0:
                continue
            bin_start = bins.iloc[i]['start']
            bin_end = bins.iloc[i]['end']
            f.write(f"{chrom}\t{bin_start}\t{bin_end}\t{val:.4f}\n")
            
    log(f"  Saved Virtual 4C track: {bg_path}")

def main():
    os.makedirs(CP8_DIR, exist_ok=True)
    os.makedirs(TRACKS_DIR, exist_ok=True)
    log("Started Phase H6 - Task 18 (CP8): Case Study Selection & Virtual 4C")
    
    scores_path = os.path.join(CP7_DIR, "candidate_gene_scores.tsv")
    if not os.path.exists(scores_path):
        log(f"Error: {scores_path} not found.")
        sys.exit(1)
        
    df = pd.read_csv(scores_path, sep='\t')
    log(f"Loaded {len(df)} scored interactions.")
    
    # 1. Select Top Candidates
    # Requirements: max score (5 or 6), no multiple plausible targets flag
    max_score = df['total_score'].max()
    log(f"Max observed score is {max_score}.")
    
    candidates = df[(df['total_score'] == max_score) & 
                    (df['flag_multiple_plausible_targets'] == False)]
                    
    # Drop duplicates by target gene to get unique loci
    candidates = candidates.drop_duplicates(subset=['target_gene']).copy()
    
    # Sort by some deterministic criteria, e.g., TSS
    candidates = candidates.sort_values('gene_tss')
    
    top_5 = candidates.head(5)
    log(f"Selected {len(top_5)} pristine loci for Case Studies.")
    
    out_table = os.path.join(CP8_DIR, "case_study_selection_table.tsv")
    top_5.to_csv(out_table, sep='\t', index=False)
    log(f"Saved selection table to {out_table}")
    
    # 2. Extract Virtual 4C
    cooler_uri = f"{MCOOL_PATH}::/resolutions/{RESOLUTION}"
    log(f"Opening Cooler: {cooler_uri}")
    try:
        clr = cooler.Cooler(cooler_uri)
    except Exception as e:
        log(f"Failed to open cooler: {e}")
        sys.exit(1)
        
    for _, row in top_5.iterrows():
        qtl_name = row['qtl_name']
        t_gene = row['target_gene'].replace('|', '_')
        c = row['gene_chrom']
        s = row['gene_tss']
        extract_virtual_4c(clr, c, s, s, s, t_gene, TRACKS_DIR)
        
    log("Task 18 completed successfully.")

if __name__ == "__main__":
    main()
