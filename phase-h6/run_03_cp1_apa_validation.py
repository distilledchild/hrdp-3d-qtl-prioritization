#!/usr/bin/env python3
# run_03_cp1_apa_validation.py
import os
import sys
import glob
import datetime
import subprocess
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# Phase H3 matrices in Dropbox
PHASE_H3_DIR = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic2/phase-h3"
PHASE_H6_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP1_DIR = os.path.join(PHASE_H6_DIR, "CP1_loop_input")

BEDPE_IN = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085.tsv")
APA_DIR = os.path.join(CP1_DIR, "apa")
LOG_PATH = os.path.join(CP1_DIR, "logs", "run_03_cp1_apa_validation.log")

RESOLUTION = 25000
FLANK = 250000

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def run_cmd(cmd):
    log(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"Command failed with exit code {result.returncode}")
        log(f"STDERR: {result.stderr}")
    return result.returncode == 0

def compute_apa_score(npz_path):
    """
    Calculate a simple APA enrichment score (center / corners).
    """
    try:
        data = np.load(npz_path, allow_pickle=True)
        # 'S' is usually the sum, 'count' is the number of valid pixels
        if 'S' in data:
            mtx = data['S']
        else:
            mtx = data[data.files[0]]
            
        n = mtx.shape[0]
        c = n // 2
        
        # Center 3x3
        center = np.nanmean(mtx[c-1:c+2, c-1:c+2])
        
        # Corners 5x5
        corner_tl = np.nanmean(mtx[0:5, 0:5])
        corner_tr = np.nanmean(mtx[0:5, -5:])
        corner_bl = np.nanmean(mtx[-5:, 0:5])
        corner_br = np.nanmean(mtx[-5:, -5:])
        corners = np.nanmean([corner_tl, corner_tr, corner_bl, corner_br])
        
        if corners == 0 or np.isnan(corners) or np.isnan(center):
            return 0.0
        
        return float(center / corners)
    except Exception as e:
        log(f"Failed to calculate APA score for {npz_path}: {e}")
        return 0.0

def main():
    log("Started Phase H6 - Task 3: APA Validation")
    
    if not os.path.exists(BEDPE_IN):
        log(f"Error: {BEDPE_IN} not found. Run Task 1 first.")
        sys.exit(1)
    
    # Validate BEDPE column names match cooltools expectations
    df_check = pd.read_csv(BEDPE_IN, sep='\t', nrows=2)
    required = {'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'}
    if not required.issubset(set(df_check.columns)):
        log(f"Error: BEDPE missing required columns. Found: {list(df_check.columns)}")
        sys.exit(1)
    log(f"BEDPE columns validated: {list(df_check.columns)}")
        
    mcool_files = glob.glob(os.path.join(PHASE_H3_DIR, "*.mcool"))
    if not mcool_files:
        log("Error: No .mcool files found in phase-h3 directory. APA requires cool files.")
        sys.exit(1)
        
    log(f"Found {len(mcool_files)} mcool files.")
    
    # Run cooltools pileup for each sample
    scores = []
    
    for mcool in mcool_files:
        sample = os.path.basename(mcool).split('.')[0]
        out_prefix = os.path.join(APA_DIR, f"{sample}_apa_{RESOLUTION}")
        # cooltools pileup auto-appends .npz to --out prefix
        out_npz = f"{out_prefix}.npz"
        
        # cooltools pileup command (quotes around version constraints for shell safety)
        cmd = [
            "uv", "run", "--python", "3.10",
            "--with", "cooltools==0.7.1",
            "--with", "pandas<2.2.0",
            "--with", "numpy<2.0.0",
            "cooltools", "pileup",
            f"{mcool}::/resolutions/{RESOLUTION}",
            BEDPE_IN,
            "--features-format", "bedpe",
            "--flank", str(FLANK),
            "--out", out_prefix,
            "-p", "4"
        ]
        
        if not os.path.exists(out_npz):
            success = run_cmd(cmd)
            if not success:
                continue
        else:
            log(f"Pileup already exists for {sample}, skipping calculation.")
            
        # Compute Score
        score = compute_apa_score(out_npz)
        log(f"Sample: {sample} | APA Score: {score:.3f}")
        scores.append({'sample': sample, 'apa_score': score})
        
    # Save scores
    df_scores = pd.DataFrame(scores)
    scores_path = os.path.join(APA_DIR, "apa_scores_per_sample.tsv")
    df_scores.to_csv(scores_path, sep='\t', index=False)
    log(f"Saved APA scores to {scores_path}")
    
    # We will pool the scores (average) to determine which loops are supported.
    # Note: A real APA score per-loop requires `--by-strand` or `--by-feature` which is slow.
    # For now, we calculate global APA per sample as required by CP1.
    
    log("Task 3 completed successfully.")

if __name__ == "__main__":
    main()
