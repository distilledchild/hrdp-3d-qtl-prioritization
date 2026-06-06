#!/usr/bin/env python3
# run_03b_apa_resolution_matched.py
# Resolution-matched APA: pileup each loop subset at its calling resolution
import os
import sys
import glob
import datetime
import subprocess
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_H3_DIR = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic2/phase-h3"
CP1_DIR = os.path.join(REPO_ROOT, "phase-h6", "CP1_loop_input")

BEDPE_IN = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085.tsv")
APA_DIR = os.path.join(CP1_DIR, "apa")
LOG_PATH = os.path.join(CP1_DIR, "logs", "run_03b_apa_resolution_matched.log")

RESOLUTIONS = [5000, 10000, 25000]
FLANK_MULTIPLIER = 10  # flank = resolution × 10

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def compute_apa_score(npz_path):
    try:
        data = np.load(npz_path, allow_pickle=True)
        mtx = data['S'] if 'S' in data else data[data.files[0]]
        n = mtx.shape[0]
        c = n // 2

        # Center 3×3
        center = np.nanmean(mtx[c-1:c+2, c-1:c+2])

        # Corner 5×5 (all four)
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
    log("Started Resolution-Matched APA Validation")

    if not os.path.exists(BEDPE_IN):
        log(f"Error: {BEDPE_IN} not found.")
        sys.exit(1)

    df_all = pd.read_csv(BEDPE_IN, sep='\t')
    log(f"Loaded {len(df_all)} total loops.")

    mcool_files = sorted(glob.glob(os.path.join(PHASE_H3_DIR, "*.mcool")))
    if not mcool_files:
        log("Error: No .mcool files found.")
        sys.exit(1)
    log(f"Found {len(mcool_files)} mcool files.")

    os.makedirs(APA_DIR, exist_ok=True)
    all_scores = []

    for res in RESOLUTIONS:
        df_sub = df_all[df_all['resolution'] == res].copy()
        if len(df_sub) == 0:
            log(f"No loops at resolution {res}. Skipping.")
            continue

        # Write temp BEDPE for this resolution subset
        sub_bedpe = os.path.join(APA_DIR, f"loops_res_{res}.tsv")
        df_sub.to_csv(sub_bedpe, sep='\t', index=False)
        log(f"Resolution {res}: {len(df_sub)} loops")

        flank = res * FLANK_MULTIPLIER

        for mcool in mcool_files:
            sample = os.path.basename(mcool).split('.')[0]
            out_prefix = os.path.join(APA_DIR, f"{sample}_apa_matched_{res}")
            out_npz = f"{out_prefix}.npz"

            if not os.path.exists(out_npz):
                cmd = [
                    "uv", "run", "--python", "3.10",
                    "--with", "cooltools==0.7.1",
                    "--with", "pandas<2.2.0",
                    "--with", "numpy<2.0.0",
                    "cooltools", "pileup",
                    f"{mcool}::/resolutions/{res}",
                    sub_bedpe,
                    "--features-format", "bedpe",
                    "--flank", str(flank),
                    "--out", out_prefix,
                    "-p", "4"
                ]
                log(f"Running pileup: {sample} @ {res}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    log(f"FAILED: {sample} @ {res}: {result.stderr[:200]}")
                    continue
            else:
                log(f"Already exists: {sample} @ {res}")

            score = compute_apa_score(out_npz)
            log(f"{sample} @ {res}: APA = {score:.3f}")
            all_scores.append({
                'sample': sample,
                'resolution': res,
                'n_loops': len(df_sub),
                'apa_score': round(score, 3)
            })

    # Save combined results
    df_scores = pd.DataFrame(all_scores)
    scores_path = os.path.join(APA_DIR, "apa_scores_resolution_matched.tsv")
    df_scores.to_csv(scores_path, sep='\t', index=False)
    log(f"Saved to {scores_path}")

    # Print summary pivot
    if not df_scores.empty:
        pivot = df_scores.pivot_table(index='sample', columns='resolution', values='apa_score')
        log(f"\n=== APA Score Summary (Resolution-Matched) ===\n{pivot.to_string()}")

    log("Resolution-Matched APA completed.")

if __name__ == "__main__":
    main()
