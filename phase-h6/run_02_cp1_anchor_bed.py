#!/usr/bin/env python3
# run_02_cp1_anchor_bed.py
import os
import sys
import datetime
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP1_DIR = os.path.join(PHASE_DIR, "CP1_loop_input")

BEDPE_IN = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085.tsv")
ANCHOR_DIR = os.path.join(CP1_DIR, "loop_anchor_bed")
LOG_PATH = os.path.join(CP1_DIR, "logs", "run_02_cp1_anchor_bed.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def main():
    log("Started Phase H6 - Task 2: Anchor BED Extraction")
    
    if not os.path.exists(BEDPE_IN):
        log(f"Error: Required input {BEDPE_IN} not found. Run Task 1 first.")
        sys.exit(1)

    df_bedpe = pd.read_csv(BEDPE_IN, sep='\t')
    log(f"Loaded {len(df_bedpe)} loops from BEDPE.")

    # Extract all anchors
    anchors = []
    for _, row in df_bedpe.iterrows():
        # Left anchor
        anchors.append({
            'chrom': row['chrom1'],
            'start': int(row['start1']),
            'end': int(row['end1']),
            'anchor_id': f"{row['loop_id']}_left",
            'score': '.',
            'strand': '.',
            'loop_id': row['loop_id'],
            'type': 'unknown' # Placeholder for CP4 refinement
        })
        # Right anchor
        anchors.append({
            'chrom': row['chrom2'],
            'start': int(row['start2']),
            'end': int(row['end2']),
            'anchor_id': f"{row['loop_id']}_right",
            'score': '.',
            'strand': '.',
            'loop_id': row['loop_id'],
            'type': 'unknown'
        })

    df_anchors = pd.DataFrame(anchors)
    
    # Save all unique anchors
    df_anchors_unique = df_anchors.drop_duplicates(subset=['chrom', 'start', 'end'])
    log(f"Extracted {len(df_anchors_unique)} unique loop anchors.")

    # 1. all_loop_anchors.bed
    all_bed_path = os.path.join(ANCHOR_DIR, "all_loop_anchors.bed")
    df_anchors_unique[['chrom', 'start', 'end', 'anchor_id', 'score', 'strand', 'type']].to_csv(
        all_bed_path, sep='\t', index=False, header=False
    )
    log(f"Saved {all_bed_path}")

    # 2. anchor_midpoints_1bp.bed
    df_midpoints = df_anchors_unique.copy()
    df_midpoints['midpoint'] = ((df_midpoints['start'] + df_midpoints['end']) / 2).astype(int)
    df_midpoints['start'] = df_midpoints['midpoint']
    df_midpoints['end'] = df_midpoints['midpoint'] + 1
    mid_bed_path = os.path.join(ANCHOR_DIR, "anchor_midpoints_1bp.bed")
    df_midpoints[['chrom', 'start', 'end', 'anchor_id', 'score', 'strand', 'type']].to_csv(
        mid_bed_path, sep='\t', index=False, header=False
    )
    log(f"Saved {mid_bed_path}")

    # 3. anchor_windows_5kb.bed (±2.5kb around midpoint → 5kb window)
    df_win = df_anchors_unique.copy()
    df_win['midpoint'] = ((df_win['start'] + df_win['end']) / 2).astype(int)
    df_win['start'] = (df_win['midpoint'] - 2500).clip(lower=0)
    df_win['end'] = df_win['midpoint'] + 2500
    win_bed_path = os.path.join(ANCHOR_DIR, "anchor_windows_5kb.bed")
    df_win[['chrom', 'start', 'end', 'anchor_id', 'score', 'strand', 'type']].to_csv(
        win_bed_path, sep='\t', index=False, header=False
    )
    log(f"Saved {win_bed_path}")

    log("Task 2 completed successfully.")

if __name__ == "__main__":
    main()
