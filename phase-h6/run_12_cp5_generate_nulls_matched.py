#!/usr/bin/env python3
# run_12_cp5_generate_nulls_matched.py
# CP5: Generate Covariate-Matched Null Intervals (TSS density + Loop Anchor density)
import os
import sys
import datetime
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP1_DIR = os.path.join(PHASE_DIR, "CP1_loop_input")
CP2_DIR = os.path.join(PHASE_DIR, "CP2_QTL_Universe")
CP3_DIR = os.path.join(PHASE_DIR, "CP3_Positional_Baseline")
CP5_DIR = os.path.join(PHASE_DIR, "CP5_Null_Model")
LOG_PATH = os.path.join(CP5_DIR, "logs", "run_12_cp5_generate_nulls.log")

PERMUTATIONS = 1000
BIN_RES = 10000  # 10 kb bins for cumsum

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def build_density_map(df, chrom_col, pos_col):
    log(f"Building Density Cumulative Sum array for {pos_col}...")
    density_map = {}
    
    for chrom, grp in df.groupby(chrom_col):
        max_pos = grp[pos_col].max()
        if pd.isna(max_pos): continue
        
        num_bins = int(max_pos // BIN_RES) + 2
        counts = np.zeros(num_bins, dtype=np.int32)
        
        # bin indices
        idxs = (grp[pos_col] // BIN_RES).astype(int).values
        np.add.at(counts, idxs, 1)
        
        density_map[chrom] = {
            'cumsum': np.cumsum(counts),
            'max_pos': max_pos,
            'num_bins': num_bins
        }
    return density_map

def count_items(chrom, start, end, density_map):
    if chrom not in density_map: return 0
    c_map = density_map[chrom]
    
    start_bin = max(0, int(start // BIN_RES))
    end_bin = min(c_map['num_bins'] - 1, int(end // BIN_RES))
    
    if start_bin >= c_map['num_bins']: return 0
    
    count_end = c_map['cumsum'][end_bin]
    count_start = c_map['cumsum'][start_bin - 1] if start_bin > 0 else 0
    return count_end - count_start

def generate_nulls(q_chrom, q_start, q_end, tss_map, anchor_map, num_perms=1000):
    if q_chrom not in tss_map or q_chrom not in anchor_map:
        return []
        
    t_c_map = tss_map[q_chrom]
    max_pos = t_c_map['max_pos']
    q_size = q_end - q_start
    
    if q_size >= max_pos:
        return []
        
    orig_tss = count_items(q_chrom, q_start, q_end, tss_map)
    orig_anc = count_items(q_chrom, q_start, q_end, anchor_map)
    
    # Define tolerance: +/- 20% or +/- 1 item
    min_tss = max(0, min(orig_tss - 1, int(orig_tss * 0.8)))
    max_tss = max(1, max(orig_tss + 1, int(orig_tss * 1.2)))
    
    min_anc = max(0, min(orig_anc - 1, int(orig_anc * 0.8)))
    max_anc = max(1, max(orig_anc + 1, int(orig_anc * 1.2)))
    
    null_intervals = []
    
    chunk_size = 20000
    
    # Safety breakout to avoid infinite loops on very sparse chromosomes
    max_attempts = 1000000 
    attempts = 0
    
    while len(null_intervals) < num_perms and attempts < max_attempts:
        r_starts = np.random.randint(0, int(max_pos - q_size), size=chunk_size)
        r_ends = r_starts + q_size
        
        for s, e in zip(r_starts, r_ends):
            attempts += 1
            c_tss = count_items(q_chrom, s, e, tss_map)
            if min_tss <= c_tss <= max_tss:
                c_anc = count_items(q_chrom, s, e, anchor_map)
                if min_anc <= c_anc <= max_anc:
                    null_intervals.append((s, e, c_tss, c_anc))
                    if len(null_intervals) == num_perms:
                        break
                    
    return null_intervals

def main():
    os.makedirs(CP5_DIR, exist_ok=True)
    log("Started Phase H6 - Task 12 (CP5): Generate STRICT Covariate-Matched Null Intervals")
    
    # Load TSS
    tss_path = os.path.join(CP3_DIR, "reference", "rn7_protein_coding_tss.bed")
    df_tss = pd.read_csv(tss_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    df_tss['center'] = (df_tss['start'] + df_tss['end']) // 2
    
    # Load Loops for Anchor density
    loop_path = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085_with_roles.tsv")
    if not os.path.exists(loop_path):
        log(f"Error: {loop_path} not found.")
        sys.exit(1)
        
    df_loops = pd.read_csv(loop_path, sep='\t')
    
    # Extract all anchor centers
    a1 = df_loops[['chrom1', 'start1', 'end1']].copy()
    a1['chrom'] = a1['chrom1']
    a1['center'] = (a1['start1'] + a1['end1']) // 2
    
    a2 = df_loops[['chrom2', 'start2', 'end2']].copy()
    a2['chrom'] = a2['chrom2']
    a2['center'] = (a2['start2'] + a2['end2']) // 2
    
    df_anchors = pd.concat([a1[['chrom', 'center']], a2[['chrom', 'center']]], ignore_index=True)
    
    tss_map = build_density_map(df_tss, 'chrom', 'center')
    anchor_map = build_density_map(df_anchors, 'chrom', 'center')
    
    primary_bed = os.path.join(CP2_DIR, "qtl_filtered_primary.bed")
    df_qtl = pd.read_csv(primary_bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    log(f"Processing {len(df_qtl)} Primary QTLs for {PERMUTATIONS} permutations each...")
    
    all_nulls = []
    
    for i, row in df_qtl.iterrows():
        nulls = generate_nulls(row['chrom'], row['start'], row['end'], tss_map, anchor_map, num_perms=PERMUTATIONS)
        for perm_id, (ns, ne, nc_tss, nc_anc) in enumerate(nulls):
            all_nulls.append({
                'qtl_name': row['name'],
                'chrom': row['chrom'],
                'perm_id': perm_id + 1,
                'null_start': ns,
                'null_end': ne,
                'null_tss_count': nc_tss,
                'null_anchor_count': nc_anc
            })
            
        if (i+1) % 50 == 0:
            log(f"  Processed {i+1} / {len(df_qtl)} QTLs.")
            
    df_nulls = pd.DataFrame(all_nulls)
    
    out_path = os.path.join(CP5_DIR, "null_intervals_density_matched.tsv")
    df_nulls.to_csv(out_path, sep='\t', index=False)
    
    log(f"Successfully generated {len(df_nulls)} STRICT null intervals.")
    log(f"Saved to {out_path}")
    log("Task 12 completed successfully.")

if __name__ == "__main__":
    main()
