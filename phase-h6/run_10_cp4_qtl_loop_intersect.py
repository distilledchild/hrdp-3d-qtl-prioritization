#!/usr/bin/env python3
# run_10_cp4_qtl_loop_intersect.py
# CP4: Map QTLs to target genes via 3D loops (Directional Mapping)
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
CP4_DIR = os.path.join(PHASE_DIR, "CP4_3D_Mapping")
LOG_PATH = os.path.join(CP4_DIR, "logs", "run_10_cp4_qtl_loop_intersect.log")

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

def run_3d_mapping(qtl_bed_path, df_loops, out_path):
    if not os.path.exists(qtl_bed_path):
        log(f"Warning: {qtl_bed_path} not found.")
        return
        
    df_qtl = pd.read_csv(qtl_bed_path, sep='\t', header=None, names=['q_chrom', 'q_start', 'q_end', 'q_name', 'q_score', 'q_strand'])
    log(f"Loaded {len(df_qtl)} QTLs from {os.path.basename(qtl_bed_path)}")
    
    df_a1 = df_loops[['chrom1', 'start1', 'end1', 'loop_id']].rename(columns={'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'})
    df_a2 = df_loops[['chrom2', 'start2', 'end2', 'loop_id']].rename(columns={'chrom2': 'chrom', 'start2': 'start', 'end2': 'end'})
    
    log("  Intersecting QTLs with Anchor 1...")
    qtl_a1 = get_overlaps(df_qtl, df_a1, 'q_', '')
    
    log("  Intersecting QTLs with Anchor 2...")
    qtl_a2 = get_overlaps(df_qtl, df_a2, 'q_', '')
    
    results = []
    
    # Case A: QTL overlaps Anchor 1
    # For a valid directional mapping, Anchor 2 must contain Target Genes (Promoter role).
    # This includes P-E (if Anchor 2 is Promoter) and P-P interactions.
    for q_idx, l_idx in qtl_a1:
        loop_row = df_loops.iloc[l_idx]
        if loop_row['interaction_type'] in ['P-E', 'P-P'] and loop_row['anchor2_role'] == 'Promoter_Anchor':
            # Target is Anchor 2
            genes = str(loop_row['anchor2_genes']).split(',')
            for g in genes:
                if g.strip():
                    results.append({
                        'qtl_name': df_qtl.iloc[q_idx]['q_name'],
                        'loop_id': loop_row['loop_id'],
                        'qtl_hit_anchor': 'Anchor1',
                        'qtl_hit_role': loop_row['anchor1_role'],
                        'target_gene': g.strip(),
                        'anchor1_chrom': loop_row['chrom1'],
                        'anchor1_start': loop_row['start1'],
                        'anchor1_end': loop_row['end1'],
                        'anchor2_chrom': loop_row['chrom2'],
                        'anchor2_start': loop_row['start2'],
                        'anchor2_end': loop_row['end2']
                    })
                    
    # Case B: QTL overlaps Anchor 2
    for q_idx, l_idx in qtl_a2:
        loop_row = df_loops.iloc[l_idx]
        if loop_row['interaction_type'] in ['P-E', 'P-P'] and loop_row['anchor1_role'] == 'Promoter_Anchor':
            # Target is Anchor 1
            genes = str(loop_row['anchor1_genes']).split(',')
            for g in genes:
                if g.strip():
                    results.append({
                        'qtl_name': df_qtl.iloc[q_idx]['q_name'],
                        'loop_id': loop_row['loop_id'],
                        'qtl_hit_anchor': 'Anchor2',
                        'qtl_hit_role': loop_row['anchor2_role'],
                        'target_gene': g.strip(),
                        'anchor1_chrom': loop_row['chrom1'],
                        'anchor1_start': loop_row['start1'],
                        'anchor1_end': loop_row['end1'],
                        'anchor2_chrom': loop_row['chrom2'],
                        'anchor2_start': loop_row['start2'],
                        'anchor2_end': loop_row['end2']
                    })
                    
    df_res = pd.DataFrame(results)
    
    # Optional: Deduplicate
    df_res.drop_duplicates(inplace=True)
    
    df_res.to_csv(out_path, sep='\t', index=False)
    
    qtl_count = df_res['qtl_name'].nunique() if not df_res.empty else 0
    gene_count = df_res['target_gene'].nunique() if not df_res.empty else 0
    
    log(f"  Found 3D targets for {qtl_count} QTLs.")
    log(f"  Mapped to {gene_count} unique genes across {len(df_res)} True Directional interactions.")
    log(f"  Saved {out_path}")

def main():
    os.makedirs(CP4_DIR, exist_ok=True)
    log("Started Phase H6 - Task 10 (CP4): 3D QTL-to-Gene Directional Mapping")
    
    # Load Loops with Roles
    loop_path = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085_with_roles.tsv")
    if not os.path.exists(loop_path):
        log(f"Error: Loops with roles file not found at {loop_path}")
        sys.exit(1)
    df_loops = pd.read_csv(loop_path, sep='\t')
    # Filter only P-E and P-P loops (exclude E-E)
    df_loops = df_loops[df_loops['interaction_type'].isin(['P-E', 'P-P'])].copy()
    df_loops.reset_index(drop=True, inplace=True)
    log(f"Loaded {len(df_loops)} P-E / P-P Loops for directional mapping.")
    
    primary_bed = os.path.join(CP2_DIR, "qtl_filtered_primary.bed")
    primary_out = os.path.join(CP4_DIR, "qtl_loop_target_genes_primary.tsv")
    run_3d_mapping(primary_bed, df_loops, primary_out)
    
    secondary_bed = os.path.join(CP2_DIR, "qtl_broad_secondary.bed")
    secondary_out = os.path.join(CP4_DIR, "qtl_loop_target_genes_secondary.tsv")
    run_3d_mapping(secondary_bed, df_loops, secondary_out)
    
    log("Task 10 completed successfully.")

if __name__ == "__main__":
    main()
