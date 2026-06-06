#!/usr/bin/env python3
# run_13_cp5_permutation_test.py
# CP5: Compute Empirical P-values by intersecting 263,000 null intervals with 3D Loops.
import os
import sys
import datetime
import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP1_DIR = os.path.join(PHASE_DIR, "CP1_loop_input")
CP3_DIR = os.path.join(PHASE_DIR, "CP3_Positional_Baseline")
CP4_DIR = os.path.join(PHASE_DIR, "CP4_3D_Mapping")
CP5_DIR = os.path.join(PHASE_DIR, "CP5_Null_Model")
LOG_PATH = os.path.join(CP5_DIR, "logs", "run_13_cp5_permutation_test.log")

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
    os.makedirs(CP5_DIR, exist_ok=True)
    log("Started Phase H6 - Task 13 (CP5): Empirical P-value Computation")
    
    # 1. Load Real Data Counts
    path_real = os.path.join(CP4_DIR, "qtl_3d_vs_baseline_comparison_primary.tsv")
    df_real = pd.read_csv(path_real, sep='\t')
    real_counts = df_real.set_index('qtl_name')['num_3d_targets'].to_dict()
    log(f"Loaded real target counts for {len(real_counts)} primary QTLs.")
    
    # 2. Load Null Intervals
    nulls_path = os.path.join(CP5_DIR, "null_intervals_density_matched.tsv")
    df_nulls = pd.read_csv(nulls_path, sep='\t')
    log(f"Loaded {len(df_nulls)} null intervals.")
    
    # Prefix columns for get_overlaps
    df_nulls_fmt = df_nulls.rename(columns={'null_start': 'n_start', 'null_end': 'n_end'})
    df_nulls_fmt['n_chrom'] = df_nulls_fmt['chrom']
    
    # 3. Load Loops & TSS
    loop_path = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085_with_roles.tsv")
    df_loops = pd.read_csv(loop_path, sep='\t')
    df_loops = df_loops[df_loops['interaction_type'].isin(['P-E', 'P-P'])].copy()
    df_loops.reset_index(drop=True, inplace=True)
    
    df_a1 = df_loops[['chrom1', 'start1', 'end1', 'loop_id']].rename(columns={'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'})
    df_a2 = df_loops[['chrom2', 'start2', 'end2', 'loop_id']].rename(columns={'chrom2': 'chrom', 'start2': 'start', 'end2': 'end'})
    
    tss_path = os.path.join(CP3_DIR, "reference", "rn7_protein_coding_tss.bed")
    df_tss = pd.read_csv(tss_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    df_tss_renamed = df_tss.rename(columns={'chrom': 't_chrom', 'start': 't_start', 'end': 't_end', 'name': 't_name'})
    # Expand TSS by 200kb for null mapping as well (to match CP4 logic)
    df_tss_renamed['t_start'] = np.maximum(0, df_tss_renamed['t_start'] - 200000)
    df_tss_renamed['t_end'] = df_tss_renamed['t_end'] + 200000
    
    # 4. Map TSS to Anchors (Reuse from CP4)
    log("Mapping TSS to Loop Anchors...")
    tss_a1 = get_overlaps(df_tss_renamed, df_a1, 't_', '')
    tss_a2 = get_overlaps(df_tss_renamed, df_a2, 't_', '')
    
    a1_to_tss = {}
    for t_idx, l_idx in tss_a1:
        a1_to_tss.setdefault(l_idx, []).append(t_idx)
        
    a2_to_tss = {}
    for t_idx, l_idx in tss_a2:
        a2_to_tss.setdefault(l_idx, []).append(t_idx)
        
    # 5. Map Nulls to Anchors
    log("Intersecting 263,000 Null Intervals with Loop Anchors (this may take a minute)...")
    null_a1 = get_overlaps(df_nulls_fmt, df_a1, 'n_', '')
    null_a2 = get_overlaps(df_nulls_fmt, df_a2, 'n_', '')
    
    # 6. Accumulate Null Target Counts
    log("Calculating Target Genes per Null Interval...")
    # Dictionary to hold count of unique genes for each null (idx)
    null_genes = {i: set() for i in range(len(df_nulls))}
    
    # Directional mapping logic for nulls
    for n_idx, l_idx in null_a1:
        loop_row = df_loops.iloc[l_idx]
        if loop_row['anchor2_role'] == 'Promoter_Anchor':
            if l_idx in a2_to_tss:
                for t_idx in a2_to_tss[l_idx]:
                    null_genes[n_idx].add(t_idx)
                
    for n_idx, l_idx in null_a2:
        loop_row = df_loops.iloc[l_idx]
        if loop_row['anchor1_role'] == 'Promoter_Anchor':
            if l_idx in a1_to_tss:
                for t_idx in a1_to_tss[l_idx]:
                    null_genes[n_idx].add(t_idx)
                
    # Convert sets to lengths
    df_nulls['num_3d_targets_null'] = [len(null_genes[i]) for i in range(len(df_nulls))]
    
    # 7. Compute P-values
    log("Computing Empirical P-values...")
    results = []
    
    for qtl_name, group in df_nulls.groupby('qtl_name'):
        real_count = real_counts.get(qtl_name, 0)
        
        null_counts = group['num_3d_targets_null'].values
        num_perms = len(null_counts)
        
        # Empirical p-value = (N_better_or_equal + 1) / (N_perms + 1)
        # We test if the real QTL has significantly MORE 3D targets than chance
        # So we count how many nulls have >= real_count
        better_or_equal = np.sum(null_counts >= real_count)
        
        p_val = (better_or_equal + 1) / (num_perms + 1)
        
        # Odds ratio / Enrichment Score approximation
        # (Real Count + 1) / (Mean Null Count + 1)
        mean_null = np.mean(null_counts)
        enrichment = (real_count + 1) / (mean_null + 1)
        
        results.append({
            'qtl_name': qtl_name,
            'real_target_count': real_count,
            'mean_null_count': mean_null,
            'enrichment_score': enrichment,
            'empirical_pval': p_val
        })
        
    df_res = pd.DataFrame(results)
    
    # Compute FDR (Benjamini-Hochberg)
    _, fdr, _, _ = smm.multipletests(df_res['empirical_pval'], method='fdr_bh')
    df_res['fdr'] = fdr
    
    out_path = os.path.join(CP5_DIR, "permutation_results_primary.tsv")
    df_res.to_csv(out_path, sep='\t', index=False)
    
    sig_count = (df_res['fdr'] < 0.05).sum()
    log(f"Analyzed {len(df_res)} QTLs.")
    log(f"Found {sig_count} QTLs with FDR < 0.05 for 3D Target Enrichment.")
    log(f"Saved to {out_path}")
    log("Task 13 completed successfully.")

if __name__ == "__main__":
    main()
