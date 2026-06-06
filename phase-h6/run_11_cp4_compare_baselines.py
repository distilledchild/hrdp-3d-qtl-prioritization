#!/usr/bin/env python3
# run_11_cp4_compare_baselines.py
# CP4: Compare 3D-assigned target genes with 1D positional baselines to identify novel long-range targets.
import os
import sys
import datetime
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP3_DIR = os.path.join(PHASE_DIR, "CP3_Positional_Baseline")
CP4_DIR = os.path.join(PHASE_DIR, "CP4_3D_Mapping")
LOG_PATH = os.path.join(CP4_DIR, "logs", "run_11_cp4_compare_baselines.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def run_comparison(prefix):
    log(f"Running comparison for {prefix} QTLs...")
    
    # 1. Load 3D Targets
    path_3d = os.path.join(CP4_DIR, f"qtl_loop_target_genes_{prefix}.tsv")
    if not os.path.exists(path_3d):
        log(f"Warning: 3D targets not found for {prefix}")
        return
        
    df_3d = pd.read_csv(path_3d, sep='\t')
    # Group by QTL -> Set of genes
    # df_3d has 'qtl_name', 'target_gene'
    target_3d = df_3d.groupby('qtl_name')['target_gene'].apply(lambda x: set(x.dropna())).to_dict()
    
    # 2. Load Interval Targets
    path_interval = os.path.join(CP3_DIR, f"interval_genes_{prefix}.tsv")
    target_interval = {}
    if os.path.exists(path_interval):
        df_int = pd.read_csv(path_interval, sep='\t')
        target_interval = df_int.groupby('qtl_name')['gene_name'].apply(lambda x: set(x.dropna())).to_dict()
        
    # 3. Load Nearest TSS
    path_nearest = os.path.join(CP3_DIR, f"nearest_tss_{prefix}.tsv")
    target_nearest = {}
    if os.path.exists(path_nearest):
        df_near = pd.read_csv(path_nearest, sep='\t')
        # df_near has 'name' as qtl_name and 'nearest_gene'
        for _, row in df_near.iterrows():
            target_nearest[row['name']] = str(row['nearest_gene'])
            
    # 4. Compare
    results = []
    
    all_qtls = set(target_3d.keys()).union(set(target_interval.keys())).union(set(target_nearest.keys()))
    
    for q in sorted(list(all_qtls)):
        genes_3d = target_3d.get(q, set())
        genes_int = target_interval.get(q, set())
        near_gene = target_nearest.get(q, "None")
        
        # Novel = in 3D but not in 1D interval
        genes_novel = genes_3d - genes_int
        
        nearest_in_3d = near_gene in genes_3d
        nearest_in_int = near_gene in genes_int
        
        res = {
            'qtl_name': q,
            'num_3d_targets': len(genes_3d),
            'num_interval_targets': len(genes_int),
            'num_novel_targets': len(genes_novel),
            'novelty_ratio': len(genes_novel) / len(genes_3d) if len(genes_3d) > 0 else 0,
            'is_nearest_tss_in_3d': nearest_in_3d,
            'is_nearest_tss_in_interval': nearest_in_int
        }
        results.append(res)
        
    df_res = pd.DataFrame(results)
    
    out_path = os.path.join(CP4_DIR, f"qtl_3d_vs_baseline_comparison_{prefix}.tsv")
    df_res.to_csv(out_path, sep='\t', index=False)
    
    log(f"  Processed {len(df_res)} QTLs.")
    log(f"  Average Novelty Ratio: {df_res['novelty_ratio'].mean():.3f}")
    log(f"  QTLs with Nearest TSS captured by 3D loop: {df_res['is_nearest_tss_in_3d'].sum()}")
    log(f"  Saved {out_path}")

def main():
    log("Started Phase H6 - Task 11 (CP4): Compare 3D Assignments vs Baselines")
    
    run_comparison("primary")
    run_comparison("secondary")
    
    log("Task 11 completed successfully.")

if __name__ == "__main__":
    main()
