#!/usr/bin/env python3
# run_17_cp7_candidate_scoring.py
# CP7: Candidate Gene Scoring & Prioritization

import os
import sys
import datetime
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP4_DIR = os.path.join(PHASE_DIR, "CP4_3D_Mapping")
CP6_DIR = os.path.join(PHASE_DIR, "CP6_Support_Evidence")
CP7_DIR = os.path.join(PHASE_DIR, "CP7_Scoring")
LOG_PATH = os.path.join(CP7_DIR, "logs", "run_17_cp7_scoring.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def main():
    os.makedirs(CP7_DIR, exist_ok=True)
    log("Started Phase H6 - Task 17 (CP7): Candidate Gene Scoring")
    
    # 1. Load Primary Targets
    targets_path = os.path.join(CP4_DIR, "qtl_loop_target_genes_primary.tsv")
    df = pd.read_csv(targets_path, sep='\t')
    log(f"Loaded {len(df)} primary interactions.")
    
    # Base Score: 3D Mapping (+2)
    df['score_3d_base'] = 2
    
    # APA Score
    # We remove the blanket +1 score because an aggregate APA score of 1.174~1.407 
    # does not imply that every individual loop is significantly APA-supported.
    df['score_apa'] = 0
    
    # 2. Load 1D Baseline comparison (Novel distal)
    baseline_path = os.path.join(PHASE_DIR, "CP3_Positional_Baseline", "nearest_tss_primary.tsv")
    if os.path.exists(baseline_path):
        df_nearest = pd.read_csv(baseline_path, sep='\t')
        # Create a set of (qtl_name, target_gene) that are nearest
        nearest_set = set(zip(df_nearest['name'], df_nearest['nearest_gene']))
        
        # +1 if target is NOT the nearest gene
        df['score_novel_distal'] = df.apply(
            lambda r: 0 if (r['qtl_name'], r['target_gene']) in nearest_set else 1, axis=1
        )
    else:
        df['score_novel_distal'] = 0
        
    # 3. Load Epigenomic ATAC support
    atac_path = os.path.join(CP6_DIR, "candidate_epigenomic_support.tsv")
    if os.path.exists(atac_path):
        df_atac = pd.read_csv(atac_path, sep='\t')
        if 'qtl_id' in df_atac.columns: df_atac = df_atac.rename(columns={'qtl_id': 'qtl_name'})
        atac_dict = df_atac.set_index(['qtl_name', 'target_gene'])['atac_support'].to_dict()
        df['score_atac'] = df.apply(lambda r: atac_dict.get((r['qtl_name'], r['target_gene']), 0), axis=1)
    else:
        df['score_atac'] = 0

    # 4. Load Ortholog support
    ortho_path = os.path.join(CP6_DIR, "candidate_ortholog_support.tsv")
    if os.path.exists(ortho_path):
        df_ortho = pd.read_csv(ortho_path, sep='\t')
        if 'qtl_id' in df_ortho.columns: df_ortho = df_ortho.rename(columns={'qtl_id': 'qtl_name'})
        # +1 if human ortholog exists
        ortho_dict = df_ortho.set_index(['qtl_name', 'target_gene'])['human_ortholog'].notna().astype(int).to_dict()
        df['score_ortholog'] = df.apply(lambda r: ortho_dict.get((r['qtl_name'], r['target_gene']), 0), axis=1)
    else:
        df['score_ortholog'] = 0
        
    # Total Score
    score_cols = ['score_3d_base', 'score_apa', 'score_novel_distal', 'score_atac', 'score_ortholog']
    df['total_score'] = df[score_cols].sum(axis=1)
    
    # Tiers
    # Tier 1: >= 5
    # Tier 2: 3-4
    # Exploratory: 1-2
    def get_tier(s):
        if s >= 5: return 'Tier 1'
        if s >= 3: return 'Tier 2'
        return 'Exploratory'
    df['tier'] = df['total_score'].apply(get_tier)
    
    # Flags
    df['flag_public_epigenomics_not_matched'] = (df['score_atac'] == 0)
    
    # Check for multiple plausible targets per QTL
    qtl_counts = df.groupby('qtl_name').size()
    multiple_qtls = qtl_counts[qtl_counts >= 3].index
    df['flag_multiple_plausible_targets'] = df['qtl_name'].isin(multiple_qtls)
    
    df['flag_no_cortex_expression_support'] = True # Skipping eQTL
    df['flag_human_only_support'] = (df['score_ortholog'] == 1) & (df['score_atac'] == 0)
    
    log("Score Summary:")
    tier_counts = df['tier'].value_counts()
    for t, c in tier_counts.items():
        log(f"  {t}: {c} interactions")
        
    out_scores = os.path.join(CP7_DIR, "candidate_gene_scores.tsv")
    df.to_csv(out_scores, sep='\t', index=False)
    log(f"Saved full scoring to {out_scores}")
    
    tier1 = df[df['tier'] == 'Tier 1']
    out_tier1 = os.path.join(CP7_DIR, "tier1_candidates.tsv")
    tier1.to_csv(out_tier1, sep='\t', index=False)
    log(f"Saved {len(tier1)} Tier 1 candidates to {out_tier1}")
    
    log("Task 17 completed successfully.")

if __name__ == "__main__":
    main()
