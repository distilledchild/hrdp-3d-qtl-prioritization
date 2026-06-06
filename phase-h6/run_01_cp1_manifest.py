#!/usr/bin/env python3
# run_01_cp1_manifest.py
import os
import sys
import datetime
import pandas as pd
import numpy as np

# Define paths
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP1_DIR = os.path.join(PHASE_DIR, "CP1_loop_input")

# Input files
LOOP_SOURCE_CSV = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/r_files/figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv"
TSS_BED = os.path.join(PHASE_DIR, "CP3_Positional_Baseline", "reference", "rn7_protein_coding_tss.bed")

# Output directories
DIRS = [
    os.path.join(CP1_DIR, "metadata"),
    os.path.join(CP1_DIR, "paper1_input"),
    os.path.join(CP1_DIR, "apa"),
    os.path.join(CP1_DIR, "loop_anchor_bed"),
    os.path.join(CP1_DIR, "logs")
]

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    with open(os.path.join(CP1_DIR, "logs", "run_01_cp1_manifest.log"), "a") as f:
        f.write(full_msg + "\n")

def main():
    # 1. Create Directories
    for d in DIRS:
        os.makedirs(d, exist_ok=True)
    
    log("Started Phase H6 - Task 1: CP1 Manifest Generation (Refactored)")
    
    if not os.path.exists(TSS_BED):
        log(f"Error: TSS bed not found: {TSS_BED}")
        sys.exit(1)
        
    # Load TSS into pybedtools and pad by 200kb
    # Pybedtools window intersect will be used directly
    
    # 2. Parse Paper 1 Final PE-loops
    log(f"Reading source PE-loop table: {LOOP_SOURCE_CSV}")
    try:
        df_source = pd.read_csv(LOOP_SOURCE_CSV)
    except Exception as e:
        log(f"Failed to read source CSV: {e}")
        sys.exit(1)

    log(f"Found {len(df_source)} loops.")
    
    rows = []
    for _, row in df_source.iterrows():
        loop_id = str(row['loop.id'])
        parts = loop_id.split('_')
        if len(parts) == 7:
            chr1, x1, x2, chr2, y1, y2, res = parts
        else:
            log(f"Warning: Unexpected loop.id format: {loop_id}")
            continue
            
        rows.append({
            'loop_id': loop_id,
            'chrom1': chr1, 'start1': int(x1), 'end1': int(x2),
            'chrom2': chr2, 'start2': int(y1), 'end2': int(y2),
            'resolution': int(res),
            'category': row.get('category', 'unknown')
        })
    
    df_bedpe = pd.DataFrame(rows)
    df_bedpe['distance'] = (df_bedpe['start2'] + df_bedpe['end2'])/2 - (df_bedpe['start1'] + df_bedpe['end1'])/2
    
    # 3. Intersect Anchors with TSS +/- 200kb to assign Roles
    # Native Pandas logic to avoid bedtools system dependency
    df_tss = pd.read_csv(TSS_BED, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    # Expand TSS by 200kb
    df_tss['start_200'] = np.maximum(0, df_tss['start'] - 200000)
    df_tss['end_200'] = df_tss['end'] + 200000
    
    a1_genes = {}
    a2_genes = {}
    
    for chrom in set(df_bedpe['chrom1']).intersection(set(df_tss['chrom'])):
        # Anchor 1
        b_sub = df_bedpe[df_bedpe['chrom1'] == chrom]
        t_sub = df_tss[df_tss['chrom'] == chrom]
        
        b_starts = b_sub['start1'].values[:, None]
        b_ends = b_sub['end1'].values[:, None]
        t_starts = t_sub['start_200'].values[None, :]
        t_ends = t_sub['end_200'].values[None, :]
        
        mask = (b_starts <= t_ends) & (b_ends >= t_starts)
        b_idx_rel, t_idx_rel = np.where(mask)
        
        b_idx_abs = b_sub.index[b_idx_rel].values
        t_idx_abs = t_sub.index[t_idx_rel].values
        
        for b_i, t_i in zip(b_idx_abs, t_idx_abs):
            loop_id = df_bedpe.at[b_i, 'loop_id']
            gene_info = df_tss.at[t_i, 'name']
            a1_genes.setdefault(loop_id, set()).add(gene_info)
            
        # Anchor 2
        b_sub2 = df_bedpe[df_bedpe['chrom2'] == chrom]
        
        b_starts2 = b_sub2['start2'].values[:, None]
        b_ends2 = b_sub2['end2'].values[:, None]
        
        mask2 = (b_starts2 <= t_ends) & (b_ends2 >= t_starts)
        b_idx_rel2, t_idx_rel2 = np.where(mask2)
        
        b_idx_abs2 = b_sub2.index[b_idx_rel2].values
        t_idx_abs2 = t_sub.index[t_idx_rel2].values
        
        for b_i, t_i in zip(b_idx_abs2, t_idx_abs2):
            loop_id = df_bedpe.at[b_i, 'loop_id']
            gene_info = df_tss.at[t_i, 'name']
            a2_genes.setdefault(loop_id, set()).add(gene_info)
            
    # Assign Roles to Dataframe
    def get_role_and_genes(loop_id, anchor_genes_map):
        genes = anchor_genes_map.get(loop_id, set())
        if len(genes) > 0:
            return "Promoter_Anchor", ",".join(sorted(list(genes)))
        else:
            return "Distal_Anchor", ""

    df_bedpe['anchor1_role'] = ""
    df_bedpe['anchor1_genes'] = ""
    df_bedpe['anchor2_role'] = ""
    df_bedpe['anchor2_genes'] = ""
    
    for idx, row in df_bedpe.iterrows():
        lid = row['loop_id']
        r1, g1 = get_role_and_genes(lid, a1_genes)
        r2, g2 = get_role_and_genes(lid, a2_genes)
        
        df_bedpe.at[idx, 'anchor1_role'] = r1
        df_bedpe.at[idx, 'anchor1_genes'] = g1
        df_bedpe.at[idx, 'anchor2_role'] = r2
        df_bedpe.at[idx, 'anchor2_genes'] = g2
        
    # Add interaction_type
    def get_type(r1, r2):
        if r1 == "Promoter_Anchor" and r2 == "Promoter_Anchor":
            return "P-P"
        elif r1 == "Distal_Anchor" and r2 == "Distal_Anchor":
            return "E-E"
        else:
            return "P-E"
            
    df_bedpe['interaction_type'] = df_bedpe.apply(lambda x: get_type(x['anchor1_role'], x['anchor2_role']), axis=1)

    out_bedpe_path = os.path.join(CP1_DIR, "paper1_input", "loops_paper1_pe_15085_with_roles.tsv")
    df_bedpe.to_csv(out_bedpe_path, sep='\t', index=False)
    log(f"Saved parsed BEDPE with roles to {out_bedpe_path}")
    
    # 4. Create loop_qc_summary.tsv
    qc_path = os.path.join(CP1_DIR, "metadata", "loop_qc_summary.tsv")
    qc_data = {
        'step': [
            'Total PE Loops Parsed',
            'P-E (Promoter-Enhancer) Loops',
            'P-P (Promoter-Promoter) Loops',
            'E-E (Enhancer-Enhancer) Loops'
        ],
        'count': [
            len(df_bedpe),
            len(df_bedpe[df_bedpe['interaction_type'] == 'P-E']),
            len(df_bedpe[df_bedpe['interaction_type'] == 'P-P']),
            len(df_bedpe[df_bedpe['interaction_type'] == 'E-E'])
        ]
    }
    pd.DataFrame(qc_data).to_csv(qc_path, sep='\t', index=False)
    log(f"Saved QC summary to {qc_path}")
    
    log("Task 1 completed successfully.")

if __name__ == "__main__":
    main()
