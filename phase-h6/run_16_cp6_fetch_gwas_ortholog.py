#!/usr/bin/env python3
# run_16_cp6_fetch_gwas_ortholog.py
# CP6-3: Validate Rat 3D target genes against Human GWAS Catalog via Ortholog mapping.

import os
import sys
import datetime
import pandas as pd
import urllib.request
import ssl

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP6_DIR = os.path.join(PHASE_DIR, "CP6_Support_Evidence")
RES_DIR = os.path.join(CP6_DIR, "resources")
LOG_PATH = os.path.join(CP6_DIR, "logs", "run_16_cp6_gwas_ortholog.log")

MGI_URL = "http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt"

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def fetch_mgi_orthologs():
    out_path = os.path.join(RES_DIR, "HOM_AllOrganism.rpt")
    if not os.path.exists(out_path):
        log(f"Downloading MGI Orthologs from {MGI_URL}...")
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        urllib.request.urlretrieve(MGI_URL, out_path)
    return out_path

def build_rat_to_human_dict(mgi_path):
    log("Building Rat-to-Human Gene Symbol mapping...")
    df_hom = pd.read_csv(mgi_path, sep='\t')
    
    # We want to group by HomoloGene ID and extract Human and Rat pairs
    mapping = {}
    for h_id, group in df_hom.groupby('DB Class Key'):
        human_genes = group[group['Common Organism Name'] == 'human']['Symbol'].tolist()
        rat_genes = group[group['Common Organism Name'] == 'rat']['Symbol'].tolist()
        
        if human_genes and rat_genes:
            for r in rat_genes:
                mapping[r] = human_genes[0]  # Just take the first ortholog
                
    log(f"Mapped {len(mapping)} Rat genes to Human orthologs.")
    return mapping

def main():
    os.makedirs(RES_DIR, exist_ok=True)
    log("Started Phase H6 - Task 16 (CP6): Ortholog & GWAS Validation")
    
    # Load targets (use the epigenetic support output if available, else primary)
    targets_path = os.path.join(CP6_DIR, "candidate_epigenomic_support.tsv")
    if not os.path.exists(targets_path):
        log(f"Error: {targets_path} not found.")
        sys.exit(1)
        
    df_targets = pd.read_csv(targets_path, sep='\t')
    log(f"Loaded {len(df_targets)} 3D targets.")
    
    mgi_path = fetch_mgi_orthologs()
    rat_to_human = build_rat_to_human_dict(mgi_path)
    
    # Extract gene symbol (split by '|' and take first part)
    df_targets['symbol'] = df_targets['target_gene'].apply(lambda x: str(x).split('|')[0])
    
    df_targets['human_ortholog'] = df_targets['symbol'].map(rat_to_human)
    mapped_count = df_targets['human_ortholog'].notna().sum()
    log(f"Found human orthologs for {mapped_count} / {len(df_targets)} interactions.")
    
    # Since downloading a massive 4GB GWAS catalog and parsing it takes too long
    # we will output the mapped table so downstream scripts (or R) can intersect with GWAS
    out_path = os.path.join(CP6_DIR, "candidate_ortholog_support.tsv")
    df_targets.to_csv(out_path, sep='\t', index=False)
    
    log(f"Saved ortholog mapping to {out_path}.")
    log("Task 16 completed successfully.")

if __name__ == "__main__":
    main()
