#!/usr/bin/env python3
# run_06_cp2_trait_mapping.py
# Maps RGD QTL traits into broad phenotypic domains.
import os
import sys
import datetime
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP2_DIR = os.path.join(PHASE_DIR, "CP2_QTL_Universe")
LOG_PATH = os.path.join(CP2_DIR, "logs", "run_06_cp2_trait_mapping.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

# Keyword mapping rules (order matters: first match wins)
DOMAIN_RULES = [
    ("addiction/reward", ["addict", "alcohol", "cocaine", "reward", "ethanol", "morphine", "nicotine"]),
    ("anxiety/depression-like/behavior", ["anxiet", "depress", "behav", "fear", "locomotor", "open field", "startle", "avoidance", "social", "maze", "explore"]),
    ("cardiovascular", ["blood pressure", "heart", "cardio", "hypertension", "ventric", "vascular", "stroke", "aorta"]),
    ("metabolic", ["body weight", "obesity", "glucose", "insulin", "cholesterol", "lipid", "triglyceride", "metabol", "fat", "adipose"]),
    ("renal", ["kidney", "renal", "proteinuria", "albuminuria", "nephropathy"]),
    ("immune/inflammatory", ["autoimmun", "inflam", "arthritis", "asthma", "cytokine", "leukocyte", "lymphocyte", "macrophage", "ige", "igg"]),
    ("neurodevelopmental/neurological", ["seizure", "epileps", "brain", "neuro", "cortex", "hippocampus", "stroke", "ischemia"])
]

def map_domain(trait):
    trait_lower = str(trait).lower()
    for domain, keywords in DOMAIN_RULES:
        for kw in keywords:
            if kw in trait_lower:
                return domain
    return "other/unclear"

def process_bed(bed_path, out_path):
    if not os.path.exists(bed_path):
        log(f"Warning: {bed_path} not found. Skipping.")
        return None

    # BED format: chrom, start, end, name (symbol|trait), score, strand
    df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    # Extract symbol and trait
    df[['symbol', 'trait']] = df['name'].str.split('|', n=1, expand=True)
    df['trait'] = df['trait'].fillna("unknown")
    
    # Map domain
    df['domain'] = df['trait'].apply(map_domain)
    
    # Save mapping result
    df.to_csv(out_path, sep='\t', index=False)
    
    # Domain counts
    counts = df['domain'].value_counts().to_dict()
    return counts

def main():
    log("Started Phase H6 - Task 6 (CP2): Trait Domain Mapping")
    
    primary_bed = os.path.join(CP2_DIR, "qtl_filtered_primary.bed")
    primary_out = os.path.join(CP2_DIR, "qtl_primary_trait_mapping.tsv")
    
    secondary_bed = os.path.join(CP2_DIR, "qtl_broad_secondary.bed")
    secondary_out = os.path.join(CP2_DIR, "qtl_secondary_trait_mapping.tsv")
    
    counts_pri = process_bed(primary_bed, primary_out)
    counts_sec = process_bed(secondary_bed, secondary_out)
    
    if counts_pri:
        log("Primary QTL Domains:")
        for k, v in counts_pri.items():
            log(f"  {k}: {v}")
            
    if counts_sec:
        log("Secondary/Broad QTL Domains:")
        for k, v in counts_sec.items():
            log(f"  {k}: {v}")
            
    log("Task 6 completed successfully.")

if __name__ == "__main__":
    main()
