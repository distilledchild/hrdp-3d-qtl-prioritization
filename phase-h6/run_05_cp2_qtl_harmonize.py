#!/usr/bin/env python3
# run_05_cp2_qtl_harmonize.py
# Parse RGD QTLs, extract mRatBN7.2 (rn7) coordinates, and split into primary/secondary BEDs.
import os
import sys
import datetime
import pandas as pd
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP2_DIR = os.path.join(PHASE_DIR, "CP2_QTL_Universe")
RAW_TXT = os.path.join(CP2_DIR, "metadata", "QTLS_RAT_raw.txt")
LOG_PATH = os.path.join(CP2_DIR, "logs", "run_05_cp2_qtl_harmonize.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def main():
    log("Started Phase H6 - Task 5 (CP2): QTL Harmonization & Filtering")
    
    if not os.path.exists(RAW_TXT):
        log(f"Error: {RAW_TXT} not found.")
        sys.exit(1)

    # RGD Columns mapping based on raw header comments
    # 0: QTL_RGD_ID, 2: QTL_SYMBOL, 14: TRAIT_NAME, 15: MEASUREMENT_TYPE
    # 46: 7.2_MAP_POS_CHR, 47: 7.2_MAP_POS_START, 48: 7.2_MAP_POS_STOP
    # The header line is the last line starting with # before data.
    # It's safer to load data using header=None and assign indices
    
    log("Parsing QTLS_RAT_raw.txt")
    
    # Read raw lines to bypass inconsistent header format
    with open(RAW_TXT, 'r', encoding='latin-1') as f:
        data_lines = [l.strip().split('\t') for l in f if not l.startswith('#') and l.strip()]
        
    df = pd.DataFrame(data_lines)
    
    # Check if we have at least 50 columns
    if len(df.columns) < 50:
        log("Error: Expected at least 50 columns in RGD file.")
        sys.exit(1)
        
    # Extract mapped columns
    try:
        df_clean = pd.DataFrame({
            'qtl_id': df[0],
            'symbol': df[2],
            'trait': df[14],
            'measurement': df[15],
            'chrom': df[46],
            'start': df[47],
            'end': df[48]
        })
    except KeyError as e:
        log(f"Error mapping columns: {e}")
        sys.exit(1)
        
    log(f"Initial QTL count: {len(df_clean)}")
    
    # Filter empty rn7 coordinates
    df_clean = df_clean.replace('', np.nan)
    df_mapped = df_clean.dropna(subset=['chrom', 'start', 'end']).copy()
    log(f"QTLs with mRatBN7.2 (rn7) coordinates: {len(df_mapped)}")
    
    # Clean coordinates
    df_mapped['chrom'] = df_mapped['chrom'].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x))
    df_mapped['start'] = pd.to_numeric(df_mapped['start'], errors='coerce').astype('Int64')
    df_mapped['end'] = pd.to_numeric(df_mapped['end'], errors='coerce').astype('Int64')
    
    df_mapped = df_mapped.dropna(subset=['start', 'end']).copy()
    
    # Calculate size and filter
    df_mapped['size_mb'] = (df_mapped['end'] - df_mapped['start']) / 1e6
    
    # Remove single-point or negative QTLs
    df_mapped = df_mapped[df_mapped['size_mb'] > 0].copy()
    
    # Split into primary (<10Mb) and secondary (>=10Mb)
    primary = df_mapped[df_mapped['size_mb'] < 10].copy()
    secondary = df_mapped[df_mapped['size_mb'] >= 10].copy()
    
    log(f"Primary QTLs (<10Mb): {len(primary)}")
    log(f"Broad/Secondary QTLs (>=10Mb): {len(secondary)}")
    
    # Format as BED: chrom, start, end, name, score, strand
    primary_bed = pd.DataFrame({
        'chrom': primary['chrom'],
        'start': primary['start'],
        'end': primary['end'],
        'name': primary['symbol'] + "|" + primary['trait'],
        'score': primary['size_mb'].round(2),
        'strand': '.'
    })
    
    secondary_bed = pd.DataFrame({
        'chrom': secondary['chrom'],
        'start': secondary['start'],
        'end': secondary['end'],
        'name': secondary['symbol'] + "|" + secondary['trait'],
        'score': secondary['size_mb'].round(2),
        'strand': '.'
    })
    
    # Save files
    primary_path = os.path.join(CP2_DIR, "qtl_filtered_primary.bed")
    secondary_path = os.path.join(CP2_DIR, "qtl_broad_secondary.bed")
    
    primary_bed.to_csv(primary_path, sep='\t', index=False, header=False)
    secondary_bed.to_csv(secondary_path, sep='\t', index=False, header=False)
    
    log(f"Saved {primary_path}")
    log(f"Saved {secondary_path}")
    
    # QC Summary
    qc_path = os.path.join(CP2_DIR, "metadata", "qtl_liftover_qc.tsv")
    qc_data = pd.DataFrame({
        'metric': ['Total RGD QTLs', 'With rn7 coordinates', 'Valid intervals', 'Primary (<10Mb)', 'Secondary (>=10Mb)'],
        'count': [len(df_clean), len(df_clean.dropna(subset=['chrom'])), len(df_mapped), len(primary), len(secondary)]
    })
    qc_data.to_csv(qc_path, sep='\t', index=False)
    log(f"Saved QC summary to {qc_path}")
    
    log("Task 5 completed successfully.")

if __name__ == "__main__":
    main()
