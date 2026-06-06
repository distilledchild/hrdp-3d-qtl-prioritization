#!/usr/bin/env python3
# run_07_cp3_fetch_genes.py
# CP3: Fetch Ensembl rn7 GTF, filter protein-coding genes, and create Gene/TSS BED files.
import os
import sys
import gzip
import urllib.request
import datetime
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP3_DIR = os.path.join(PHASE_DIR, "CP3_Positional_Baseline")
REF_DIR = os.path.join(CP3_DIR, "reference")
LOG_PATH = os.path.join(CP3_DIR, "logs", "run_07_cp3_fetch_genes.log")

GTF_URL = "http://ftp.ensembl.org/pub/release-112/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gtf.gz"
GTF_GZ = os.path.join(REF_DIR, "Rattus_norvegicus.mRatBN7.2.112.gtf.gz")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def parse_gtf_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(';'):
        item = item.strip()
        if not item: continue
        parts = item.split(' ', 1)
        if len(parts) == 2:
            key = parts[0].strip()
            val = parts[1].strip().strip('"')
            attrs[key] = val
    return attrs

def main():
    os.makedirs(REF_DIR, exist_ok=True)
    log("Started Phase H6 - Task 7 (CP3): Fetch and Parse Reference Genes")

    # 1. Download
    if not os.path.exists(GTF_GZ):
        log(f"Downloading GTF from {GTF_URL}...")
        try:
            urllib.request.urlretrieve(GTF_URL, GTF_GZ)
            log("Download completed.")
        except Exception as e:
            log(f"Download failed: {e}")
            sys.exit(1)
    else:
        log("GTF already downloaded.")

    # 2. Parse GTF for protein coding genes
    log("Parsing GTF for protein-coding genes...")
    genes = []
    with gzip.open(GTF_GZ, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) != 9: continue
            
            chrom, source, feature, start, end, score, strand, phase, attr_str = fields
            if feature != 'gene': continue
            
            # Filter chromosome names (keep 1-20, X, Y, MT)
            if not chrom.isalnum() and not chrom in ['MT']:
                continue
                
            attrs = parse_gtf_attributes(attr_str)
            gene_biotype = attrs.get('gene_biotype', '')
            
            # Ensembl GTF uses 'protein_coding' for biotype
            if gene_biotype == 'protein_coding':
                gene_id = attrs.get('gene_id', 'unknown')
                gene_name = attrs.get('gene_name', gene_id)
                # Ensure UCSC-style chr prefix
                chrom_name = f"chr{chrom}" if chrom != 'MT' else 'chrM'
                
                genes.append({
                    'chrom': chrom_name,
                    'start': int(start) - 1, # 0-indexed BED
                    'end': int(end),
                    'name': f"{gene_name}|{gene_id}",
                    'score': 0,
                    'strand': strand
                })

    df_genes = pd.DataFrame(genes)
    log(f"Extracted {len(df_genes)} protein-coding genes.")

    # 3. Create Gene BED
    gene_bed_path = os.path.join(REF_DIR, "rn7_protein_coding_genes.bed")
    df_genes.to_csv(gene_bed_path, sep='\t', index=False, header=False)
    log(f"Saved Gene BED: {gene_bed_path}")

    # 4. Create TSS BED (1bp size)
    # If strand is '+', TSS is start to start+1. If '-', TSS is end-1 to end.
    df_tss = df_genes.copy()
    
    is_pos = df_tss['strand'] == '+'
    df_tss.loc[is_pos, 'end'] = df_tss.loc[is_pos, 'start'] + 1
    
    is_neg = df_tss['strand'] == '-'
    df_tss.loc[is_neg, 'start'] = df_tss.loc[is_neg, 'end'] - 1
    
    tss_bed_path = os.path.join(REF_DIR, "rn7_protein_coding_tss.bed")
    df_tss.to_csv(tss_bed_path, sep='\t', index=False, header=False)
    log(f"Saved TSS BED: {tss_bed_path}")

    log("Task 7 completed successfully.")

if __name__ == "__main__":
    main()
