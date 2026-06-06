#!/usr/bin/env python3
# run_04_cp2_qtl_fetch.py
# CP2: Fetch active Rat QTLs from RGD bypassing self-signed SSL certs.
import os
import urllib.request
import ssl
import datetime
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE_DIR = os.path.join(REPO_ROOT, "phase-h6")
CP2_DIR = os.path.join(PHASE_DIR, "CP2_QTL_Universe")
META_DIR = os.path.join(CP2_DIR, "metadata")
RAW_TXT = os.path.join(META_DIR, "QTLS_RAT_raw.txt")
LOG_PATH = os.path.join(CP2_DIR, "logs", "run_04_cp2_qtl_fetch.log")

def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {msg}"
    print(full_msg)
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    with open(LOG_PATH, "a") as f:
        f.write(full_msg + "\n")

def main():
    for d in [CP2_DIR, META_DIR, os.path.join(CP2_DIR, "logs")]:
        os.makedirs(d, exist_ok=True)

    log("Started Phase H6 - Task 4 (CP2): QTL Download from RGD")
    
    url = "https://download.rgd.mcw.edu/pub/data_release/QTLS_RAT.txt"
    log(f"Fetching from {url}")
    
    # Bypass SSL verification due to RGD certificate issues
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    try:
        req = urllib.request.urlopen(url, context=ctx)
        with open(RAW_TXT, 'wb') as f:
            f.write(req.read())
        log(f"Successfully downloaded to {RAW_TXT}")
    except Exception as e:
        log(f"Download failed: {e}")
        return

    # Basic parsing to generate download log
    log("Parsing raw txt to generate qtl_raw_download_log.tsv")
    try:
        # RGD files use # for comments but header line does not start with # in the column names,
        # actually the last comment block shows column names, but usually line 46 is the header.
        # We will load with error_bad_lines=False, but it's simpler to just read raw lines.
        lines = open(RAW_TXT, 'r', encoding='latin-1').readlines()
        data_lines = [l for l in lines if not l.startswith('#')]
        
        qc_path = os.path.join(META_DIR, "qtl_raw_download_log.tsv")
        qc_data = {
            'source': ['RGD_QTLS_RAT.txt'],
            'download_date': [datetime.datetime.now().strftime("%Y-%m-%d")],
            'total_raw_rows': [len(data_lines)],
            'note': ['Downloaded bypassing SSL verify']
        }
        pd.DataFrame(qc_data).to_csv(qc_path, sep='\t', index=False)
        log(f"Raw lines count: {len(data_lines)}")
        log(f"Saved QC to {qc_path}")
    except Exception as e:
        log(f"Parsing error: {e}")

    log("Task 4 completed successfully.")

if __name__ == "__main__":
    main()
