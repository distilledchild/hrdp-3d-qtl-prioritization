import os
import subprocess
import sys
import datetime

# =======================================================================
# Configuration - Local Mac Environment
# =======================================================================
# Automatically recognize the folder where the script is located as the working directory
WORK_DIR = os.path.dirname(os.path.abspath(__file__))

# Set H3 matrix path (phase-h3 located in the parent directory of phase-h4 folder)
H3_MCOOL_DIR = os.path.join(os.path.dirname(WORK_DIR), "phase-h3")

# GC track file downloaded from HPC and placed in WORK_DIR
GC_FILE = os.path.join(WORK_DIR, "rn7_gc_1Mb.tsv")

RESOLUTION  = 1000000 # 1Mb
NUM_THREADS = 4
SAMPLES     = ["607", "DA21A", "DBA9A", "DA68A", "DE8BA", "D765A", "592BB", "DA08A", "A2DB", "74AA"]

# Log file setup (generated including time)
LOG_FILE = os.path.join(WORK_DIR, f"phase_h4_compartment_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
# =======================================================================

def print_log(msg):
    """Print to terminal and record in log file simultaneously"""
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(msg + "\n")

def run_cmd(cmd):
    """Execute command and capture output to terminal and log file simultaneously"""
    print_log(f"Running: {cmd}")
    try:
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            sys.stdout.write(line)
            with open(LOG_FILE, "a") as f:
                f.write(line)
        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)
    except Exception as e:
        print_log(f"❌ Command execution failed: {cmd}\nError message: {str(e)}")
        sys.exit(1)

def main():
    # Move to working directory (script location)
    os.chdir(WORK_DIR)
    
    # Initialize log file
    with open(LOG_FILE, "w") as f:
        f.write(f"=== Phase H4 Compartment Analysis Log ({datetime.datetime.now()}) ===\n")
        
    print_log(f"📂 Working Directory: {WORK_DIR}")
    print_log(f"📝 Log File: {LOG_FILE}")

    # Check if GC file exists first
    if not os.path.exists(GC_FILE):
        print_log(f"🚨 Error: GC track file downloaded from HPC is missing!")
        print_log(f"Please place the file exactly at the following path: {GC_FILE}")
        sys.exit(1)

    # H4 Compartment calculation loop for each sample
    for sample in SAMPLES:
        mcool_path = os.path.join(H3_MCOOL_DIR, f"{sample}.mcool")
        if not os.path.exists(mcool_path):
            print_log(f"\n⚠️ Warning: Cannot find mcool file for {sample} ({mcool_path}). Skipping.")
            continue
            
        print_log(f"\n================ [ {sample} ] Starting Compartment Analysis ================")
        
        expected_file      = f"{sample}_expected_1Mb.tsv"
        compartment_prefix = f"{sample}_compartment"
        saddle_prefix      = f"{sample}_saddle"

        # [Step A] Calculate Expected-cis
        if not os.path.exists(expected_file):
            print_log(f"[{sample}] Step A: Calculating expected-cis...")
            run_cmd(f"uv run cooltools expected-cis \"{mcool_path}::/resolutions/{RESOLUTION}\" -p {NUM_THREADS} -o \"{expected_file}\"")
        else:
            print_log(f"[{sample}] Step A: Already completed.")

        # [Step B] Calculate Eigenvectors (Oriented with GC)
        vecs_file = f"{compartment_prefix}.cis.vecs.tsv"
        if not os.path.exists(vecs_file):
            print_log(f"[{sample}] Step B: Calculating Eigenvectors (A/B Compartment)...")
            run_cmd(f"uv run cooltools eigs-cis \"{mcool_path}::/resolutions/{RESOLUTION}\" --phasing-track \"{GC_FILE}\" --phasing-col GC -o \"{compartment_prefix}\"")
        else:
            print_log(f"[{sample}] Step B: Already completed.")

        # [Step C] Calculate Saddle Plot and Strength matrix
        saddle_out = f"{saddle_prefix}.saddledump.npz"
        if not os.path.exists(saddle_out):
            print_log(f"[{sample}] Step C: Generating Saddle Plot...")
            run_cmd(f"uv run cooltools saddle \"{mcool_path}::/resolutions/{RESOLUTION}\" \"{expected_file}\" \"{vecs_file}\" --contact-type cis --qrange 0.02 0.98 --n-bins 50 --fig html -o \"{saddle_prefix}\"")
        else:
            print_log(f"[{sample}] Step C: Already completed.")

    print_log("\n✅ All H4 Compartment analysis scripts have finished execution for all samples!")

if __name__ == "__main__":
    main()
