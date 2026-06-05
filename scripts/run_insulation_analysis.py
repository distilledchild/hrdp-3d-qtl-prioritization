import os
import subprocess
import sys
import datetime

# =======================================================================
# Configuration - Local Mac Environment
# =======================================================================
# Automatically recognize the folder where the script is located as the working directory
WORK_DIR = os.path.dirname(os.path.abspath(__file__))

# Set H3 matrix path (phase-h3 located in the parent directory of phase-h5 folder)
H3_MCOOL_DIR = os.path.join(os.path.dirname(WORK_DIR), "phase-h3")

SAMPLES = ["607", "DA21A", "DBA9A", "DA68A", "DE8BA", "D765A", "592BB", "DA08A", "A2DB", "74AA"]

# Window sizes for 100kb resolution (typically 3x, 5x, 10x of resolution)
RES_100K     = 100000
WINDOWS_100K = ["300000", "500000", "1000000"]

# Window sizes for 250kb resolution
RES_250K     = 250000
WINDOWS_250K = ["750000", "1250000", "2500000"]

NUM_THREADS = 4

# Log file setup (generated including time)
LOG_FILE = os.path.join(WORK_DIR, f"phase_h5_insulation_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
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
    os.makedirs(WORK_DIR, exist_ok=True)
    os.chdir(WORK_DIR)
    
    # Initialize log file
    with open(LOG_FILE, "w") as f:
        f.write(f"=== Phase H5 Insulation Analysis Log ({datetime.datetime.now()}) ===\n")
        
    print_log(f"📂 Working Directory: {WORK_DIR}")
    print_log(f"📝 Log File: {LOG_FILE}")

    # TAD/Insulation calculation loop for each sample
    for sample in SAMPLES:
        mcool_path = os.path.join(H3_MCOOL_DIR, f"{sample}.mcool")
        if not os.path.exists(mcool_path):
            print_log(f"\n⚠️ Warning: Cannot find mcool file for {sample} ({mcool_path}). Skipping.")
            continue
            
        print_log(f"\n================ [ {sample} ] Insulation & TAD Boundary Analysis ================")
        
        out_100k = f"{sample}_insulation_100kb.tsv"
        out_250k = f"{sample}_insulation_250kb.tsv"

        # [Step A] 100kb resolution analysis
        if not os.path.exists(out_100k):
            print_log(f"[{sample}] Step A: 100kb resolution analysis in progress... (Windows: {', '.join(WINDOWS_100K)})")
            cmd_100k = f"uv run cooltools insulation \"{mcool_path}::/resolutions/{RES_100K}\" {' '.join(WINDOWS_100K)} -p {NUM_THREADS} -o \"{out_100k}\""
            run_cmd(cmd_100k)
        else:
            print_log(f"[{sample}] Step A: 100kb resolution analysis is already completed.")

        # [Step B] 250kb resolution analysis
        if not os.path.exists(out_250k):
            print_log(f"[{sample}] Step B: 250kb resolution analysis in progress... (Windows: {', '.join(WINDOWS_250K)})")
            cmd_250k = f"uv run cooltools insulation \"{mcool_path}::/resolutions/{RES_250K}\" {' '.join(WINDOWS_250K)} -p {NUM_THREADS} -o \"{out_250k}\""
            run_cmd(cmd_250k)
        else:
            print_log(f"[{sample}] Step B: 250kb resolution analysis is already completed.")

    print_log("\n✅ Individual calculations for Phase H5 (Insulation/TAD) have been completed for all samples!")
    print_log("You can proceed with integrated analysis to find Constitutive/Variable Boundaries by collecting the calculated TSV files.")

if __name__ == "__main__":
    main()
