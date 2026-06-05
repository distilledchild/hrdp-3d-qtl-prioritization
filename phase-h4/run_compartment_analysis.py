import os
import subprocess
import sys
import datetime

# =======================================================================
# Configuration - Local Mac Environment
# =======================================================================
# Dropbox output directory
DROPBOX_H4_DIR = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic2/phase-h4"
os.makedirs(DROPBOX_H4_DIR, exist_ok=True)

# Set H3 matrix path (Dropbox location where phase-h3 outputs are saved)
H3_MCOOL_DIR = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic2/phase-h3"

# GC track file downloaded from HPC and placed in Dropbox H4 dir
GC_FILE = os.path.join(DROPBOX_H4_DIR, "rn7_gc_1Mb.tsv")

RESOLUTION  = 1000000 # 1Mb
NUM_THREADS = 4
SAMPLES     = ["607", "DA21A", "DBA9A", "DA68A", "D765A", "592BB", "DA08A", "A2DB", "74AA", "DE8BA"]

# Log file setup (generated including time)
LOG_FILE = os.path.join(DROPBOX_H4_DIR, f"phase_h4_compartment_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
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
    # Move to Dropbox output directory
    os.chdir(DROPBOX_H4_DIR)
    
    # Initialize log file
    with open(LOG_FILE, "w") as f:
        f.write(f"=== Phase H4 Compartment Analysis Log ({datetime.datetime.now()}) ===\n")
        
    print_log(f"📂 Output Directory (Dropbox): {DROPBOX_H4_DIR}")
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
            run_cmd(f"uv run --python 3.10 --with cooltools==0.7.1 --with \"pandas<2.2.0\" --with \"numpy<2.0.0\" cooltools expected-cis \"{mcool_path}::/resolutions/{RESOLUTION}\" -p {NUM_THREADS} -o \"{expected_file}\"")
        else:
            print_log(f"[{sample}] Step A: Already completed.")

        # [Step B] Calculate Eigenvectors (Oriented with GC)
        vecs_file = f"{compartment_prefix}.cis.vecs.tsv"
        if not os.path.exists(vecs_file):
            print_log(f"[{sample}] Step B: Calculating Eigenvectors (A/B Compartment)...")
            run_cmd(f"uv run --python 3.10 --with cooltools==0.7.1 --with \"pandas<2.2.0\" --with \"numpy<2.0.0\" cooltools eigs-cis \"{mcool_path}::/resolutions/{RESOLUTION}\" --phasing-track \"{GC_FILE}::GC\" -o \"{compartment_prefix}\"")
        else:
            print_log(f"[{sample}] Step B: Already completed.")

        # [Step C] Calculate Saddle data (no --fig; PNG is generated separately in Step D)
        saddle_out = f"{saddle_prefix}.saddledump.npz"
        if not os.path.exists(saddle_out):
            print_log(f"[{sample}] Step C: Computing Saddle data...")
            run_cmd(f"uv run --python 3.10 --with cooltools==0.7.1 --with \"pandas<2.2.0\" --with \"numpy<2.0.0\" cooltools saddle \"{mcool_path}::/resolutions/{RESOLUTION}\" \"{vecs_file}::E1\" \"{expected_file}\" --contact-type cis --qrange 0.02 0.98 --n-bins 50 -o \"{saddle_prefix}\"")
        else:
            print_log(f"[{sample}] Step C: Already completed.")

        # [Step D] Generate publication-quality Saddle Plot PNG
        # Follows the official cooltools tutorial: https://cooltools.readthedocs.io/en/latest/notebooks/compartments_and_saddles.html
        saddle_png = f"{saddle_prefix}.png"
        if not os.path.exists(saddle_png):
            print_log(f"[{sample}] Step D: Rendering Saddle Plot PNG...")
            vecs_file_path = f"{compartment_prefix}.cis.vecs.tsv"
            digitized_file = f"{saddle_prefix}.digitized.tsv"
            plot_script = f'''
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib import ticker

# --- Load data ---
data = np.load("{saddle_out}")
S = data["saddledata"]

# Load eigenvector track (raw E1 values)
vecs = pd.read_csv("{vecs_file_path}", sep="\\t")
# Load digitized track (bin assignments)
digi = pd.read_csv("{digitized_file}", sep="\\t")

# --- Compute groupmean (mean E1 per digitized bin) ---
e1_col = "E1"  # explicitly reference E1, avoid column index mistakes
digi_col = digi.columns[3]
groupmean = vecs[e1_col].groupby(digi[digi_col]).mean()

# --- Quantile bin edges ---
Q_LO, Q_HI, N_BINS = 0.02, 0.98, 50
binedges = np.linspace(Q_LO, Q_HI, N_BINS + 1)

# Trim outlier flanking rows/cols from saddle matrix
n = S.shape[0]
C = S[1:-1, 1:-1] if (n - N_BINS) == 2 else S

# Extract only the N_BINS valid groups (indices 1..N_BINS)
gm = groupmean.reindex(range(1, N_BINS + 1)).values
bw = np.diff(binedges)
be = binedges[:-1]

X, Y = np.meshgrid(binedges, binedges)
vmin, vmax = 0.5, 2.0

# --- Layout ---
fig = plt.figure(figsize=(5, 5))
gs = GridSpec(nrows=3, ncols=3,
              width_ratios=[0.2, 1, 0.08],
              height_ratios=[0.2, 1, 0.08],
              wspace=0.05, hspace=0.05)

# --- Heatmap ---
ax = plt.subplot(gs[4])
norm = LogNorm(vmin=vmin, vmax=vmax)
img = ax.pcolormesh(X, Y, C, cmap="coolwarm", norm=norm, rasterized=True)
ax.yaxis.set_visible(False)
ax.set_xlim(Q_LO, Q_HI)
ax.set_ylim(Q_HI, Q_LO)
ax.grid(False)
ax.set_xlabel("E1 quantiles")

ticks = [0.2, 0.4, 0.6, 0.8]
ax.set_xticks(ticks)
ax.set_xticklabels(ticks)

# --- Left margin (sharey with heatmap) ---
ax_left = plt.subplot(gs[3], sharey=ax)
ax_left.barh(be, height=bw, width=gm, align="edge",
             edgecolor="k", facecolor=None, linewidth=1)
ax_left.set_xlim(ax_left.get_xlim()[1], ax_left.get_xlim()[0])  # flip
ax_left.set_ylim(Q_HI, Q_LO)
ax_left.set_yticks(ticks)
for sp in ["top","bottom","left"]: ax_left.spines[sp].set_visible(False)
ax_left.xaxis.set_visible(False)
ax_left.set_ylabel("E1 quantiles")
ax_left.axvline(0, color="gray", linestyle="-", linewidth=0.5)

# --- Top margin (sharex with heatmap) ---
ax_top = plt.subplot(gs[1], sharex=ax)
ax_top.bar(be, width=bw, height=gm, align="edge",
           edgecolor="k", facecolor=None, linewidth=1)
ax_top.set_xlim(Q_LO, Q_HI)
for sp in ["top","right","left"]: ax_top.spines[sp].set_visible(False)
ax_top.xaxis.set_visible(False)
ax_top.yaxis.set_visible(False)
ax_top.axhline(0, color="gray", linestyle="-", linewidth=0.5)

# --- Colorbar ---
class MinOneMaxFormatter(ticker.LogFormatter):
    def set_locs(self, locs=None):
        self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])
    def __call__(self, x, pos=None):
        if x not in [vmin, 1, vmax]: return ""
        else: return f"{{x:g}}"

ax_cbar = plt.subplot(gs[5])
cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=ax_cbar,
                  fraction=0.8, label="average observed/expected\\ncontact frequency")
cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

fig.savefig("{saddle_png}", dpi=300, bbox_inches="tight")
plt.close()
'''
            run_cmd(f"uv run --python 3.10 --with cooltools==0.7.1 --with \"pandas<2.2.0\" --with \"numpy<2.0.0\" --with matplotlib python3 -c '{plot_script}'")
        else:
            print_log(f"[{sample}] Step D: Already completed.")

    print_log("\n✅ All H4 Compartment analysis scripts have finished execution for all samples!")

if __name__ == "__main__":
    main()
