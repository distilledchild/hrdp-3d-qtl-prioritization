import pandas as pd
import gzip
import subprocess
import os

input_gz = "resources/GSE134935_ATAC_4samples_ucsctracks.bed.gz"
temp_rn6 = "resources/temp_rn6.bed"
out_rn7 = "resources/rn7_GSE134935_ATAC.bed"
unmapped = "resources/unmapped.bed"

print("Parsing GSE134935 ATAC summits...")
records = []
with gzip.open(input_gz, 'rt') as f:
    for line in f:
        if line.startswith("track"):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            # ATAC summits are 1bp or 0bp. Expand by +/- 250bp to represent a 500bp peak.
            start = max(0, start - 250)
            end = end + 250
            records.append([chrom, start, end])

df = pd.DataFrame(records, columns=['chrom', 'start', 'end'])
# Remove duplicates (since there are 4 samples, many summits overlap)
df = df.drop_duplicates()
df.to_csv(temp_rn6, sep='\t', index=False, header=False)
print(f"Wrote {len(df)} peaks to {temp_rn6}")

print("Running liftOver to rn7...")
cmd = [
    "./resources/liftOver",
    temp_rn6,
    "resources/rn6ToRn7.over.chain.gz",
    out_rn7,
    unmapped
]
subprocess.run(cmd, check=True)

print("Sorting and merging lifted over peaks...")
# Reload rn7 and sort
df_rn7 = pd.read_csv(out_rn7, sep='\t', header=None, names=['chrom', 'start', 'end'])
df_rn7 = df_rn7.sort_values(by=['chrom', 'start'])

# Basic merging of overlapping peaks using pure pandas
# This prevents double counting
def merge_intervals(group):
    # Sort by start (already sorted, but safe)
    starts = group['start'].values
    ends = group['end'].values
    merged = []
    
    current_start = starts[0]
    current_end = ends[0]
    
    for st, en in zip(starts[1:], ends[1:]):
        if st <= current_end:
            current_end = max(current_end, en)
        else:
            merged.append([group.name, current_start, current_end])
            current_start = st
            current_end = en
    merged.append([group.name, current_start, current_end])
    return pd.DataFrame(merged, columns=['chrom', 'start', 'end'])

df_merged = df_rn7.groupby('chrom').apply(merge_intervals).reset_index(drop=True)
df_merged.to_csv(out_rn7, sep='\t', index=False, header=False)

print(f"Final merged ATAC-seq peaks (rn7): {len(df_merged)}")
print(f"Saved to {out_rn7}")

# Cleanup temp files
os.remove(temp_rn6)
os.remove(unmapped)
