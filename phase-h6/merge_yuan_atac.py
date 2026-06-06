import pandas as pd
import glob
import os

peaks_dir = "/Users/pete/Desktop/playground/epigenetic-analysis-2026/ATAC-seq_2021_Yuan_et_al/out_atac_rn7_pe"
peak_files = glob.glob(f"{peaks_dir}/*/peaks/*.narrowPeak")

dfs = []
for pf in peak_files:
    # MACS2 narrowPeak format: chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak
    df = pd.read_csv(pf, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    dfs.append(df)

df_all = pd.concat(dfs, ignore_index=True)
df_all = df_all.sort_values(by=['chrom', 'start'])

print(f"Total raw Yuan peaks: {len(df_all)}")

# Merge overlapping
def merge_intervals(group):
    starts = group['start'].values
    ends = group['end'].values
    merged = []
    
    if len(starts) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
        
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

df_merged = df_all.groupby('chrom').apply(merge_intervals).reset_index(drop=True)
print(f"Merged Yuan peaks: {len(df_merged)}")

out_path = "/Users/pete/Desktop/playground/hrdp-3d-qtl-prioritization/phase-h6/CP6_Support_Evidence/resources/rn7_Yuan2021_ATAC.bed"
df_merged.to_csv(out_path, sep='\t', index=False, header=False)
print(f"Saved to {out_path}")
