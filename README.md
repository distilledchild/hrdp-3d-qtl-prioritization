# HRDP Frontal Cortex Hi-C: 3D QTL-to-Gene Prioritization

## Overview
This repository contains the workflows, pipelines, and intermediate outputs for the HRDP Frontal Cortex Hi-C project. The primary goal of this research is to leverage 3D chromatin contacts from 10 Hybrid Rat Diversity Panel (HRDP) strains to refine and prioritize the target genes of disease-model QTLs, moving beyond conventional 1D positional annotations.

By incorporating high-confidence promoter-enhancer (PE) loops, public epigenomic datasets (ATAC-seq, H3K27ac), and rigorous matched null-model permutations, this framework establishes biologically plausible disease-candidate genes for translational research.

## Pipeline Architecture
The analysis is structured into multiple robust phases, managed via `cooltools` and `Juicer` processing steps:

- **Phase H1:** Metadata curation and raw input auditing.
- **Phase H2:** Hi-C Quality Control (QC) and PE-Loop calling using HiCCUPS across multiple resolutions (5k, 10k, 25k).
- **Phase H3:** Multi-resolution `.mcool` matrix generation with ICE balancing. Built with custom local SSD caching to handle massive I/O operations smoothly over cloud-storage directories.
- **Phase H4:** A/B Compartment analysis (`cooltools eigs-cis`) at 1Mb resolution, phased by GC content to assess genome-wide architecture context.
- **Phase H5:** TAD and Insulation analysis (`cooltools insulation`) at 100kb and 250kb resolutions.
- **Phase H6 (Core Workflow):** The central QTL-to-Gene assignment framework. This includes Aggregate Peak Analysis (APA) validation of loops, anchoring QTL intervals to 3D contacts, and establishing a statistically rigorous matched null model to validate trait-domain enrichments.

## Repository Structure
- `HICCUPS/`: HiCCUPS loop calling outputs (`merged_loops.bedpe`) for each of the 10 HRDP strains.
- `QC_with_sb_options/`: Juicer quality control and pipeline metrics (e.g., `inter_30.txt`).
- `phase-h3/`: Contains `run_matrix_gen.py` for `.mcool` generation.
- `phase-h4/`: Contains `run_compartment_analysis.py` for A/B compartment calculations.
- `phase-h5/`: Contains `run_insulation_analysis.py` for TAD boundary calling.

## Key Considerations
- Replicates are not available; thus, all comparative and differential architecture claims are conservatively framed using "candidate" language.
- The pipeline utilizes `mRatBN7.2 (rn7)` as the core reference genome assembly.
- The primary claim of this research focuses on **QTL-to-gene prioritization**, not merely strain-specific architectural differences.
