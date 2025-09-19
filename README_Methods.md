# Methods
## Bhavna Mahi
## Miga Lab

This README lists every task, command, descriptions, and file type for everything pertaining to my research. This section below lists all methods for bSat analysis.

### 1. Identify and Extract Repeat Units
```
bedtools getfasta -fi genome.fa -bed input.bed -fo output.fa -s
```
Extract sequences from a BED file into FASTA (strand-aware).

Input: BED + FASTA genome

Output: FASTA

```
cat file1.fa file2.fa > combined_seqs.fa
```
Concatenate multiple FASTA files into one.

Input: FASTA

Output: FASTA

### 2. Build and Search HMMs (HMMER)
```
hmmbuild myModel.hmm mySeq.alignment.fa
```
Build an HMM from an aligned FASTA file.

Input: aligned FASTA

Output: HMM

```
hmmbuild myModel1.hmm mySeq1.alignment.fa
hmmbuild myModel2.hmm mySeq2.alignment.fa
cat *hmm > modelName
hmmpress modelName
```
Build multiple HMMs, combine them, and prepare the database.

Input: aligned FASTA

Output: HMM, pressed database

```
time nhmmer --cpu 12 --notextw --noali --tblout genome.model.out modelName /path/to/genome.fa
```
Search genome FASTA against HMM model database.

Input: HMM + genome FASTA

Output: nhmmer table (.out)

```
awk -v th=0.5 -f hmmertblout2bed.awk genome.out > genome.model.bed
```
Convert nhmmer output into a BED file for viewing.

Input: nhmmer output

Output: BED

### 3. Filter BED and FASTA Sequences
```
awk '{print $3 - $2}' input.bed | sort -n
```
Print sizes of regions in BED file, sorted.

Input: BED

Output: list of lengths (stdout)

```
awk '($3 - $2) > 445' input.bed > filtered.bed
```
Filter out short regions (<445 bp).

Input: BED

Output: BED

```
awk '{ if(($3 - $2) >= 65 && ($3 - $2) <= 70) print }' input.bed > filtered.bed
```
Keep only sequences between 65–70 bp.

Input: BED

Output: BED

**Deduplication**
```
conda install bioconda::seqkit
seqkit rmdup -s < input.fa > dedup.fa
```
Remove duplicate sequences from a FASTA file.

Input: FASTA

Output: FASTA

### 4. Align Sequences
**Install MUSCLE**
```
conda install -c bioconda muscle=5.1
```
Install MUSCLE via bioconda.

**MUSCLE**
```
muscle -super combined_seqs.fa -output combined_seqs_aligned.fa
```
Align two sequences in a combined FASTA.

Input: FASTA

Output: aligned FASTA

```
muscle -super5 input.fa -output aligned.fa
```
Align many sequences with MUSCLE super5.

Input: FASTA

Output: aligned FASTA

```
muscle -super5 chr1bsatResults5_filtered_rmdup.fa -output chr1bsatResults5_filtered_rmdup_MUSCLE_aligned.fa
```
Align full deduplicated BSR set with MUSCLE.

Input: FASTA

Output: aligned FASTA

**MAFFT**
```
"/path/to/mafft" --auto --inputorder input.fa > aligned.fa
```
Align sequences with MAFFT (auto mode).

Input: FASTA

Output: aligned FASTA

**Clustal Omega**

(Submitted via EBI webserver — results downloaded)

Input: FASTA

Output: aligned FASTA

### 5. Subsampling (for testing trees)
```
conda install bioconda::seqtk
seqtk sample -s100 chr1bsatResults5_filtered_rmdup.fa 200 > chr1bsatResults5_filtered_rmdup_200.fa
```
Randomly subsample 200 sequences for faster test runs.

Input: FASTA

Output: FASTA

### 6. Build Phylogenetic Trees & Clusters
**With custom scripts**
```
python align_and_tree.py input.fa
```
Generate alignment, distance matrix, and trees from FASTA.

Input: FASTA

Output: aligned FASTA, distance TSV, Newick

```
python clustering_script.py input_dist.tsv
```
Cluster sequences based on distance matrix (default linkage).

Input: distance TSV

Output: cluster CSV, dendrogram PNG

```
python clustering_script_MP.py input.afa
```
Cluster sequences using Maximum Parsimony.

Input: aligned FASTA

Output: cluster CSV, dendrogram PNG

```
python clustering_script_ML.py input.afa
```
Cluster sequences using Maximum Likelihood.

Input: aligned FASTA

Output: cluster CSV, dendrogram PNG

**UPGMA / NJ Trees**
```
python tree_with_clusters.py --tree MUSCLE_UPGMA_newick.nwk --colors cluster_colors.tsv --k 20 --out_prefix UPGMA_tree
python tree_with_clusters.py --tree MUSCLE_NJ_newick.nwk --colors cluster_colors.tsv --k 20 --out_prefix NJ_tree
```
Build UPGMA or NJ trees, color by cluster, and export.

Input: Newick tree + colors TSV

Output: cluster CSV, dendrogram PNG

### 7. Clustering (Centrolign + Percent Identity)
```
source /private/groups/migalab/bmahi/ModDotPlot/mdpenv/bin/activate
```
Activate ModDotPlot/centrolign conda environment.

```
centrolign_automater.sh input.fa output_dir
```
Run pairwise alignments for all sequences in a FASTA.

Input: FASTA

Output: alignment TXT files

```
python calculate_centrolign_identity.py
```
Calculate percent identity from centrolign alignments.

Input: alignment TXT files

Output: percent identity TXT

```
python perid_clustering.py
```
Cluster sequences based on pairwise percent identity.

Input: percent identity TXT

Output: cluster CSV, cluster TXT, dendrogram PNG

### 8. Plotting Alignments
**Fedor’s plot**
```
python3 plot_alignment.py -a aligned.fa
```
Visualize alignment as heatmap-like plot.

Input: aligned FASTA

Output: PNG

**Cluster-based alignment plots**
```
python ../../plot_alignment_COPY.py --alignment aligned.fa --pid_table pairwise_percent_identity.txt --output_prefix prefix --max_clusters 20 --verbose
```
Plot alignment with cluster coloring (max 20 clusters).

Input: aligned FASTA + percent identity TXT

Output: PNG

```
python ../../plot_alignment_COPY_2.py --alignment aligned.fa --clusters clusters.csv --colors colors.tsv --sortby clusters --output output.png
```
Plot alignment sorted by cluster assignment.

Input: aligned FASTA + clusters CSV + colors TSV

Output: PNG

```
python ../../plot_alignment_COPY_2_wDendrogram.py --alignment aligned.fa --clusters clusters.csv --colors colors.tsv --sortby clusters --pid_table pairwise_percent_identity.txt --output_prefix output_prefix
```
Plot alignment + dendrogram on y-axis.

Input: aligned FASTA + clusters CSV + colors TSV + percent identity TXT

Output: PNG

**Positional plots**
```
python3 plot_alignment_copy.py --alignment aln.fa --output aln.png --sortby coords
```
Plot alignment ordered by genomic coordinates.

Input: aligned FASTA

Output: PNG

**With strand bar**
```
python3 plot_alignment_str.py --alignment aln.fa --output aln_str.png --sortby coords --bed input.bed
```
Plot alignment with strand-orientation bar.

Input: aligned FASTA + BED

Output: PNG

### 9. NTRprism (Higher-Order Repeat Discovery)
```
perl ../../NTRprism_ProcessFasta_v0.3e.pl input.fa defaultParams 20000 5 6 1 0 0
```
Run NTRprism on FASTA with standard span.

Input: FASTA

Output: TXT (k-mer tables, TopHits, bin files)

```
perl ../../NTRprism_ProcessFasta_v0.3e.pl input.fa defaultParams 500000 5 6 1 0 0
```
Run NTRprism with extended span (covers full arrays).

Input: FASTA

Output: TXT (k-mer tables, TopHits, bin files)

**Plot spectra**
```
Rscript NTRprism_PlotSpectrum.r --args input.txt input.txt 20000 output_spectrum
```
Plot repeat spectra.

Input: TXT from NTRprism

Output: PNG/PDF

### 10. ModDotPlot
```
moddotplot interactive -f input.fa
```
Interactive dotplot visualization of FASTA sequences.

Input: FASTA

Output: PNG/SVG

### 11. Alignment Trimming (BMGE)
```
conda create -c bioconda bmge -n bmge
conda activate bmge
```
Install and activate BMGE environment.

```
bmge -i input.aligned.fa -t DNA -of output.bmge.fa -oh output.bmge.html
```
Trim alignment with BMGE defaults.

Input: aligned FASTA

Output: trimmed FASTA, HTML report

```
bmge -i input.fa -t DNA -h 0.8 -of output.fa -oh output.html
```
Relax entropy threshold to keep more sites.

Input: aligned FASTA

Output: trimmed FASTA, HTML report

```
bmge -i input.fa -t DNA -h 0.8 -m DNAPAM250:4 -of output.fa -oh output.html
```
Trim with permissive DNA substitution matrix.

Input: aligned FASTA

Output: trimmed FASTA, HTML report

### 12. Miscellaneous Utilities
```
awk '{print $3 - $2}' input.bed | sort -n
```
Print BED region lengths.

Input: BED

Output: sizes (stdout)

```
awk '($3 - $2) > 445' input.bed > filtered.bed
```
Filter out short BED regions.

Input: BED

Output: BED

### 13. Bedtools (Sequence Extraction & Utilities)

annotate → annotate coverage of features

closest → find nearest interval

cluster → group nearby intervals

getfasta → extract FASTA sequences from BED

intersect → find overlapping intervals

nuc → compute nucleotide content

Input: BED + FASTA

Output: FASTA / BED / summary tables

### 14. Gepard (Dotplots)

Installed and run locally with GUI; parameters set in interactive window.

Parameters: `Word length: 75`, `Window size: 0`, `Matrix: DNA`

Input: FASTA

Output: PNG/SVG dotplots

### 15. Centrolign Alignments

Executable:
`/private/groups/migalab/bmahi/centrolign/build/centrolign`

Config file:
`/private/groups/migalab/bmahi/centrolign/build/config.yaml`

Parameters edited each run:

- `minimum_segment_score: 6500`

- `max_count: 1500`

- `max_num_match_pairs: 625000`

- `constraint_method: 0`

Input: FASTA

Output: alignment TXT files

### 16. Plotting Alignments with Dotplots

Command:

`python plot_dotplot_alignment.py <fasta1> <fasta2> <alignmentFile> <svg_outfile>`

Parameters:

- `min_length: 100`

- `mem_line_width: 3.5`

- `mum_line_width: 3.5`

Input: FASTA + Centrolign alignment

Output: SVG (dotplot with alignment overlay)

### 17. ModDotPlot (Interactive Dotplots)

Run in interactive mode:

`moddotplot interactive -f <combinedFASTA> --port 8080 --compare-only`

Port forwarding (run locally in separate terminal):

`ssh -L 8080:127.0.0.1:8080 bmahi@emerald`

Input: FASTA

Output: Interactive browser-based dotplot

--------------------------------------------------------------------

The methods below this point relate to HSat analysis.

### 18. Run Centrolign on HSat Arrays
`centrolign --config config.yaml > output.cigar.txt`

Run centrolign pairwise alignment with parameters defined in config.

Input: FASTA, config.yaml

Output: CIGAR alignment TXT

### 19. Plot Dotplots with Alignments
```
python3 ../plot_dotplot_alignment_COPY.py \
    chrArrayFASTA \
    chrArrayFASTA \
    cigar \
    output.svg \
    1
```

Generate dotplot with alignment overlay.

Input: two FASTA files + CIGAR alignment

Output: SVG

### 20. Extract Sequences from Multi-FASTA
`awk '/^>sampleHeader/ {p=1; print; next} /^>/ {p=0} p' chr9_HSat3_bSat.merged.fa > output.fa`

Extract one sample (by header) into its own FASTA file.

Input: multi-FASTA

Output: FASTA

### 21. Pairwise Alignment of Sample FASTAs
`./centrolign --config config.yaml > output.cigar.txt`

Run centrolign on two combined sample FASTAs.

Input: pairwise FASTA, config.yaml

Output: CIGAR alignment TXT

```
python3 ../plot_dotplot_alignment_COPY.py \
    sampleFASTA1 \
    sampleFASTA2 \
    cigar \
    output.svg \
    1
```

Plot pairwise alignment between samples.

Input: 2 FASTA + CIGAR

Output: SVG

### 22. Prune Tree for Clade Analysis
```
nw_prune -v -f \
  HPRC_chr9_48354832_49055551_76500000_76850000_het141_m_hprc_dgp_upgma.nwk \
  ../FAs/clade1_samples.txt \
> pruned_clade1_HPRC_chr9_48354832_49055551_76500000_76850000_het141_m_hprc_dgp_upgma.nwk
```

Prune Newick tree to a subset of clade samples.

Input: Newick tree + list of samples

Output: pruned Newick tree

### 23. Generate GFA from Centrolign
`./centrolign --config config.yaml > ../../GuideTrees/clade1_msa.gfa`

Produce a GFA graph alignment file for clade sequences.

Input: FASTA, config.yaml

Output: GFA