# HSat Alignment with centrolign
### Name: Bhavna Mahi
### PI: Dr. Karen Miga
### Mentor: Hailey Loucks
### Institution: UCSC Genomics Institute - Miga Lab
### Working Directory: `bmahi@emerald:/private/groups/migalab/bmahi/hsat_alignment`
-------------------------------------------------------------------------------------------------

## Project Summary
Try running centrolign on hsat arrays to see which are alignable. 

## Mentor Meeting (05/23/2025)
- Try putting 200 seqs into Clustal Omega and create Fedor plot
- In RepeatMasker the repeat monomer for BSR is a dimer
    - Known SNP in original BSR unit definition (Waye Willard paper)
    - See if Waye/Willard SNP at position 25 exists on alignment plots
- Double check that alignments are using MUSCLE alignment file
    - Subset based on variants and see if Fedor plot aligns better
- In NTR TopHits of BSR seqs
    - We see 68 hit, then 137, then variable number
- Plot spectra for BSR TopHits
- HSat centrolign
    - HSat3 on HG002 on chr9
    - Run centrolign on 2 haplotypes first and then on the entire HPRC
    - Might be able to get graph from Sasha to find 2 closely related haplotypes
    - Try it first on arrays on CHM13 and HG002
    - Create alignment plot with alignment line going through (Jordan's script)
    - Hailey will send chr9 bed files for hsat3 for CHM13 and HG002 and another smaller array that should be alignable
- Be prepared to discuss BSats or HSats at Wed CenSat meeting

## Run centrolign on chr9/chr10 hsat3 arrays of chm13 and hg002
- Here are the arrays for chm13 and hg002 on chr9:
    - chr9: 48476975-76694045 (28,717,070 bp)
    - chr9_MATERNAL: 48305096-69801493 (21,496,397 bp)
    - chr9_PATERNAL: 46250756-57492248 (11,241,492 bp)
- Here are the arrays for chm13 and hg002 on chr10:
    - chr10: 43171895-43199449 (27,554 bp)
    - chr10_PATERNAL: 44205188-44232421 (27,233 bp)
    - chr10_MATERNAL: 43608334-43635884 (27,550 bp)
- First create paired FASTA files
  - CHM13 + HG2_Mat: chm13_hg002_mat_chr9.fasta and chm13_hg002_mat_chr10.fasta
  - CHM13 + HG2_Pat: chm13_hg002_pat_chr9.fasta and chm13_hg002_pat_chr10.fasta
  - HG2_Mat + HG2_Pat: hg002_mat_pat_chr9.fasta and hg002_mat_pat_chr10.fasta
- Run centrolign on all pairwise alignment files: `centrolign --config config.yaml > output.cigar.txt`
- Run plotting script:
```
python3 ../plot_dotplot_alignment_COPY.py \
    chrArrayFASTA \
    chrArrayFASTA \
    cigar \
    output.svg \
    1
```
/private/groups/migalab/bmahi/hsat_alignment/hsat_2_3_analysis/FAs/chm13_hg002_mat_hsat_2_3.fasta
## Run centrolign on chr10 HSat2 + HSat3 array
- For the chromosome 10 arrays that we have sathap trees for here are the array coordinates:
    - chm13 chr10: 42,063,500-42,985,858 (922,358 bp)
    - HG002 chr10_MATERNAL: 42,578,486-43,418,931 (840,445 bp)
    - chr10_PATERNAL: 42,352,179-44,023,315 (1,671,136 bp)
- First create paired FASTA files
- Run centrolign on all pairwise alignment files: `centrolign --config config.yaml > output.cigar.txt`
- Run plotting script:
```
python3 ../plot_dotplot_alignment_COPY.py \
    chrArrayFASTA \
    chrArrayFASTA \
    cigar \
    output.svg \
    1
```

# Mentor Meeting (05/30/2025) w/ Hailey + Julian L.
- Hsat project:
    - centrolign doesn't show the details of the alignment
    - We need to get closely related HSats and rerun this before we share any results
    - On the SatHap tree, find the closely related clades and run all vs all pairwise alignment within clades

# Mentor Meeting (06/06/2025):
- HSat stuff
    - Hailey will get me fasta files to be able to create gfas
    - Find closely related nodes to do pairwise alignment (look for the closest nodes on the trees)
        - Hailey will make list of interesting pairwise alignments with the HSats
        - Recreate bandage plots similar to Julian's
    - Analysis of cigar strings to find types of differences
        - Differences between more similar and more different haps
        - Calculate transitions vs transversions
        - Insertions vs deletions
        - Single nucleotide differences (what types of mutations are we seeing)
            - What areas are more likely to mutate
        - Run centrolign to get GFA
            - Analyze GFAs
                - Do we want intermediate GFAs???

# More pairwise alignments
- Using Sasha's chr9 sat tree
- Do pariwise alignment of these:
    - HG01943.1 and HG01192.2
    - HG02922.2 and HG02055.1
    - Maybe more complicated but the little clade with HG02965.1, HG01960.1, HG01884.2 and HG03540.2 could be interested to do pairwise on and try to estimate the accumulation of different types of vars between closer/more distant samples
- Steps: 
    1. Extract samples from Sasha's fasta: `awk '/^>sampleHeader/ {p=1; print; next} /^>/ {p=0} p' chr9_HSat3_bSat.merged.fa > output.fa`
    2. Combine each pair of samples into one fasta (replace all #s with _s)
    3. Update the centrolign config to the combined pair fasta file
    4. Run centrolign on the paired fasta: `./centrolign --config config.yaml > output.cigar.txt`
    5. Plot the pairwise alignment on the dot plot:
    ```
    python3 ../plot_dotplot_alignment_COPY.py \
    sampleFASTA1 \
    sampleFASTA2 \
    cigar \
    output.svg \
    1
    ```

# Mentor Meeting (06/20/2025)
- Take subset tree from Sasha's tree (pick interesting but small clade)
    - Do an alignment of all pairs in the tree (use centrolign automation script)
    - Do GFA of clade
        - Trim down tree for specific clade
        - Ask Mira about an easy way to parse the newick
    - In the pairwise alignments:
        - For every 100 mbp how many insertions do we find
        - Rate of transitions vs tranversions
        - Long term goal: discover which regions in the array are more susceptible to mutations
            - Look for tools that can analyze cigar string data
            - Parsing locations of events to map back onto the sequence
                - Look at Jordan's plotting script for parsing cigar tips

# Clade analysis on Sasha's tree
- Samples:
    1. HG02391.2*
    2. HG02071.1*
    3. HG04157.1*
    4. HG00609.1*
    5. NA18565.2*
    6. HG01943.1*
    7. HG01192.2*
    8. HG00639.2*
    9. HG01981.2*
    10. HG02015.1*
    11. HG03710.2*
    12. HG01109.1*
    13. HG00741.2*
    14. HG04228.2*
    15. HG03239.1*
    16. HG00658.2*
- Prune newick file: 
```
nw_prune -v -f \
  HPRC_chr9_48354832_49055551_76500000_76850000_het141_m_hprc_dgp_upgma.nwk \
  ../FAs/clade1_samples.txt \
> pruned_clade1_HPRC_chr9_48354832_49055551_76500000_76850000_het141_m_hprc_dgp_upgma.nwk
```
- Run centrolign to create gfa: `./centrolign --config config.yaml > ../../GuideTrees/clade1_msa.gfa`

# Mentor Meeting (06/27/2025)
- What are we doing with the gfa
    - Align samples back onto gfa
        - Try mapping sister clades to assess similarity
    - Read through vgteam documentation to learn how to map seqs back onto gfa
        - May need to create a psuedo-fastq 
        - Ask vg channel for help if needed