# Methods
## Bhavna Mahi
## Miga Lab

This README details every file on my server account including data, scripts (and their functions), and results. 

Top-level directory: `/private/groups/migalab/bmahi`

1. bsatLong12 - This directory contains fasta files of the longest 12 beta satellite arrays across CHM13 version 2.0.
    - 12 different fasta files

-------------------------------------------------------------------------

2. centrolign - This is an original copy of Jordan's centrolign tool copied over on 01/01/2025. There have been changes to the script since the copy date.
    - build
        - **centrolign** - Main executable file
        - CMakeCache.txt
        - CMakeFiles
        - cmake_install.cmake
        - **config.yaml** - Configuration file
        - libcentrolign.so
        - Makefile
        - version.cpp
    - CMakeLists.txt
    - include
    - LICENSE
    - README.md
    - src

Note: I didn't really mess with anything past `centrolign/build` so I did not list all the contents of the subdirectories here.

-------------------------------------------------------------------------

3. centrolign_automater.sh - Automates pairwise alignment using Centrolign for all regional FASTA files and aligns each region against the full chromosome 1 beta satellite array.

-------------------------------------------------------------------------

4. chm13v2.0.fa - The entire CHM13 version 2.0 genome.

-------------------------------------------------------------------------

5. chm13v2.0.fa.fai - The entire CHM13 version 2.0 genome indexed.

-------------------------------------------------------------------------

6. chr1_bsat_project - This directory contains all data and scripts related to beta satellite analysis on chromosome 1. 
    - alignmentPlots
        - arrayFASTAS - Contains fasta files for all the beta satellite array regions I identified (14 total). Also contains scripts to automate the process of running ModDotPlot on all of these pairs of fastas. 
        - calculate_centrolign_identity.sh - Computes percent identity using only match counts from centrolign output.
        - centrolign_alignments - Contains all the cigar output of all array regions aligned against each other.
        - chr1BsatArraySelf - Contains bed file for entire array and versions of self/self alignment plots.
        - chr1BsatArrayvReg - Contains several subdirectories which all contain alignment data for each array region versus the entire array itself.
        - invertedRegPlots - Contains several subdirectories which all contain alignment data for each array region versus itself and each other array region.
        - percentIdentity.txt - Text file that lists the percent identity calculations of each array region against each other.
    - arrayRegions.bed - BED file of the 14 array regions I determined.
    - chm13v2.0_chr1.fa - CHM13 version 2 entire chromosome 1 fasta.
    - chm13v2.0_chr1.fa.fai - CHM13 version 2 entire chromosome 1 fasta indexed.
    - hg002v1.1_mat_chr1.fasta - HG002 maternal entire chromosome 1 fasta.
    - hg002v1.1_mat_chr1.fasta.fai - HG002 maternal entire chromosome 1 fasta indexed.
    - hg002v1.1_pat_chr1.fasta - HG002 paternal entire chromosome 1 fasta.
    - hg002v1.1_pat_chr1.fasta.fai - HG002 paternal entire chromosome 1 fasta indexed.
    - hmmertblout2bed.awk - Script that converts HMMER table output to a BED file.
    - mfa2cons.py - Script that I copied over from Fedor that finds a consensus sequence between a multiple sequence alignment.
    - model1 - Contains all data pertaining to model 1 run on chromosome 1 beta satellite array.
        - chr1bsatModel.hmm - First model hmm file.
        - chr1_bsat_repeat.bed - BED file for the model 1 repeat unit.
        - chr1_bsat_repeat.fa - Fasta file for the model 1 repeat unit.
        - chr1bsatResults.bed - BED file of HMMER output from model 1 run.
        - chr1bsatResults.out - HMMER table output from model 1 run.
    - model2 - Contains all data pertaining to model 2 run on chromosome 1 beta satellite array.
        - all_subclusters.bed - BED file of all subclusters identified in this model. 
        - chr1bsatModel2.hmm - Second model hmm file.
        - chr1_bsat_repeat_2.bed - BED file for the model 2 repeat unit.
        - chr1_bsat_repeat_2.fa - Fasta file for the model 2 repeat unit.
        - chr1bsatResults2.bed - BED file of HMMER output from model 2 run.
        - chr1bsatResults2_filtered_aligned_alignment_plot.png - Alignment plot of size filtered model 2 results.
        - chr1bsatResults2_filtered_aligned.fa - All the model 2 run output seqs filtered by size and aligned.
        - chr1bsatResults2_filtered.bed - BED file of all the model 2 output seqs filtered by size. 
        - chr1bsatResults2_filtered.fa - Fasta file of all the model 2 output seqs filtered by size. 
        - chr1bsatResults2.out - HMMER table output from model 2 run.
        - peridClusters - Directory with all data pertaining to model 2 results clustering by percent identity. 
        - sizeClusteres -  Directory with all data pertaining to model 2 results clustering by sequence size. 
    - model3 - Contains all data pertaining to model 3 run on chromosome 1 beta satellite array.
        - chr1bsatModel3.hmm - Third model hmm file.
        - chr1_bsat_repeat_3.bed - BED file for the model 3 repeat unit.
        - chr1_bsat_repeat_3.fa - Fasta file for the model 3 repeat unit.
        - chr1bsatResults3.bed - BED file of HMMER output from model 3 run.
        - chr1bsatResults3.out - HMMER table output from model 3 run.
    - model4 - Contains all data pertaining to model 4 run on chromosome 1 beta satellite array.
        - chr1bsatModel4.hmm - Fourth model hmm file.
        - chr1_bsat_repeat_2.fa - Fasta file for the model 2 repeat unit.
        - chr1_bsat_repeat_BSR.fa - Fasta file for the signature BSR repeat unit.
        - chr1bsatResults4.bed - BED file of HMMER output from model 4 run.
        - chr1bsatResuts4.out - HMMER table output from model 4 run.
        - combined_seqs.fa - Combined fasta with model 2 repeat unit and signature BSR repeat unit.
    - model5 - Contains all data pertaining to model 5 run on chromosome 1 beta satellite array.
        - alignedFASTA_PIcalc.py - Calculates percent identity between aligned sequences.
        - bsr_sequence_lengths.png  - BSR sequence lengths from model 5 output.
        - calculate_centrolign_identity.py - Calculates percent identity of sequence pairs from centrolign alignment output.
        - centrolign_automater.sh - Copied over centrolign automation script.
        - chr1bsatModel5.hmm - Fifth model hmm file.
        - chr1bsatResults5.bed - BED file of HMMER output from model 5 run.
        - chr1bsatResults5.fa - Fasta file of HMMER output seqs from model 5 run.
        - chr1bsatResults5_filtered.bed - BED file of HMMER output from model 5 run filtered by size.
        - chr1bsatResults5_filtered.fa - Fasta file of HMMER output from model 5 run filtered by size.
        - chr1bsatResults5_filtered_rmdup.fa - Fasta file of HMMER output from model 5 run filtered by size and duplicate sequences removed.
        - chr1bsatResults5.out - Hmmer table output from model 5 run.
        - chr1_bsr_repeat_5.fa - Fasta file for model 5 repeat unit. 
        - model5_ClustalOmegaAlignment - Directory with all data pertaining to model 5 output sequence alignment with Clustal Omega alignment tool. 
        - model5_MAFFTalignment - Directory with all data pertaining to model 5 output sequence alignment with MAFFT alignment tool. 
        - model5_MUSCLEalignment - Directory with all data pertaining to model 5 output sequence alignment with MUSCLE alignment tool. 
        - model5_NEW_MUSCLEalignment - Directory with all data pertaining to model 5 output sequence alignment with MUSCLE alignment tool but redone on select sequences. 
            - **align_and_tree.py** - Uses MUSCLE alignment to try and create a tree of sequences (janky!)
            - **clustering_script** - Attempts at trying to create clusters using different methods from MUSCLE alignment.
        - model5_NTRPrism - Directory with all data pertaining to NTR Prism runs with model 5 output sequences.
        - model5_PosAln - Directory with all data pertaining to positional alignment attempts with model 5 output sequences.
            - **regional_switch_logo_plot.py** - Attempt at trying to create a logo plot that shows dominant nucleotide at switch sites on array (janky!).
        - perid_clustering_Clustal.py - Script that clusters sequences from Clustal Omega alignment based on percent identity (janky!).
        - perid_clustering.py - Script that clusters sequences from an alignment based on percent identity (janky!).
        - perid_clusters_to_bed.py - Generates a new BED file of sequences with color codes based on each sequence's cluster.
        - plot_bsat_lengths.py - Creates histogram of sequence sizes (used to create bsr_sequence_lengths.png).
    - plot_alignment_copy.py - A copy of Jordan's plotting script with adjustments for specific tasks.
    - plot_alignment_original.py - Another copy of Jordan's plotting script with adjustments for specific tasks.

    Note: There is a lot of data in my model 5 analysis. A lot of it is redundant or negligable. Due to this reason, I have not listed the entirety of the content in this directory, but all of the analysis that is most crucial. 
-------------------------------------------------------------------------

7. chr21_chr22_bsat1 - This directory contains all data and scripts related to beta satellite analysis on chromosomes 21 and 22. 
    - chr21_bsat1_self.fa
        - Self-self fasta file of chromosome 21 beta satellite 1.
    - chr21_bsat1_seq.fa
        - Chromosome 21 beta satellite 1 fasta. 
    - chr21_bsat1_seq.fa.fai
        - Chromosome 21 beta satellite 1 fasta indexed.
    - chr21_chr22_centrolign.txt
        - Chromosome 21 and 22 beta satellite 1 centrolign cigar.
    - chr21_chr22_combined.fa
        - Chromosome 21 and 22 beta satellite 1 combined fasta.
    - chr21_chr22_combined.fa.fai
        - Chromosome 21 and 22 beta satellite 1 combined fasta indexed.
    - chr21_chr22_dotplot.svg
        Chromosome 21 and 22 beta satellite 1 alignment dotplot.
    - chr21_self_centrolign.txt
        - Chromosome 21 self alignment cigar.
    - chr21_self_dotplot.svg
        - Chromosome 21 self alignment dotplot.
    - chr22_bsat1_self.fa
        - Self-self fasta file of chromosome 22 beta satellite 1.
    - chr22_bsat1_seq.fa
        - Chromosome 22 beta satellite 1 fasta. 
    - chr22_bsat1_seq.fa.fai
        - Chromosome 22 beta satellite 1 fasta indexed.
    - chr22_self_centrolign.txt
        - Chromosome 22 self alignment cigar.
    - chr22_self_dotplot.svg
        - Chromosome 22 self alignment dotplot.

-------------------------------------------------------------------------

8. hsat_alignment - This directory contains all data and scripts related to human satellite analysts. 
    - BEDs - Contains all the BED files for chromosome 9 and 10 hsat arrays across CHM 13 version 2 and HG002 maternal and paternal copies.
    - centrolignCOPY - Just a copy of `/private/groups/migalab/bmahi/centrolign` with different parameters set in the config.yaml.
    - centrolignOutput - This directory contains all the txt and svg files outputted from centrolign alignments of all vs all hsat arrays from CHM13 and HG002 haplotypes.
    - chr9_chr10_FAs - This directory contains all the fasta files (both indexed and non-indexed) of the target chromosome 9 and 10 HSat arrays across CHM 13 version 2.0 and HG002 haplotypes.
    - chr9_satHap_alns - This directory contains HSat sat haps fastas that were aligned against each other using centrolign. 
        - clade1 - Contains copies of centrolign and the `centrolign_automater.sh` script.
            - FAs - Contains all fasta files of the sat haps we're targeting in this first 'clade'. 
            - GuideTrees - Contains newick files for sat hap guide trees and the gfa that was produced for clade 1.
    - FAs - This directory contains all the fastas for the entirety of chromosomes 9 and 10 on CHM13 version 2.0 and HG002 haplotypes. It also contains combined fasta files for centrolign runs.
    - GuideTrees - This directory contains a guide tree for chr 10.
    - hsat_2_3_analysis - This directory contains all data and results for an HSat 2 and 3 analysis on chromosome 10 on CHM13 version 2 and HG002 haplotypes.
        - BEDs - Contains bed files for target HSat arrays.
        - centrolignCOPY - Copy of centrolign with modified parameters for this run. 
        - FAs - Contains fasta files for target HSat arrays.
        - plot_dotplot_alignment_COPY.py - Copy of Jordan's alignment ploting script for this analysis. 
    - plot_dotplot_alignment_COPY.py - Copy of Jordan's alignment ploting script.

-------------------------------------------------------------------------

9. ModDotPlot - This is the ModDotPlot package to run alignments. Not sure which version this is. 
    - build
    - CITATION.cff
    - config
    - images
    - LICENSE
    - **mdpenv** - Virtual environment I used to run ModDotPlot.
    - pyproject.toml
    - README.md
    - sequences
    - src

Note: I ran ModDotPlot in my virtual environment without making any modifications to the config yaml so I will not list the contents of these subdirectories.

-------------------------------------------------------------------------

10. plot_dotplot_alignment.py - This is an original copy of Jordan's dotplot alignment script copied over on 01/01/2025. I'm not sure if he's made changes to the script since then. 