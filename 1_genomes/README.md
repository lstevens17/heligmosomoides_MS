# Section 1: Chromosome-level reference genomes for H. bakeri and H. polygyrus

[**cumulative_length**](<https://github.com/lstevens17/heligmosomoides_MS/tree/main/1_genomes/cumulative_length>)

Contains code to generate cumulative length curves shown in Figure 1B. TSVs for each assembly (three single worm H. bakeri assemblies, two H. polygyrus assemblies, and the two previous H. bakeri reference genomes) contain the order of contigs (ordered by descending length) and the length of each contig. These values are plotted using an R script (plot_cumulative_lengths.R).

[**nigon_hcontortus**](<https://github.com/lstevens17/heligmosomoides_MS/tree/main/1_genomes/nigon_hcontortus>)

Contains code and files to generate Nigon and _H. contortus_ paintings shown in Figure 1D and Figure S3. The locations of BUSCOs and their Nigon/_H. contortus_ labels are found in TSV files, which are plotted using an R script (plot_nigon_distribution.R).

[**strongylomorph_phylogeny**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/1_genomes/strongylomorph_phylogeny)

Contains the concatenated alignment (supermatrix.fa) and the resulting newick file (derived from ASTRAL, with branch lengths from IQ-TREE) associated with the strongylomorph phylogeny shown in Figure 1A. 

[**telomeres**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/1_genomes/telomeres)

Contains code and files to generate telomere density plots shown in Figure S4. TSVs contain density of nematode telomeric repeats (TTAGGC) in non-overlapping 1kb windows derived from seqkit locate. These are plotted using an R script (plot_telomere_density.R). 

[**gc_coverage_bias**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/1_genomes/gc_coverage_bias)

Contains code and files to generate the GC coverage files in Figure S5. BED files containing GC% and coverage are plotted using an R script (plot_nxHelBake1_ngHelPoly1_gc_coverage_bias.R).