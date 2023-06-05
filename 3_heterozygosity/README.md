# Section 3: Heterozygosity is concentrated in distinct regions of the H. bakeri genome

[**benchmark_divergent_haplotype_calling**](<https://github.com/lstevens17/heligmosomoides_MS/tree/main/3_heterozygosity/benchmark_divergent_haplotype_calling>)

Contains code and files used to benchmark hyper-divergent haplotype calling parameters in C. elegans. The hyper-divergent haplotypes calls for C. elegans are found in the [GitHub](https://github.com/AndersenLab/Ce-328pop-div/blob/master/Processed_Data/divergent_classification.RData) associated with the Lee et al. 2021 manuscript. A Rscript (Ce328_extract_divergent_region_coordinates.R) loads those data and generates TSVs for each strain that has a PacBio assembly (\*_divergent_regions.tsv). BED formatted versions of these files (\*_divergent_regions.bed) are used to calculate specificity and sensitivity for hyper-divergent haplotype calls across a range of parameter values. 

[**divergent_regions**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/3_heterozygosity/divergent_regions)

Contains code and files used to call hyper-divergent haplotypes in H. bakeri. A bash script (run_analyses.sh) runs seven steps: 1) finds alignment gaps between alignments of a given length and percent identity using BEDtools complement, 2) calculates SNP density, coverage of 1kb bins that contain SNPs, and coverage of alignments of any identity in each alignment gap, 3) find hyper-divergent haplotypes using the specified parameters 4) removes hyper-divergent calls in X chromosome in alternate assemblies, 5) removes hyper-divergent calls in non-chromosomal scaffolds, 6) uses BEDtools to identify protein-coding genes that are overlapped 50% by a hyper-divergent haplotype.

The locations of the identified hyper-divergent haplotypes are plotted alongside paftools-derived SNP densities using an R script (plot_SNP_density_divergent_regions.R). 

[**go_term_enrichment**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/3_heterozygosity/go_term_enrichment)

Contains code and files for GO term enrichment using topGo. The IDs of transcripts/mRNAs derived from genes that are overlapped at least 50% by a hyper-divergent haplotype are provided to topGO (all_divergent_mRNA_IDs.txt) along with the GO terms derived from InterProScan. GO term enrichment is then run using all three ontologies using an R script (run_topgo.R). The enrichment results are output in one TSV file per ontology, summarised in Table S5 (tableS5_goterm.csv), and plotted using an R script (plot_go_terms.R).

[**snp_density**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/3_heterozygosity/snp_density)

Contains code and files to plot SNP density in three H. bakeri individuals and two H. polygyrus individuals relative to the reference genome. A BED file for each individual contains SNP density in 10 kb windows, which are plotted by two R scripts (plot_SNP_density_Fig3.R and plot_SNP_density_other_individuals.R).