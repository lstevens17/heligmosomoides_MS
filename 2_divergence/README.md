# Section 2: The H. bakeri and H. polygyrus genomes are highly divergent 

[**barcode_trees**](<https://github.com/lstevens17/heligmosomoides_MS/tree/main/2_divergence/barcode_trees>)

FASTA files, alignments, and newick files for the ITS2 and COI trees. The divergence between each sequence is caculated using a python script (dist_matrix.py), which ignore gaps calculate identity as number of identical sites divide by the total number of aligned sites (this is analagous to the sequence identities calculated by clustal omega).

[**busco_minimap2_divergence**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/2_divergence/busco_minimap2_divergence)

Contains code and files to plot divergence between H. bakeri and H. polygyrus and between two Caenorhabditis sister species pairs. BED file contains minimap2 alignments of the H. polygyrus genome to the H. bakeri reference genome and the nucleotide divergence. This file is plotted using an R script (plot_nuc_aa_alignment.R). TSV files contain   to the H. bakeri reference genome and the TSV files contain sequence IDs for each reference species (H. bakeri, C. briggsae, and C. remanei, respectively), their location, and their protein sequence identity (calculated using dist_matrix.py). These TSVs are plotted using two R scripts (plot_nuc_aa_alignment.R and plot_comparative_aa_divergence.R).

[**synteny**](https://github.com/lstevens17/heligmosomoides_MS/tree/main/2_divergence/synteny)

Contains code and files to generate the Oxford plot show gene order between H. bakeri and H. polygyrus. A TSV contains the locations of BUSCO genes in each species along with the Nigon assignment. This TSV is plotted by an R scipt (plot_nigon_Oxford.R).

