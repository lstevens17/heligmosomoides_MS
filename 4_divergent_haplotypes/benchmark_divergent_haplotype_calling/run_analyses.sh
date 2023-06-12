# get alignment gaps
cat strains.txt | while read strain; do for id in 94 94.5 95 95.5 96 96.5 97; do for len in 2000 3000 4000 5000 6000; do echo "bedtools complement -i <(cat ${strain}_vs_ref.filtered.bed | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$3-\$2}' | awk '\$4 >= $id && \$5 >= ${len}' | sort -k1,1 -k2,2n) -g <(sort -k1,1 -k2,2n c_elegans.PRJNA13758.WS286.genomic.fa.fai) >${strain}_alignment_gaps_${id}_${len}.bed"; done; done; done >bedtools_complement_commands.txt
parallel -j 128 :::: bedtools_complement_commands.txt

# get cov files (takes some time)
cat strains.txt | while read strain; do for file in ${strain}*alig*gaps*bed; do bedtools coverage -a $file -b ${strain}_vs_ref.sorted.biallelic.vcf | bedtools coverage -a - -b <(awk '$7 > 0' ${strain}_vs_ref.sorted.biallelic.vcf_1kb_bins.bed) | bedtools coverage -a - -b <(awk '$4 <= 95' ${strain}_vs_ref.filtered.bed) | cut -f1,2,3,7,11,15 | awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4"\t"$5"\t"$6}' >`echo $file | sed 's/.bed/.coverage.bed/g'`; echo $file; done; done

# filter at various values 
cat strains.txt | while read strain; do for len in 5000 10000 15000; do for dens in 0.001 0.002 0.003 0.004 0.005; do for file in ${strain}*coverage.bed; do echo "awk '\$4 >= ${len} && \$5 >= ${dens} && \$6 >= 0 && \$7 >= 0' $file >`echo ${file} | sed 's/.coverage.bed//'`_${len}_${dens}_0_0_divergent.bed"; done; done; done; done >filter_commands.txt
parallel -j 128 :::: filter_commands.txt

# for each divergent region call, get total number of bases classified as divergent, the number of overlapping bases between prev calls and total number of bases in previous call
cat strains.txt | while read strain; do for file in ${strain}*divergent.bed; do echo -e "$file\t`cut -f4 $file | paste -sd+ | bc -l`\t`bedtools coverage -a <(cut -f1,2,3 $file) -b ${strain}_divergent_regions.bed | cut -f5 | paste -sd+ | bc -l`\t`awk '{print $3 - $2}' ${strain}_divergent_regions.bed | paste -sd+ | bc -l`"; done; done >benchmark_values.tsv

# combine to summary for each strain
cut -f1 benchmark_values.tsv | cut -f2- -d'_' | sort | uniq  | while read param; do echo -e $param"\t"`grep $param data_matrix.tsv | cut -f2 | paste -sd+ | bc -l`"\t"`grep $param data_matrix.tsv | cut -f3 | paste -sd+ | bc -l`"\t"`grep $param data_matrix.tsv | cut -f4 | paste -sd+ | bc -l`; done >benchmark_values.combined.tsv
