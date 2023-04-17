# set parameters
id=95
flen=2000
len=10000
dens=0.002

# find alignment gaps in all
for file in *filtered.bed; do echo "bedtools complement -i <(awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$3-\$2}' $file | awk '\$4 >= $id && \$5 >= $flen' | sort -k1,1 -k2,2n) -g <(sort -k1,1 -k2,2n nxHelBake1.1.primary.fa.fai) >`echo $file | cut -f1 -d'_'`_${id}_${flen}_alignment_gaps.bed"; done | bash

# generate coverage files
for file in *gaps.bed; do prefix=`echo $file | cut -f1 -d'_'`; bedtools coverage -a $file -b ${prefix}_vs_ref.sorted.biallelic.vcf | bedtools coverage -a - -b <(awk '$7 > 0' ${prefix}_vs_ref.sorted.biallelic.vcf_1kb_bins.bed) | bedtools coverage -a - -b <(awk '$4 <= 100' ${prefix}_vs_ref.filtered.bed) | cut -f1,2,3,7,11,15 | awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4"\t"$5"\t"$6}' >`echo $file | sed 's/.bed/.coverage.bed/g'`; echo $file; done

# call divergent haplotypes
for file in *coverage.bed; do echo "awk '\$4 >= $len && \$5 >= $dens && \$6 >= 0.3 && \$7 >= 0.3' $file >`echo ${file} | sed 's/.coverage.bed//'`_${len}_${dens}_0.3_0.3_divergent.bed"; done | bash

# remove X calls in alternate assemblies
for file in *alternate*divergent.bed; do grep -v X $file >tmp && mv tmp $file; done

# remove non-chromosomal scaffolds
for file in *divergent.bed; do grep -v scaffold $file | grep -v unloc >tmp && mv tmp $file; done

# get divergent genes 
for file in *divergent.bed; do bedtools intersect -a <(grep -P "\tmRNA\t" nxHelBake1.1.primary.final_annotations.gff3_longest_isoforms) -b $file -f 0.5 | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' >${file}_genes.tsv; echo $file; done
