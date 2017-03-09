#!/usr/bin/env bash

# Download the 1000 Genomes autosome vcf files from the website
for i in {1..22}; do
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done

# Download the expression data

wget ftp://jungle.unige.ch/veqtl_paper/geuvadis.linc.protein.cov.50PC.bed

# Subset the vcf files to individuals with genotype data, select variants with MAF > 0.05 and concatenate all chromosomes together
bsub -o out.vcf "bcftools concat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Ou | bcftools view --force-samples -S Samples -Ou | bcftools view --min-af 0.05:minor -Oz -o Genotypes.individuals_filtered.maf0.05.vcf.gz && bcftools index Genotypes.individuals_filtered.maf0.05.vcf.gz"

# Delete the original files
rm ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Run the v-eQTL mapping, split into 490 jobs, each analysing 29 genes ### MAY NEED TO ADAPT THIS LINE ACCORDING TO CLUSTER TYPE
bsub -o out.veqtl -J"veqtl[1-490]" "veqtl_mapper --bed geuvadis.linc.protein.cov.50PC.bed --vcf Genotypes.individuals_filtered.maf0.05.vcf.gz --genes 29 --job-number \$LSB_JOBINDEX --verbose --perm 10000,4 --out results.veqtl\$LSB_JOBINDEX"
##339.17 CPU hours

# Concatenate the 490 results files, keeping the header
awk 'FNR>1||NR==1' results.veqtl* > veqtl.results

# Extract the most significant asssociation per gene
(head -1 veqtl.results; tail -n+2 veqtl.results |  awk '{if (best[$1]=="" || p_best[$1] > $7) {best[$1] = $0; p_best[$1]=$7}} END{for (var in best) print best[var]}') > veqtl.results.best

# Delete the original results
rm results.veqtl*

# Filter the genotype file to only contain variants with at least 50 individuals in each of the three genotype groups for parent of origin using python script
bsub -o out.filter python3 filter_snps.py

# Repeat steps for parent of origin analysis
bsub -o out.het -J"het[1-490]" "veqtl_mapper --bed geuvadis.linc.protein.cov.50PC.bed --vcf Genotypes.individuals_filtered.parent.of.origin.vcf.gz --genes 29 --job-number \$LSB_JOBINDEX --verbose --het --perm 10000,4 --out results.poeqtl\$LSB_JOBINDEX"
##101.21 CPU hours

awk 'FNR>1||NR==1' results.poeqtl* > poveqtl.results

(head -1 poveqtl.results; tail -n+2 poveqtl.results |  awk '{if (best[$1]=="" || p_best[$1] > $7) {best[$1] = $0; p_best[$1]=$7}} END{for (var in best) print best[var]}') > poveqtl.results.best

rm results.poeqtl*

# These results files can be analysed to produce the results and Figures from the paper using the R script Figures.R 
