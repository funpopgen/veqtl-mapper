#VEQM :  Mapping genetic variance expression quantitative trait loci.

##Introduction:

TODO
##NAME

VEQM  -  Mapping genetic variance expression quantitative trait loci.
##SYNOPSIS

VEQM [options]
##DESCRIPTION

TO DO
##OPTIONS

SYNOPSIS
       VEQM [options]

OPTIONS

       --help    Display help information.

       --version Display version information.

       --bed CHAR
                 Phenotype file [last argument].

       --vcf CHAR
                 Genotype file.

       --out, --o CHAR
                 Output file [stdout].

       --job-number INT
                 Split  the  analysis  into  a  number of smaller runs which can run in parallel on a cluster. This option
                 specifies which of the sub-analyses should be run.

       --genes INT
                 This specifies the number of genes to be analysed in each job.

       --window INT
                 The size in base pairs of the cis window around the transcription start site of the gene [1,000,000].

       --perm INT(,INT)
                 Calculated permuted p values, one following number indicates the number of permutations, two comma  sepa‐
                 rated numbers gives the number of permutations and the seed.

       --het     Tests for variance differences in the heterozygous relative to homozygous groups, an indication of parent
                 of origin effects.

       --bcftools
                 Use the bcftools plugin to calculate dosage from the GT field.

       --noheader
                 Suppress writing of header line.

       --nocheck Do not use the header lines to match genotype and phenotype, assume samples already match.

##FILE FORMATS

   INPUT FILE FORMATS

       Phenotype
              BED file.

       Genotype
              VCF format. Format should be one DS field unless bcftools flag is given. Dosages
			  should  be either  0, 1 or 2. If the bcftools flag is given, the bcftools plugin
			  dosage will be used to calulate dosage from GT field.

   OUTPUT FILE FORMAT

       Output contains the chromosome, location, reference and  alternate  alleles  of  the
       SNP, followed by the correlation, P value, Permutation  P  value,  then adjusted for 
       multiple testing P values by permutations and beta approximation.

VEQM-1.0.0                              8th March 2015                               VEQM(1)
