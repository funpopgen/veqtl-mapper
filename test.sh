#!/bin/bash
if [[ $( ./bin/veqtl-mapper --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 | sha1sum | awk {'print toupper($1)'}) == "26F6F2F43BB3543FE7A8E060F0F9A088A74EACA1" ]]; then
    echo "Passed: variance test."
else
    echo "Failed: variance test."
    exit 1
fi

if [[ $( ./bin/veqtl-mapper --het --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 | sha1sum | awk {'print toupper($1)'}) == "600B468E56A76D97B2B6CB1775A431F47ED4D65E" ]]; then
    echo "Passed: heterozygote test."
else
    echo "Failed: heterozygote test."
    exit 1
fi

if [[ $( ./bin/veqtl-mapper --normal --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 | sha1sum | awk {'print toupper($1)'}) == "194CF372ACA1414B3927238EEA35324E9C40BE16" ]]; then
    echo "Passed: quantile normalisation test."
else
    echo "Failed: quantile normalisation test."
    exit 1
fi

echo "All tests completed successfully."
