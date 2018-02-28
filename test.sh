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

if [[ $( ./bin/veqtl-mapper --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 --eqtl data/eQTL | sha1sum | awk {'print toupper($1)'}) == "FB61D77C8BFAFEBEDAA7DE48921BCDDD50068901" ]]; then
    echo "Passed: eQTL test."
else
    echo "Failed: eQTL test."
    exit 1
fi

if [[ $( ./bin/veqtl-mapper --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 --eqtl data/eQTL --cov data/cov | sha1sum | awk {'print toupper($1)'}) == "5BE0E57D82CFA2B29319B12A89257D192D9150DD" ]]; then
    echo "Passed: single covariate test."
else
    echo "Failed: single covariate test."
    exit 1
fi

if [[ $( ./bin/veqtl-mapper --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 10000,4 --eqtl data/eQTL --cov data/cov_full | sha1sum | awk {'print toupper($1)'}) == "59C4BAF6199F5B015B13C710A63DC6AA9DA850CE" ]]; then
    echo "Passed: multiple covariates test."
else
    echo "Failed: multiple covariates test."
    exit 1
fi

echo "All tests completed successfully."
