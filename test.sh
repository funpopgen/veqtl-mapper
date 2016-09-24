#!/bin/bash
make
if [[ $( ./bin/VEQM --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype_veqm.vcf.gz --perm 10000,4 | sha1sum | awk {'print toupper($1)'}) == "166F762354D571A6E44DDC7C04318C5FE9B11866" ]]; then
    echo "Passed: variance test."
else
    echo "Failed: variance test."
    exit 1
fi

if [[ $( ./bin/VEQM --het --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype_veqm.vcf.gz --perm 10000,4 | sha1sum | awk {'print toupper($1)'}) == "5317CB6C7EDE41C8B06986FF039B5DBA369859F3" ]]; then
    echo "Passed: heterozygote test.."
else
    echo "Failed: heterozygote test.."
    exit 1
fi

echo "All tests completed successfully."
