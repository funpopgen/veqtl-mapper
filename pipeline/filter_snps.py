#! /usr/bin/env python3
from __future__ import print_function
import argparse
import gzip
import os.path
import subprocess
import sys


def check_args(args):
    """Checks that commandline options are valid:
Miminum count in genotype groups is greater than zero.
Input file exists.
Output file, compressed output file and index for output file do not
exist and so won't be overwritten.
"""
    if args.threshold < 0:
        print("Minimum count must be greater than zero.", file=sys.stderr)
        sys.exit()
    if not os.path.isfile(args.input):
        print("Input file,", args.input, "doesn't exist.", file=sys.stderr)
        sys.exit()
    if os.path.isfile(args.output):
        print(args.output, "would be overwritten.", file=sys.stderr)
        sys.exit()
    if os.path.isfile(args.output + ".gz"):
        print(args.output + ".gz", "would be overwritten.", file=sys.stderr)
        sys.exit()
    if os.path.isfile(args.output + ".gz.csi"):
        print(args.output + ".gz.csi", "would be overwritten.",
              file=sys.stderr)
        sys.exit()    
    
def filter_snps(args):
    """Given an input vcf file print only snps where all 3 genotype groups have
at least a certain number of individuals"""
    threshold = args.threshold - 1
    with gzip.open(args.input, 'rt') as f:
        with open(args.output, 'w') as g:
            for line in f:
                if line[0] == "#":
                    _ = g.write(line.strip() + '\n')
                else:
                    line = line.strip().split()
                    count_hom1, count_het, count_hom2 = 0, 0, 0
                    for snp in line[9:]:
                        if snp == "0|0":
                            count_hom1 += 1
                        elif snp == "1|1":
                            count_hom2 += 1
                        else:
                            count_het += 1
                    if count_hom1 > threshold and count_hom2 > threshold and count_het > threshold:
                        _ = g.write('\t'.join(line) + '\n')
    subprocess.call(['bgzip', args.output])
    subprocess.call(['bcftools', 'index', args.output + ".gz"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filters a vcf file for SNPs with at least a certain count in all three genotype groups. Warning: only works on vcf files where genotype is specified only by the phased GT flag (such as the standard 1000 genomes downloads).")
    parser.add_argument("--input",
                        default="Genotypes.individuals_filtered.maf0.05.vcf.gz",
                        help="Name of vcf file to be filtered.")
    parser.add_argument("--output",
                        default="Genotypes.individuals_filtered.parent.of.origin.vcf",
                        help="Name of output vcf file (without .gz suffix).")
    parser.add_argument("--threshold", default=50, type=int,
                        help="Minimum count in genotype groups.")
    args = parser.parse_args()
    check_args(args)
    filter_snps(args)

        
    
