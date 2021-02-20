#!/bin/bash

# 1. for each population, subset to samples and coding variants
# 2. trim alt alleles
# 3. keep biallelic variable sites only within each population

POP=$1

echo "running for population ${POP}"

# 1.
for i in {1..22}; do echo $i; bcftools view ../../../../../Data/1000G/genotypes/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -R variants/coding_variants.chr${i}.txt -S samples/samples.${POP}.txt -v snps -Oz -o temp.${POP}.${i}.vcf.gz; bcftools index temp.${POP}.${i}.vcf.gz; done;

# 2.
echo "trimming alts"
for i in {1..22}; do echo $i; bcftools view temp.${POP}.${i}.vcf.gz --trim-alt-alleles -c 1 -Oz -o temp.${POP}.2.${i}.vcf.gz; bcftools index temp.${POP}.2.${i}.vcf.gz; done;

# 3.
echo "keeping biallelic variable sites"
for i in {1..22}; do echo $i; bcftools view temp.${POP}.2.${i}.vcf.gz -m2 -M2 -Oz -o ../vcfs/coding.${POP}.chr${i}.vcf.gz; done;

# 4.
echo "cleaning up temp.${POP} vcfs"
rm temp.${POP}.*.vcf.gz*

