#!/bin/bash

if [ \( $# -lt 3 \) -o \( $# -gt 4 \) ]
then
    echo "**********************************************************************"
    echo "Allele-specific expression and open chromatin."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <aligned.bam> <reference.fasta> <outprefix> [<input.vcf.gz>]"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=/g/funcgen/bin/:${PATH}

# CMD params
THREADS=4
BAM=${1}
HG=${2}
OP=${3}

# Input Variants
if [ $# -eq 3 ]
then
    # Call variants
    freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${HG} --genotype-qualities -b ${BAM} -v ${OP}.freebayes.vcf
    bgzip ${OP}.freebayes.vcf
    tabix ${OP}.freebayes.vcf.gz
else
    # Link variants
    bcftools view ${4} | bgzip > ${OP}.freebayes.vcf.gz
    tabix ${OP}.freebayes.vcf.gz
fi

# Phase against 1kGP
FILES=""
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    echo "Eagle2 phasing chr${CHR}"
    bcftools view ${OP}.freebayes.vcf.gz chr${CHR} | sed "s/^##contig=<ID=chr${CHR},/##contig=<ID=${CHR},/" | grep -v "^##contig=<ID=chr" | sed "s/^chr${CHR}\t/${CHR}\t/" | bcftools view -O b -o ${OP}.input.chr${CHR}.bcf -
    bcftools index ${OP}.input.chr${CHR}.bcf
    eagle --numThreads ${THREADS} --vcfRef ${BASEDIR}/../refpanel/chr${CHR}.bcf --vcfTarget ${OP}.input.chr${CHR}.bcf --geneticMapFile ${BASEDIR}/../refpanel/genetic_map_hg19_withX.txt.gz --outPrefix ${OP}.chr${CHR}.eagle2 --vcfOutFormat b --chrom ${CHR} 2>&1 | gzip -c > ${OP}.chr${CHR}.eagle2.log.gz
    bcftools index ${OP}.chr${CHR}.eagle2.bcf
    rm ${OP}.input.chr${CHR}.bcf ${OP}.input.chr${CHR}.bcf.csi
    bcftools view ${OP}.chr${CHR}.eagle2.bcf | fill-an-ac | bcftools annotate -x ^INFO/AC,INFO/AN,^FORMAT/GT - | sed "s/^##contig=<ID=${CHR},/##contig=<ID=chr${CHR},/" | sed "s/^${CHR}\t/chr${CHR}\t/" | bcftools view -O b -o ${OP}.output.chr${CHR}.bcf -
    bcftools index ${OP}.output.chr${CHR}.bcf
    rm ${OP}.chr${CHR}.eagle2.bcf ${OP}.chr${CHR}.eagle2.bcf.csi
    FILES=${FILES}" "${OP}.output.chr${CHR}.bcf
done
bcftools concat ${FILES} | grep -v "^##contig=<ID=" | bgzip > ${OP}.eagle2.vcf.gz
tabix ${OP}.eagle2.vcf.gz
bcftools view -O b -o ${OP}.eagle2.bcf ${OP}.eagle2.vcf.gz
bcftools index ${OP}.eagle2.bcf
rm ${OP}.eagle2.vcf.gz ${OP}.eagle2.vcf.gz.tbi
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    rm ${OP}.output.chr${CHR}.bcf ${OP}.output.chr${CHR}.bcf.csi
done

