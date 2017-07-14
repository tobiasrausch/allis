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
    echo "FreeBayes variant calling"
    freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${HG} --genotype-qualities -b ${BAM} -v ${OP}.freebayes.vcf
    bgzip ${OP}.freebayes.vcf
    tabix ${OP}.freebayes.vcf.gz
else
    # Link variants
    bcftools view ${4} | bgzip > ${OP}.freebayes.vcf.gz
    tabix ${OP}.freebayes.vcf.gz
fi

# Phase against 1kGP
rm -f ${OP}.rename.fwd.chrs ${OP}.rename.rev.chrs
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    echo chr${CHR} ${CHR} >> ${OP}.rename.fwd.chrs
    echo ${CHR} chr${CHR} >> ${OP}.rename.rev.chrs
done
bcftools annotate -O b -o ${OP}.input.bcf --rename-chrs ${OP}.rename.fwd.chrs ${OP}.freebayes.vcf.gz
bcftools index ${OP}.input.bcf
FILES=""
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    echo "Eagle2 phasing chr${CHR}"
    eagle --numThreads ${THREADS} --vcfRef ${BASEDIR}/../refpanel/chr${CHR}.bcf --vcfTarget ${OP}.input.bcf --geneticMapFile ${BASEDIR}/../refpanel/genetic_map_hg19_withX.txt.gz --outPrefix ${OP}.chr${CHR}.eagle2 --vcfOutFormat b --chrom ${CHR} 2>&1 | gzip -c > ${OP}.chr${CHR}.eagle2.log.gz
    bcftools index ${OP}.chr${CHR}.eagle2.bcf
    FILES=${FILES}" "${OP}.chr${CHR}.eagle2.bcf
done
rm ${OP}.input.bcf ${OP}.input.bcf.csi
bcftools concat ${FILES} | fill-an-ac | bcftools annotate -O b -o ${OP}.eagle2join.bcf -x ^INFO/AC,INFO/AN,^FORMAT/GT -
bcftools index ${OP}.eagle2join.bcf
bcftools annotate -O b -o ${OP}.ealge2.bcf --rename-chrs ${OP}.rename.rev.chrs ${OP}.eagle2join.bcf
bcftools index ${OP}.ealge2.bcf
rm ${OP}.eagle2join.bcf ${OP}.eagle2join.bcf.csi
rm ${OP}.rename.fwd.chrs ${OP}.rename.rev.chrs
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    rm ${OP}.chr${CHR}.eagle2.bcf ${OP}.chr${CHR}.eagle2.bcf.csi
done

# Run Allis
echo "Running allis"
SAMPLE=`bcftools view ${OP}.ealge2.bcf | grep -m 1 "^#CHROM" | cut -f 10`
${BASEDIR}/allis -g ${HG} -v ${OP}.ealge2.bcf -p ${OP}.h1.bam -q ${OP}.h2.bam -s ${SAMPLE} -a ${OP}.tsv ${BAM} | gzip -c > ${OP}.allis.log.gz
