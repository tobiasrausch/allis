#!/bin/bash

export PATH=/g/funcgen/bin:${PATH}

if [ $# -ne 1 ]
then
    echo "**********************************************************************"
    echo "Allele-specific expression and open chromatin."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <hg19.fa>"
    echo ""
    exit -1
fi
    
if [ ! -f ../refpanel/chr1.bcf ]
then
    echo "Please download the reference panel first!"
    exit -1
fi

for SAMPLE in HG00096
do
    # Download RNA-Seq data
    wget "http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/${SAMPLE}.1.M_111124_6.bam"
    samtools index ${SAMPLE}.1.M_111124_6.bam

    # Fetch variant calls
    FILES=""
    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
    do
	COL=`bcftools view ../refpanel/chr${CHR}.bcf | grep -m 1 "^#CHROM"  | tr '\t' '\n' | awk '{print NR"\t"$0;}' | grep -w "${SAMPLE}" | cut -f 1`
	bcftools view ../refpanel/chr${CHR}.bcf | cut -f 1-9,${COL} | grep -v -P "\t<" | grep -v -P "\tAC=[0-9]*," | grep -v '0|0' | sed 's/\t1|1/\t1\/1/' | sed 's/\t0|1/\t0\/1/' | sed 's/\t1|0/\t0\/1/' | bcftools view -O b -o ${SAMPLE}.chr${CHR}.bcf -
	bcftools index ${SAMPLE}.chr${CHR}.bcf
	FILES=${FILES}" "${SAMPLE}.chr${CHR}.bcf
    done
    bcftools concat -O b -o ${SAMPLE}.variants.bcf ${FILES}
    bcftools index ${SAMPLE}.variants.bcf

    # Rename chrs
    rm -f rename.chrs
    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
    do
	echo ${CHR} chr${CHR} >> rename.chrs
	rm ${SAMPLE}.chr${CHR}.bcf ${SAMPLE}.chr${CHR}.bcf.csi
    done
    bcftools annotate -O b -o ${SAMPLE}.var.bcf --rename-chrs rename.chrs ${SAMPLE}.variants.bcf
    bcftools index ${SAMPLE}.var.bcf
    rm ${SAMPLE}.variants.bcf ${SAMPLE}.variants.bcf.csi rename.chrs

    # Run Allis
    ../src/allis.sh ${SAMPLE}.1.M_111124_6.bam ${HG} ${SAMPLE} ${SAMPLE}.var.bcf
done
