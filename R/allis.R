library(ggplot2)
library(scales)
library(qvalue)
library(genefilter)

args=commandArgs(trailingOnly=TRUE)
as=read.table(args[1], header=T)

as=as[as$depth>20 & as$vaf>0.2 & as$vaf<0.8,]

m = nrow(as)
threshold = 0.1
# Bonferoni
as$bonferoni=(as$pvalue<threshold/m)
# p-adjust fdr
as$padj=(p.adjust(as$pvalue, method = 'fdr')<threshold)
# q-values
as$qval=(qvalue(as$pvalue)$qvalues<threshold)

print("Fraction allele-specific")
print(paste0("Bonferoni: ", mean(as$bonferoni)))
print(paste0("P-adjust FDR: ", mean(as$padj)))
print(paste0("q-values: ", mean(as$padj)))

# Plot pvalue distribution
png(paste0(args[1], ".pvalue.png"), height=600, width=600)
p1=ggplot(data=as, aes(x=pvalue))
p1=p1 + geom_histogram(aes(y=..density..), binwidth=0.05)
p1=p1 + xlab("pvalue") + ylab("Density")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("pvalue distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1
dev.off()
print(warnings())


# Plot VAF distribution
png(paste0(args[1], ".vaf.png"), height=600, width=600)
p1=ggplot(data=as, aes(x=vaf))
p1=p1 + geom_freqpoly(aes(y=..density..), binwidth=0.05)
p1=p1 + xlab("Variant Allele Frequency") + ylab("Density")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Variant Allele Frequency Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1
dev.off()
print(warnings())

