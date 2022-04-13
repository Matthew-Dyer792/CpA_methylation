# extract CpG from aligned, deduped bam file
#

CHH_LOCATION=$1
SAMPLEID=$2
OUTFOLDER=$3


# R script to do methylRaw data.
R --quiet --vanilla <<EO_RSCRIPT
library(methylKit)

files = (strsplit("${CHH_LOCATION}", ","))[[1]]
file.list = list(files[[1]], files[[2]])

ids = (strsplit("${SAMPLEID}", ","))[[1]]
sample_ids = list(ids[[1]], ids[[2]])

myobj=methRead(file.list, sample.id=sample_ids, assembly="hg38", treatment=c(0,1), context="CpG", mincov=5)

setwd("${OUTFOLDER}")

project_name = paste0(sample_ids[[1]], "_vs_", sample_ids[[2]])

for (var in 1:lengths(sample_ids)) {
  pdf(paste0(sample_ids[[var]], "_methylation.pdf"))
  getMethylationStats(myobj[[var]],plot=TRUE,both.strands=FALSE)
  dev.off()
  pdf(paste0(sample_ids[[var]], "_CpG_coverage.pdf"))
  getCoverageStats(myobj[[var]],plot=TRUE,both.strands=FALSE)
  dev.off()
}

meth=unite(myobj, destrand=FALSE, mc.cores=8)

myDiff=calculateDiffMeth(meth, mc.cores=8)

myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

pdf(paste(project_name, "chromosome_methylation.pdf", sep="_"))
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()

pdf(paste0(project_name, "_volcano_plot.pdf"))
plot(myDiff25p\$meth.diff,-log(myDiff25p\$qvalue))
abline(v=-25, col = 'red')
abline(v=25, col = 'red')
abline(h=-log(0.01), col = 'blue')
title(project_name)
dev.off()

bed <- getData(myDiff25p)
bed <- bed[c(1,2,3,4)]
write.table(bed, file = paste0(project_name, ".my_diff.bed"), sep = "\t", quote = F)

bedHypo <- getData(myDiff25p.hypo)
bedHypo <- bedHypo[c(1,2,3,4)]
write.table(bedHypo, file = paste0(project_name, ".my_diff.hypo.bed"), sep = "\t", quote = F)

bedHyper <- getData(myDiff25p.hyper)
bedHyper <- bedHyper[c(1,2,3,4)]
write.table(bedHyper, file = paste0(project_name, ".my_diff.hyper.bed"), sep = "\t", quote = F)

q();
EO_RSCRIPT

