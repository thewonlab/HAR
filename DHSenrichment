####################################################################
## 02/04/2019
## Codes by Hyejung Won: hyejung_won@med.unc.edu
## DHS enrichment analysis for HARs
####################################################################
options(stringsAsFactors=F)
library(GenomicRanges)
library(ggplot2)
load("phastCons.hg19.conservedregions.Rdata") # evolutionary conserved regions defined by phastCons score
load("DHS_celltype.rda") # coordinates for DHS in each cell type; saved as chromranges; downloaded from Roadmap Epigenomics: http://www.roadmapepigenomics.org/
bedfile = read.table("HAR_coordinate.bed")
colnames(bedfile)[1:3] = c("chrom","chromStart","chromEnd")

haranges = GRanges(seqnames=bedfile$Chr,ranges=IRanges(as.numeric(bedfile$Start),as.numeric(bedfile$End)))

olap = findOverlaps(haranges,chromranges);
harivalues = haranges[queryHits(olap)];
mcols(harivalues) = cbind(mcols(haranges[queryHits(olap)]), mcols(chromranges[subjectHits(olap)]))

olap = findOverlaps(conserved,chromranges);
conivalues = conserved[queryHits(olap)];
mcols(conivalues) = cbind(mcols(conserved[queryHits(olap)]), mcols(chromranges[subjectHits(olap)]))

chrommark = unique(chromranges$celltype)

pval = c()
enrichment = c()
harnum = length(unique(haranges))
connum = length(unique(conserved))
finaldat = data.frame(pval=rep(NA,length(chrommark)), or=rep(NA,length(chrommark)), ci1=rep(NA,length(chrommark)), ci2=rep(NA,length(chrommark)))
for(i in 1:length(chrommark)){
  harcell = unique(harivalues[harivalues$celltype==chrommark[i]])
  concell = unique(conivalues[conivalues$celltype==chrommark[i]])
  cont.mat = matrix(c(length(harcell), harnum-length(harcell), length(concell), connum-length(concell)), nrow=2) # contingency matrix
  fisherdat = fisher.test(cont.mat)
  
  finaldat[i,"pval"] = fisherdat$p.value
  finaldat[i,"or"] = fisherdat$estimate
  finaldat[i,"ci1"] = fisherdat$conf.int[1]
  finaldat[i,"ci2"] = fisherdat$conf.int[2]
  
  print(i)
}

pvalhisg = data.frame("EpigeneticMarks"=rownames(finaldat), "Pval"=-log10(finaldat$pval))
orhisg = data.frame("EpigeneticMarks"=rownames(finaldat), "OR"=finaldat$or, "ci1"=finaldat$ci1, "ci2"=finaldat$ci2)

pdf(file="HAR_enrichment_to_conserved_elements.pdf", width=15, height=4)
barpl = ggplot(pvalhisg, aes(x=EpigeneticMarks, y=Pval))
barpl + geom_bar(stat="identity",position="dodge", fill="midnightblue") + theme_bw() + labs(y="-log(P-val)", fill="", x="") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) 
barpl = ggplot(orhisg, aes(x=EpigeneticMarks, y=OR))
barpl + geom_bar(stat="identity", position="dodge", fill="midnightblue") + theme_bw() + labs(y="OR", fill="", x="") + 
  geom_errorbar(aes(ymin=ci1, ymax=ci2), width=.2, position=position_dodge(0.9)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) 
dev.off()

## phastCons score matched genomic regions
library(phastCons100way.UCSC.hg19)

phast <- phastCons100way.UCSC.hg19

phasthar = gscores(phast, haranges)
phasthar[phasthar$default==1]$default <- 0.999
phasthar$wid = width(phasthar)

phastrange = seq(0,0.9,0.1)

harbg = vector(length=10000, mode="list")
for(j in 1:10000){
  phastdefault = c()
  for(i in 1:length(phastrange)){
    harwirange = phasthar[phasthar$default>=phastrange[i] & phasthar$default<(phastrange[i]+0.1)]
    load(paste0("./phastCons/phastCons.",(i-1),".hg19.conservedregions.Rdata"))
    conserved = conserved[seqnames(conserved) %in% c(paste0("chr", c(1:22,"X")))]
    permconserved = conserved[sample(length(conserved), length(harwirange), replace=F)]
    phastdefault = c(phastdefault, permconserved)
  }
  phastdefault = unlist(GRangesList(phastdefault))
  phastdefault = unique(phastdefault)
  print(length(phastdefault))
  harbg[[j]] = phastdefault
  print(j)
}

pval = c()
enrichment = c()
fisherP = fisherOR = bgratio = data.frame(matrix(NA, nrow=10000, ncol=length(chrommark)))
colnames(fisherP) = colnames(fisherOR) = colnames(bgratio) = chrommark

olap = findOverlaps(haranges,chromranges);
dhsivalues = haranges[queryHits(olap)];
mcols(dhsivalues) = cbind(mcols(haranges[queryHits(olap)]), mcols(chromranges[subjectHits(olap)]))

haratio = c() 
for(i in 1:length(chrommark)){
  harchrom = unique(dhsivalues[dhsivalues$celltype==chrommark[i]])
  haratio = c(haratio, length(harchrom)/length(haranges))
}
names(haratio) = chrommark

for(j in 1:10000){
  hartarget = harbg[[j]]
  olap = findOverlaps(hartarget,chromranges);
  dhsibg = hartarget[queryHits(olap)];
  mcols(dhsibg) = cbind(mcols(hartarget[queryHits(olap)]), mcols(chromranges[subjectHits(olap)]))
  
  for(i in 1:length(chrommark)){
    harchrom = unique(dhsivalues[dhsivalues$celltype==chrommark[i]])
    bgchrom = unique(dhsibg[dhsibg$celltype==chrommark[i]])
    
    bgratio[j,i] = length(bgchrom)/length(hartarget)
    
    fisher2result = fisher.test(matrix(c(length(harchrom), length(haranges), length(bgchrom), length(hartarget)),2,2))
    fisherP[j,i] = fisher2result$p.value
    fisherOR[j,i] = fisher2result$estimate
  }
  
  print(j)
}

enrichmentP = c()
for(i in 1:length(chrommark)){
  bgratio.cell = bgratio[,i]
  enrichmentP = c(enrichmentP, sum(haratio[i] < bgratio.cell)/length(bgratio.cell))
}

orhisg = data.frame("EpigeneticMarks"=chrommark, "OR"=colMeans(fisherOR), "ci1"=colMeans(fisherOR)-apply(fisherOR, 2, sd), "ci2"=colMeans(fisherOR)+apply(fisherOR, 2, sd))
pdf(file="HAR_enrichment_to_conserved_matched_elements.pdf", width=15, height=4)
barpl = ggplot(orhisg, aes(x=EpigeneticMarks, y=OR))
barpl + geom_bar(stat="identity", position="dodge", fill="midnightblue") + theme_bw() + labs(y="OR", fill="", x="") + 
  geom_errorbar(aes(ymin=ci1, ymax=ci2), width=.2, position=position_dodge(0.9)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) 
dev.off()
