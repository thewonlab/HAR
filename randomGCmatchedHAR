####################################################################
## 02/04/2019
## Codes by Hyejung Won: hyejung_won@med.unc.edu
## Selecting random regions the GC-content of which is matched with HARs
####################################################################
library(GenomicRanges);
library(GenomicFeatures);
library(biomaRt);
library(BSgenome.Hsapiens.UCSC.hg19);
options(stringsAsFactors=FALSE);

##File containing HARs
bedfile = read.table("HAR_coordinate.bed")
colnames(bedfile)[1:3] = c("chrom","chromStart","chromEnd")
HAR = GRanges(bedfile$chrom, IRanges(bedfile$chromStart,bedfile$chromEnd))
HAR = sort(sortSeqlevels(HAR))

##Get GC content
peakseqs = getSeq(Hsapiens,seqnames(HAR),start(HAR),end(HAR));
GCcontent = rowSums(letterFrequency(peakseqs,c("G","C"),as.prob=TRUE));
##Add GC content to each promoter genomic range
mcols(HAR)$GCcontent = GCcontent;
seqlengths(HAR) = seqlengths(Hsapiens)[1:23];

##Loop over each promoter and choose a random interval with the same decile of GC-content
##Within the same chromosome
goodrandproms <- GRanges();
##Run the whole sampling 10 times to get enough random intervals for Hi-C
for (j in 1:10) {
  ##Now take many random samplings from the genome of the same size to determine how many have this many differentially open peaks
  for (i in 1:length(seqlevels(HAR))) {
    cat('chr',i,'\n');
    ##this chr promoter
    thischrprom = HAR[seqnames(HAR)==seqlevels(HAR)[i]];
    ##Flag to determine if all promoters within the chromosome have a GC matched random equivalent
    allpromsnotgoodflag = TRUE;
    ##Loop until all promoters within the chromosome have a GC matched random equivalent
    while (allpromsnotgoodflag) {
      cat('in redo loop\n');
      ##Get the number of ranges within this chromosome
      nranges = length(thischrprom);
      ##Get the maximum length of the beginning site on a chromosome
      maxseqlength = seqlengths(HAR)[i]-(max(width(thischrprom)));
      ##Get random interval starting point
      randstart = as.integer(runif(nranges, min=1,max=maxseqlength));
      ##Width of all promoters on this chromosome
      promwidths = width(thischrprom);
      ##Random promoters
      rand.gr = GRanges(seqnames=seqlevels(HAR)[i],IRanges(randstart,width=promwidths));
      ##Get GC content of these promoters
      randpromoterseqs = getSeq(Hsapiens, rand.gr);
      GCcontent = rowSums(letterFrequency(randpromoterseqs,c("G","C"),as.prob=TRUE));
      mcols(rand.gr)$GCcontent = GCcontent;
      ##Find where GC content is within 5% of the GC content of the original
      keepind = which(abs(mcols(rand.gr)$GCcontent - mcols(thischrprom)$GCcontent) <= 0.05);
      redoind = which(abs(mcols(rand.gr)$GCcontent - mcols(thischrprom)$GCcontent) > 0.05);
      ##The random promoters to keep
      if (length(keepind)>0)
        goodrandproms = c(goodrandproms,rand.gr[keepind]);
      ##If there are promoters which don't match GC, re-loop
      if (length(redoind)==0) {
        allpromsnotgoodflag = FALSE;
      } else {
        thischrprom = thischrprom[redoind];
        cat(length(redoind),'\n');
      }
    }
  }
}

goodrandpromsort = sort(goodrandproms)
olap = findOverlaps(goodrandpromsort,HAR)

goodrandpromsortrm = goodrandpromsort[c(-queryHits(olap))]
goodrandproms = goodrandpromsortrm

write.table(cbind(as.character(seqnames(goodrandproms)),start(goodrandproms),end(goodrandproms)),file="RandomGCMatchedHARs.bed",quote=FALSE,row.names=FALSE,col.names=FALSE);    
