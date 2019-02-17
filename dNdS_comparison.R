####################################################################
## 02/04/2019
## Codes by Hyejung Won: hyejung_won@med.unc.edu
## dN/dS comparison
####################################################################

options(stringsAsFactors=F)
FBfile = "HARgenes.txt" # Hi-C interacting genes for HARs
mousednds = "mouse_dnds.rda"

library(biomaRt)

FP = unlist(read.table(FBfile))

## Protein-coding genes from Gencode v19
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") # Using Gencode v19 annotations
geneAnno1 <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X"),mart=mart)
geneAnno1 <- geneAnno1[geneAnno1[,"gene_biotype"]=="protein_coding",]; geneAnno1 <- geneAnno1[!duplicated(geneAnno1[,"hgnc_symbol"]) & geneAnno1[,"hgnc_symbol"]!="",] ## 19163 with hgnc symbols after removing duplicates. We use gene symbols here because the RDNV data comes in gene symbol format. Otherwise, I prefer to use ENSG IDs.

## Select out protein-coding genes 
FBpc = FB[FB %in% geneAnno1$ensembl_gene_id]
FBpc = geneAnno1[geneAnno1$ensembl_gene_id %in% FB, ]
FBensg = FBpc$ensembl_gene_id
FBhgnc = FBpc$hgnc_symbol

mart.hs = mart
listmarths = listAttributes(mart.hs)
hommart = listmarths[grep("homolog", listmarths$name),1] # homologs mart

allpc = geneAnno1[!(geneAnno1$ensembl_gene_id %in% FB), ]
allensg = allpc$ensembl_gene_id
allhgnc = allpc$hgnc_symbol

FBdnds_mouse = c()
FBdiffdnds = c()
for(i in 1:length(FBensg)){
  dn = getBM(attributes = "mmusculus_homolog_dn", filters = "ensembl_gene_id", values = FBensg[i], mart = mart.hs)
  ds = getBM(attributes = "mmusculus_homolog_ds", filters = "ensembl_gene_id", values = FBensg[i], mart = mart.hs)
  if(dim(dn)[1]!=0 & dim(ds)[1]!=0){
    if(dim(dn)[1]==dim(ds)[1]){
      FBdnds_mouse = rbind(FBdnds_mouse, c(FBensg[i], dn/ds))
      print(i)
    }else{
      FBdiffdnds = c(FBdiffdnds, i)
    }
  }
}

FBDdnds_mouse = c()
FBDdiffdnds = c()
for(i in 1:length(FBDensg)){
  dn = getBM(attributes = "mmusculus_homolog_dn", filters = "ensembl_gene_id", values = FBDensg[i], mart = mart.hs)
  ds = getBM(attributes = "mmusculus_homolog_ds", filters = "ensembl_gene_id", values = FBDensg[i], mart = mart.hs)
  if(dim(dn)[1]!=0 & dim(ds)[1]!=0){
    if(dim(dn)[1]==dim(ds)[1]){
      FBDdnds_mouse = rbind(FBDdnds_mouse, c(FBDensg[i], dn/ds))
      print(i)
    }else{
      FBDdiffdnds = c(FBDdiffdnds, i)
    }
  }
}

alldnds_mouse = c()
alldiffdnds = c()
for(i in 1:length(allensg)){
  dn = getBM(attributes = "mmusculus_homolog_dn", filters = "ensembl_gene_id", values = allensg[i], mart = mart.hs)
  ds = getBM(attributes = "mmusculus_homolog_ds", filters = "ensembl_gene_id", values = allensg[i], mart = mart.hs)
  if(dim(dn)[1]!=0 & dim(ds)[1]!=0){
    if(dim(dn)[1]==dim(ds)[1]){
      alldnds_mouse = rbind(alldnds_mouse, c(allensg[i], dn/ds))
      print(i)
    }else{
      alldiffdnds = c(alldiffdnds, i)
    }
  }
}

FBdnds = as.numeric(unlist(FBdnds_mouse[,2]))
FBdnds = FBdnds[!is.na(FBdnds)]; FBdnds = FBdnds[FBdnds!=0 & FBdnds!="Inf"]

Adnds = as.numeric(unlist(alldnds_mouse[,2]))
Adnds = Adnds[!is.na(Adnds)]; Adnds = Adnds[Adnds!=0 & Adnds!="Inf"]

logFBdnds = log2(FBdnds)
logAdnds = log2(Adnds)

ks.test(logFBdnds, logAdnds) 

densFB = density(log2(Pdnds)) 
densA = density(log2(Adnds))

xlim = range(-10,2); ylim = range(0,densFB$y, densA$y)
Hcol = adjustcolor("darkorange1", 0.6)
Lcol = adjustcolor("dodgerblue1", 0.6)

pdf(file="density_of_dnds_mouse.pdf", width=6, height=2.9)
par(mfrow=c(1,2),mar=c(3,3,1,1))
plot(densFB, xlim=xlim, ylim=ylim, main ='Distribution', xlab="", ylab="",cex.main=1.1, cex.axis=1.0, xaxt="n")
polygon(densFB, density = -1, col = Hcol)
polygon(densA, density = -1, col = Lcol)
axis(1, at=seq(-10,2,by=2), las=1)
mtext(side=1, text ='log2(dN/dS)', line=2, cex=1.1)
mtext(side=2, text ='Density', line=2, cex=1.1)
legend('topleft',c('HAR','All'),
       fill = c(Hcol, Lcol), bty = 'n', border = NA, cex=1.0)
dev.off()
