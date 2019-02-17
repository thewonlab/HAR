####################################################################
## 02/04/2019
## Codes by Hyejung Won: hyejung_won@med.unc.edu
## Fisher's exact test
####################################################################
options(stringsAsFactors = FALSE)
library(biomaRt)
hargene = unlist(read.table("HARgene.txt")) # evolutionary gene sets in hgnc symbol
targetgene = unlist(read.table("Targetgene.txt")) # genes that are associated with a certain disorder/cortical layer
getinfo = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") # Using Gencode v19 annotations
geneAnno1 = getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X"),mart=mart)
allgene= geneAnno1$hgnc_symbol
exomeLength = read.table("exomelength.txt", header=T) ## Downloaded from the Supplementary Table of Parikshak et al., 2013

metaMat = exomeLength[match(rownames(metaMat),exomeLength$hgnc),"exomelength"]
rownames(metaMat) = geneAnno1$hgnc_symbol
colnames(metaMat) = "exomeLength"

metaMat = cbind(metaMat, rep(NA, nrow(metaMat)))
listname = "HAR"
colnames(metaMat)[ncol(metaMat)] = listname
matchlist = match(rownames(metaMat), hargene)
metaMat[!is.na(matchlist), listname] = 1
metaMat[is.na(matchlist), listname] = 0 

metaMat = cbind(metaMat, rep(NA, nrow(metaMat)))
listname = "Disease"
colnames(metaMat)[ncol(metaMat)] = listname
matchlist = match(rownames(metaMat), targetgene)
metaMat[!is.na(matchlist), listname] = 1
metaMat[is.na(matchlist), listname] = 0

glm.out = glm(metaMat[,2]~metaMat[,3]+metaMat[,1],family=binomial) # If exome length is used as covariates
glm.out = glm(metaMat[,2]~metaMat[,3],family=binomial) # If there's no covariates: equivalent to Fisher's exact test
