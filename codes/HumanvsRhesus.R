####################################################################
## 02/04/2019
## Codes by Hyejung Won: hyejung_won@med.unc.edu
## Human vs. Rhesus expression data
## Files and codes downloaded from: git clone https://github.com/AllenBrainAtlas/DevRhesusLMD.git. 
####################################################################
# Load libraries 
library(reshape2)
library(segmented)
library(ggplot2)
library(scatterplot3d)
library(RColorBrewer)
library(limma)

# Load functions 
source(file="../src/fReorderFactorLevels.R")
source(file="../src/fConvertPcdtoEventScore.R")

# Try loading cached fits
dev.pred.fn <- "../cache/dev_expr_species/dev.expr_dev.pred_2hs.RData" 
try(load(dev.pred.fn), silent=TRUE)

# Get genes you would like to test for
evolgene = unlist(read.table("evol_gene.txt"))

# Calc species expression variation across development 
rh.var <- apply(dev.pred[["macaque"]], 1, function(x) sd(x) / mean(x)) 
h.var <- apply(dev.pred[["human"]], 1, function(x) sd(x) / mean(x))
rat.var <- apply(dev.pred[["rat"]], 1, function(x) sd(x) / mean(x))

# Let's check whether the variance differ for the genes that you selected for
checkvar = function(genelist){
  h.var.test = h.var[names(h.var) %in% genelist]
  h.var.ntest = h.var[!(names(h.var) %in% genelist)]
  
  rh.var.test = rh.var[names(rh.var) %in% genelist]
  rh.var.ntest = rh.var[!(names(rh.var) %in% genelist)]
  
  rat.var.test = rat.var[names(rat.var) %in% genelist]
  rat.var.ntest = rat.var[!(names(rat.var) %in% genelist)]
  
  wch = wilcox.test(h.var.test, h.var.ntest)
  wcrh = wilcox.test(rh.var.test, rh.var.ntest)
  wcrat = wilcox.test(rat.var.test, rat.var.ntest)
  
  th = t.test(h.var.test, h.var.ntest)
  trh = t.test(rh.var.test, rh.var.ntest)
  trat = t.test(rat.var.test, rat.var.ntest)
  
  hrhvar.test = log(h.var.test/rh.var.test) 
  hrhvar.ntest = log(h.var.ntest/rh.var.ntest) 
  wchrh = wilcox.test(hrhvar.test, hrhvar.ntest)
  thrh = t.test(hrhvar.test, hrhvar.ntest)
  
  vartestlist = list(wch, wcrh, wcrat, wchrh, th, trh, trat, thrh)
  names(vartestlist) = c("var_human:wilcox", "var_rhesus:wilcox", "var_rat:wilcox", "var_human_vs_rhesus:wilcox", 
                         "var_human:t-test", "var_rhesus:t-test", "var_rat:t-test", "var_human_vs_rhesus:t-test")
  return(vartestlist)
  
}

evol.var = checkvar(evolgene) 

# Let's plot developmental expression trajectories for a gene list 
dev.Expr = dev.expr[dev.expr$species %in% c("human", "macaque"),]
dev.Expr$evol = as.factor(ifelse(dev.Expr$gene %in% evolgene, 1, 2))

## Let's select out the Z-scores in critical epochs
dev.expr.bound = dev.Expr
dev.expr.bound$estage = ifelse(dev.expr.bound$escore<0.5, 0, 1)
dev.expr.bound$estage = dev.expr.bound$estage + ifelse(dev.expr.bound$escore>0.5, 1, 0)
dev.expr.bound$estage = dev.expr.bound$estage + ifelse(dev.expr.bound$escore>1, 1, 0)

g2 <- ggplot(dev.expr.bound, aes(x=evol, y=exprz)) + geom_boxplot(aes(fill=species)) + facet_wrap( ~ estage) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = pal1) +
  scale_fill_manual(values = pal1) +
  labs( y="Normalized expression Z-score", 
        x="",
        title = "Fetal Brain Developmental Expression: Evol vs. non-Evol")

# Calculate statistics
dev.expr.bound1 = dev.Expr[dev.Expr$escore>0.5 & dev.Expr$escore<1, ]
dev.expr.bound1$index = paste(dev.expr.bound1$species, dev.expr.bound1$evol, sep=":")
dev.expr.bound2 = dev.Expr[dev.Expr$escore<0.5, ]
dev.expr.bound2$index = paste(dev.expr.bound2$species, dev.expr.bound2$evol, sep=":")

aov1 = aov(exprz~index, data=dev.expr.bound1)
TukeyHSD(aov1)
aov2 = aov(exprz~index, data=dev.expr.bound2)
TukeyHSD(aov2)

g3 <- ggplot(dev.expr.bound1, aes(exprz)) + geom_density(aes(fill=factor(index)), alpha=0.5) + 
  scale_fill_manual(values = pal1) +
  labs( y="Density", 
        x="Normalized expression Z-score", 
        title="Evol vs. non-Evol")

# Let's calculate Delta(Z) based on escore
dev.expr.mq = dev.expr[dev.expr$species=="macaque",]
dev.expr.hu = dev.expr[dev.expr$species=="human",]
timepoint.mq = c(0.46,0.51,0.77)
timepoint.hu = c(0.48,0.54,0.76) 
dev.expr.mq = dev.expr.mq[dev.expr.mq$escore %in% timepoint.mq, ]
dev.expr.hu = dev.expr.hu[dev.expr.hu$escore %in% timepoint.hu, ]
dev.mq = cbind(dev.expr.mq, dev.expr.hu)
sum(dev.mq[,2] != dev.mq[,8])
dev.mq$zdiff = dev.mq[,11] - dev.mq[,5]
dev.mq = dev.mq[,c(2,12,13)]
dev.mq$evol = ifelse(dev.mq$gene %in% evolgene, "evol", "nonevol")

g4 <- ggplot(dev.mq, aes(zdiff)) + geom_density(aes(fill=factor(evol)), alpha=0.5) + 
  scale_fill_manual(values = pal1) +
  labs( y="Density", 
        x=expression(Delta*" Normalized expression Z-score (Human-Macaque)"), 
        title="Evol vs. non-Evol")

pdf("Evol_expression_traits.pdf", height=5, width=7)
plot(g1)
plot(g2)
plot(g3)
plot(g4)
dev.off()

# Load functions 
source(file="../src/fReorderFactorLevels.R")
source(file="../src/fConvertPcdtoEventScore.R")

# Load genes of interest
load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
hargene = unlist(read.table("HARgene.txt"))
hugaingene = unlist(read.table("Hugaingene.txt"))
hlegene = unlist(read.table("HLEgene.txt"))
hgegene = unlist(read.table("HGEgene.txt"))

# Load breakpoints for increasing genes in all species/regions
bp.df <- read.csv(file = "../cache/dev_expr_species/species_region_bp.csv", 
                  stringsAsFactors = FALSE)  # SI Table 11 subset
bp.df$bpmin <- mapply(ConvertPcdtoEventScore, log2(bp.df$bp..pcd. - 1.96 * bp.df$bpse), 
                      bp.df$species)
bp.df$bpmax <- mapply(ConvertPcdtoEventScore, log2(bp.df$bp..pcd. + 1.96 * bp.df$bpse), 
                      bp.df$species)
bp.subset <- subset(bp.df, species %in% c("macaque", "human")) # & region == "ACG")
bp.escore <- dcast(bp.subset, region + gene ~ species, value.var = "bp..event.score.")
bp.min <- dcast(bp.subset, region + gene ~ species, value.var = "bpmin")
bp.max <- dcast(bp.subset, region + gene ~ species, value.var = "bpmax")
bp.diff <- data.frame(bp.escore, bp.min[, 3:4], bp.max[, 3:4],
                      h_m_escore = bp.escore$human - bp.escore$macaque,
                      hmin_mmax = bp.min$human - bp.max$macaque,
                      hmax_mmin = bp.max$human - bp.min$macaque)

bp.diff$sig <- apply(bp.diff, 1, function(x) (as.numeric(x["h_m_escore"]) > 0 & 
                                                min(as.numeric(x[c("hmin_mmax", "hmax_mmin")])) > 0) | 
                       (as.numeric(x["h_m_escore"]) < 0 & 
                          max(as.numeric(x[c("hmin_mmax", "hmax_mmin")])) < 0))

### This file now contains breakpoints for a given gene in different speceis
# Compare breakpoints based on used-defined groups
bp.diff$har = ifelse(bp.diff$gene %in% hargene, "HAR", "nHAR")
bp.diff$hugain = ifelse(bp.diff$gene %in% hugaingene, "HGE:FB", "nHGE")
bp.diff$hge = ifelse(bp.diff$gene %in% hgegene, "HGE:AB", "nHGE")
bp.diff$hle = ifelse(bp.diff$gene %in% hlegene, "HLE:AB", "nHLE")

wilcox.test(h_m_escore~har, data=bp.diff) 
wilcox.test(h_m_escore~hugain, data=bp.diff) 
wilcox.test(h_m_escore~hge, data=bp.diff)
wilcox.test(h_m_escore~hle, data=bp.diff)

# Summarize genes that have significantly different breakpoints between human and macaque based on event scores
diff.summary <- data.frame()

# Human early genes
for (region1 in unique(bp.subset$region)) {
  genes <- subset(bp.diff, region == region1 & sig == TRUE & h_m_escore < 0)$gene
  diff1 <- with(subset(bp.subset, region == region1 & species == "human" & gene %in% genes),
                data.frame(set = "human_early", region1,  
                           Plateau = sum(! slope2sig), 
                           Decreasing = sum(slope2sig & slope2 < 0), 
                           Increasing = sum(slope2sig & slope2 > 0)))
  diff.summary <- rbind(diff.summary, diff1)
}

humanearlygene = diffdec = diffinc = c()
for(region1 in c("V1","ACG")){
  genes <- subset(bp.diff, region==region1 & sig==TRUE & h_m_escore<0)$gene
  humanearlygene = union(humanearlygene, genes)
  diffdec <- union(diffdec, subset(bp.subset, region == region1 & species == "human" & gene %in% genes & slope2sig & slope2<0)$gene)
  diffinc <- union(diffinc, subset(bp.subset, region == region1 & species == "human" & gene %in% genes & slope2sig & slope2>0)$gene)
}

# Run Fisher's exact test
bgenes = unique(bp.diff$gene) # These genes show distinct changes in expression trajectories, or breakpoints. A total of 179 increasing and 179 decreasing genes met these criteria in all three species

harbg = intersect(hargene, bgenes)
harhuearly = intersect(humanearlygene, hargene)
fisher.test(matrix(c(length(harhuearly), length(setdiff(harbg, harhuearly)), 
                     length(setdiff(humanearlygene, harhuearly)), length(setdiff(bgenes, union(harbg, humanearlygene)))),2,2))

hugainbg = intersect(hugaingene, bgenes)
hugainhuearly = intersect(humanearlygene, hugaingene)
fisher.test(matrix(c(length(hugainhuearly), length(setdiff(hugainbg, hugainhuearly)), 
                     length(setdiff(humanearlygene, hugainhuearly)), length(setdiff(bgenes, union(hugainbg, humanearlygene)))),2,2))

# Human late genes
for (region1 in unique(bp.subset$region)) {
  genes <- subset(bp.diff, region == region1 & sig == TRUE & h_m_escore > 0)$gene
  diff1 <- with(subset(bp.subset, region == region1 & species == "human" & gene %in% genes),
                data.frame(set = "human_late", region1,  
                           Plateau = sum(! slope2sig), 
                           Decreasing = sum(slope2sig & slope2 < 0), 
                           Increasing = sum(slope2sig & slope2 > 0)))
  diff.summary <- rbind(diff.summary, diff1)
}

diff.summaryl <- melt(diff.summary, id = c("set", "region1"))
colnames(diff.summaryl)[3:4] <- c("slope_after_breakpoint", "num_genes")
num.genes <- diff.summaryl[diff.summaryl$set == "human_early", "num_genes"] 
diff.summaryl[diff.summaryl$set == "human_early", "num_genes"] <- -num.genes

diff.summaryl$region1 <- factor(diff.summaryl$region1, 
                                levels = c("STR", "AM", "HP", "V1", "ACG"))
