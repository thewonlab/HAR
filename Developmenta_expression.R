####################################################################
## 02/04/2019
## Codes by Hyejung Won: hyejung_won@med.unc.edu
## Plotting developmental expression pattern
####################################################################
options(stringsAsFactors = FALSE)
library(gridExtra)
diseasedir =c("HARgene.txt","HGEgene.txt")
diseasename = c("HAR", "HGE") 
load("BrainSpan_expression.rda") # Brain gene expression data downloaded from BrainSpan: http://www.brainspan.org/; saved as datExpr 

datExpr = scale(datExpr,center = T, scale=FALSE) + 1
exprdat = vector(mode="list", length=length(diseasedir))
for(i in 1:length(diseasedir)){
  exprdat[[i]] = apply(datExpr[match(diseasehg[[i]], rownames(datExpr)),],2,mean,na.rm=T)
}

dat = c()
for(i in 1:length(diseasedir)){
  datframe = data.frame(Group=diseasename[[i]], Region=datMeta$Region, Age=datMeta$Age_unit, Unit=datMeta$Unit, Expr=exprdat[[i]])
  dat = rbind(dat, datframe)
}

dat$Group=factor(dat$Group, levels=c("HAR","Hugain"))

pdf(file="developmental_expression_HAR_Hugain.pdf")
ggplot(dat,aes(x=Age, y= Expr, fill=Group, color=Group)) + geom_jitter(alpha=.2,width=2)  + ylab("Normalized expression") + 
  geom_smooth(span=1) + facet_grid(.~Unit, scales="free_x",space="free_x") + xlab("Prenatal Age (week)         Postnatal Age (year)") +
  ggtitle("Brain Developmental Expression Trajectory") 
dev.off()
