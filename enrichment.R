

library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
filename<- "c7.all.v7.1.entrez.gmt"
gmtfile <- system.file(filename)
c6 <- read.gmt(gmtfile)



yourEntrezIdList<- c(8924,8924,4522,84193,9697,80829,6139,23607,4288,644815,153339,55165,6009,308,5049,54476,9392,10553,7402,8924)


ImmunSigEnrich <-clusterProfiler:: enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.5)
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenes.csv")


goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Hs.eg.db, ont="ALL",pAdjustMethod="BH",pvalueCutoff = 0.5,readable= TRUE)

keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "hsa",pAdjustMethod="BH",pvalueCutoff = 0.5)
write.csv(keggEnrich,"MyKEGGRelatedGenes.csv")



dotplot(goEnrich, showCategory=30)
cnetplot(goEnrich, foldChange=geneList)
plotGOgraph(goEnrich)


mutate(goEnrich, qscore = -log(p.adjust, base=10)) %>% barplot(x="qscore")

threshold <- res_table$E$padj < padj.cutoff & abs(res_tableOE$log2FoldChange) > lfc.cutoff

