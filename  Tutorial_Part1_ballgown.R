library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
pheno_data = read.csv("MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD.csv")
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
bg
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
save(bg, file='bg.rda')

results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

write.table(results_transcripts, "MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

write.table(results_transcripts, "MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

write.table(results_genes, "MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

write.table(sig_transcripts, "MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

write.table(sig_genes, "MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
quit(save="no")

#Bash scripts after 
head MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results.tsv

#How many genes are there on this chromosome?
grep -v feature MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results.tsv | wc -l
#60676

#How many passed filter in MDA-MB-231-control or MDA-MB-231-UBE3D-KD?

grep -v feature MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_filtered.tsv | wc -l
#9178

#How many differentially expressed genes were found on this chromosome (p-value < 0.05)?
grep -v feature MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_sig.tsv | wc -l
#2220

#Display the top 20 DE genes
grep -v feature MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_sig.tsv | sort -rnk 3 | head -n 20 
#Higher abundance in   MDA-MB-231-control

grep -v feature MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_sig.tsv | sort -nk 3 | head -n 20
#Higher abundance in  MDA-MB-231-UBE3D-KD


#Save all genes with P<0.05 to a new file.
grep -v feature MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > DE_genes.txt

head DE_genes.txt






