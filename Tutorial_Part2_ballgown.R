library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
outfile="/media/esra/new/workspace/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf"
pheno_data = read.csv("MDA-MB-231-control_vs_MDA-MB-231-UBE3D-KD.csv")

pheno_data
                       ids                type
1  MDA-MB-231-control-rep1  MDA-MB-231-control
2  MDA-MB-231-control-rep2  MDA-MB-231-control
3  MDA-MB-231-control-rep3  MDA-MB-231-control
4 MDA-MB-231-UBE3D-KD-rep1 MDA-MB-231-UBE3D-KD
5 MDA-MB-231-UBE3D-KD-rep2 MDA-MB-231-UBE3D-KD
6 MDA-MB-231-UBE3D-KD-rep3 MDA-MB-231-UBE3D-KD
                                                                                     path
1  /media/esra/new/workspace/rnaseq/expression/stringtie/ref_only/MDA-MB-231-control-rep1
2  /media/esra/new/workspace/rnaseq/expression/stringtie/ref_only/MDA-MB-231-control-rep2
3  /media/esra/new/workspace/rnaseq/expression/stringtie/ref_only/MDA-MB-231-control-rep3
4 /media/esra/new/workspace/rnaseq/expression/stringtie/ref_only/MDA-MB-231-UBE3D-KD-rep1
5 /media/esra/new/workspace/rnaseq/expression/stringtie/ref_only/MDA-MB-231-UBE3D-KD-rep2
6 /media/esra/new/workspace/rnaseq/expression/stringtie/ref_only/MDA-MB-231-UBE3D-KD-rep3

# Load the ballgown object from file
load('bg.rda')

# The load command, loads an R object from a file into memory in our R session. 
# You can use ls() to view the names of variables that have been loaded
ls()
 "bg"         "outfile"    "pheno_data"

# Print a summary of the ballgown object
bg
ballgown instance with 227818 transcripts and 6 samples

# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")

# View the last several rows of the FPKM table
tail(fpkm)
       FPKM.MDA-MB-231-control-rep1 FPKM.MDA-MB-231-control-rep2
227813                     0.000000                            0
227814                     0.000000                            0
227815                     0.000000                            0
227816                     0.000000                            0
227817                     0.000000                            0
227818                     0.135814                            0
       FPKM.MDA-MB-231-control-rep3 FPKM.MDA-MB-231-UBE3D-KD-rep1
227813                            0                             0
227814                            0                             0
227815                            0                             0
227816                            0                             0
227817                            0                             0
227818                            0                             0
       FPKM.MDA-MB-231-UBE3D-KD-rep2 FPKM.MDA-MB-231-UBE3D-KD-rep3
227813                             0                      0.000000
227814                             0                      0.000000
227815                             0                      0.000000
227816                             0                      0.000000
227817                             0                      0.000000
227818                             0                      0.131232

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)


# View the last several rows of the transformed FPKM table
tail(fpkm)
FPKM.MDA-MB-231-control-rep1 FPKM.MDA-MB-231-control-rep2
227813                    0.0000000                            0
227814                    0.0000000                            0
227815                    0.0000000                            0
227816                    0.0000000                            0
227817                    0.0000000                            0
227818                    0.1837266                            0
       FPKM.MDA-MB-231-control-rep3 FPKM.MDA-MB-231-UBE3D-KD-rep1
227813                            0                             0
227814                            0                             0
227815                            0                             0
227816                            0                             0
227817                            0                             0
227818                            0                             0
       FPKM.MDA-MB-231-UBE3D-KD-rep2 FPKM.MDA-MB-231-UBE3D-KD-rep3
227813                             0                     0.0000000
227814                             0                     0.0000000
227815                             0                     0.0000000
227816                             0                     0.0000000
227817                             0                     0.0000000
227818                             0                     0.1778948



# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(as.factor(pheno_data$type))+1,las=2,ylab='log2(FPKM+1)')

I found top gene name from biomart as
Gene stable ID	Gene name
ENSG00000149968	MMP3

for (val in 1: 227818)
{
    if (geneNames(bg)[val] == "MMP3") {
	print(val)
	}
}
It gives the results
[1] 38331
[1] 38332
[1] 38333
and I applied the rest of codes for the results. You can see outputs as : Tutorial_Part2_ballgown_output38333.pdf
Tutorial_Part2_ballgown_output38332.pdf
Tutorial_Part2_ballgown_output38331.pdf

# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[38333]
            38333 
"ENST00000524478" 

ballgown::geneNames(bg)[38333]
 38333 
"MMP3"

# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
boxplot(fpkm[38333,] ~ pheno_data$type, border=c(2,3), main=paste(ballgown::geneNames(bg)[38333],' : ', ballgown::transcriptNames(bg)[38333]),pch=19, xlab="Type", ylab='log2(FPKM+1)')

# Add the FPKM values for each sample onto the plot
points(fpkm[38333,] ~ jitter(c(2,2,2,1,1,1)), col=c(2,2,2,1,1,1)+1, pch=16)


# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[38333], bg, main=c('TST in all MDA-MB-231-UBE3D-KD samples'), sample=c('MDA-MB-231-UBE3D-KD-rep1', 'MDA-MB-231-UBE3D-KD-rep2', 'MDA-MB-231-UBE3D-KD-rep3'), labelTranscripts=TRUE)



plotTranscripts(ballgown::geneIDs(bg)[38333], bg, main=c('TST in all MDA-MB-231-control samples'), sample=c('MDA-MB-231-control-rep1', 'MDA-MB-231-control-rep2', 'MDA-MB-231-control-rep3'), labelTranscripts=TRUE)

# Close the PDF device where we have been saving our plots
dev.off()









