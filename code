#Check current working directory
getwd()

#Get required packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force = TRUE)
BiocManager::install("Matrix", force = TRUE)

#Setting Rtools on PATH
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")

#Set working directory
setwd("C:/Users/Tom/Desktop/single_cell/")
#Set sample name vector and load ArchR package
SAMPLES <- c("510", "611", "993")
library(ArchR)

#Set genome reference to hg38, this is required on each start-up
addArchRGenome("hg38")

#Create Arrowfiles for compatible use of the data with ArchR
Arrowfiles <- createArrowFiles(
    inputFiles = c(paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[1], "_FC_fragments.tsv.gz"), 
                  paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[2], "_FC_fragments.tsv.gz"),
                  paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[3], "_FC_fragments.tsv.gz")),
    sampleNames = SAMPLES, 
    minTSS = 4, 
    minFrags = 1000)

#Doubletscores
doubleScores <- addDoubletScores(
    input = Arrowfiles, k = 10, knnMethod = "UMAP", LSIMethod = 1
	)

#Make ArchRProject file, "copyArrows = TRUE" makes a backup copy of the arrowfiles.
project <- ArchRProject(
    ArrowFiles = Arrowfiles, outputDirectory = "singlecell", copyArrows = TRUE
)
