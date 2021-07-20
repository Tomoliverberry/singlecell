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
saveArchRProject(project, outputDirectory = "save_proj_1", load = FALSE)

#Remove Doublets
project_nodoublets <- filterDoublets(project)
#Dimentionality Reduction with Iterative LSI
project_LSI <- addIterativeLSI(
    ArchRProj = project_nodoublets,
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

BiocManager::install("harmony", force = TRUE)

#Clustering without batch correction
project_cluster_nobatchcorrection <- addClusters(
     input = project_LSI,
     reducedDims = "IterativeLSI",
     method = "Seurat",
     name = "Clusters",
 )

#Batch correction
project_batch_corrected <- addHarmony(
     ArchRProj = project_LSI,
     reducedDims = "IterativeLSI",
     name = "Harmony",
     groupBy = "Sample"
 )
 #Clustering batch corrected project file
 project_cluster <- addClusters(
     input = project_batch_corrected,
     reducedDims = "Harmony",
     method = "Seurat",
     name = "Clusters",
 )

#Make 2 matricies, showing number of cells per cluster for each donor. Compare to view reduced number of cells from batch correction.
cM_nobatchcorrection <- confusionMatrix(paste0(project_cluster_nobatchcorrection$Clusters), paste0(project_cluster_nobatchcorrection$Sample))
cM <- confusionMatrix(paste0(project_cluster$Clusters), paste0(project_cluster$Sample))
 
#RNAseq data read in
seurat.fc <- readRDS("seurat.pfc.final.rds")
seurat.fc$cellIDs <- gsub('FC-', '', seurat.fc$cellIDs)

#Integrate atac-seq data and rna-seq data together using batch corrected data going forward
atac_and_rna <- addGeneIntegrationMatrix(
     ArchRProj = project_cluster,
     useMatrix = "GeneScoreMatrix",
     matrixName = "GeneIntegrationMatrix",
     reducedDims = "Harmony",
     seRNA = seurat.fc,
     addToArrow = FALSE,
     groupRNA = "cellIDs",
     nameCell = "predictedCell_Un",
     nameGroup = "predictedGroup_Un",
     nameScore = "predictedScore_Un"
 )
