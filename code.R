getwd()
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")

setwd("C:/Users/Tom Berry/Desktop/single_cell/")
SAMPLES <- c("510", "611", "993")
library(ArchR)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force = TRUE)
BiocManager::install("Matrix", force = TRUE)
addArchRGenome("hg38")

Arrowfiles <- createArrowFiles(
    inputFiles = c(paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[1], "_FC_fragments.tsv.gz"), 
                  paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[2], "_FC_fragments.tsv.gz"),
                  paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[3], "_FC_fragments.tsv.gz")),
    sampleNames = SAMPLES, 
    minTSS = 4, 
    minFrags = 1000)



project <- ArchRProject(
    ArrowFiles = Arrowfiles, outputDirectory = "singlecell", copyArrows = TRUE
)


saveArchRProject(project, outputDirectory = "save_proj_1", load = FALSE)

#Remove Doublets
project_doubletscore <- addDoubletScores(input = project, k = 10, knnMethod = "UMAP", LSIMethod = 1)
project_nodoublets <- filterDoublets(project_doubletscore)
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

#No batch correction
project_cluster_nobatchcorrection <- addClusters(
     input = project_LSI,
     reducedDims = "IterativeLSI",
     method = "Seurat",
     name = "Clusters",
 )
 
#Cluster counts no batch correction IterativeLSI
install.packages("pheatmap)
library(pheatmap)
cluster_counts <- as.data.frame(t(as.data.frame(as.vector((table(project_cluster_nobatchcorrection$Clusters))))))
rownames(cluster_counts) <- NULL
colnames(cluster_counts) <- names(table(project_cluster_nobatchcorrection$Clusters)) 

#Batch correction
project_cluster <- addHarmony(
     ArchRProj = project_LSI,
     reducedDims = "IterativeLSI",
     name = "Harmony",
     groupBy = "Sample"
 )
 
project_cluster <- addClusters(
     input = project_cluster,
     reducedDims = "Harmony",
     method = "Seurat",
     name = "Clusters_Harmony",
 )
 
cM_nobatchcorrection <- confusionMatrix(paste0(project_cluster_nobatchcorrection$Clusters), paste0(project_cluster_nobatchcorrection$Sample))
cM <- confusionMatrix(paste0(project_cluster$Clusters), paste0(project_cluster$Sample))
 
 #RNAseq data read in
 seurat.fc <- readRDS("seurat.pfc.final.rds")
seurat.fc$cellIDs <- gsub('FC-', '', seurat.fc$cellIDs)

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
 
cM_atac_rna <- as.matrix(confusionMatrix(atac_and_rna$Clusters, atac_and_rna$predictedGroup_Un))
preClust <- colnames(cM_atac_rna)[apply(cM_atac_rna, 1, which.max)]
cbind(preClust, rownames(cM_atac_rna))

unique(unique(atac_and_rna$predictedGroup_Un))
ExN <- paste0(c(1,3,4,7,10,12), collapse = "|")
InN <- paste0(c(2,9,13,14), collapse = "|")
RG <- paste0(c(6,5), collapse = "|")
MG <- paste0(8, collapse = "|")
