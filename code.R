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



project1 <- ArchRProject(
    ArrowFiles = Arrowfiles, outputDirectory = "singlecell", copyArrows = TRUE
)


saveArchRProject(project, outputDirectory = "save_proj_1", load = FALSE)

#Remove Doublets
doubletscore <- addDoubletScores(input = project1, k = 10, knnMethod = "UMAP", LSIMethod = 1)
project2 <- filterDoublets(doubletscore)
#Dimentionality Reduction with Iterative LSI
project2 <- addIterativeLSI(
    ArchRProj = project2,
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
project2 <- addClusters(
     input = project2,
     reducedDims = "IterativeLSI",
     method = "Seurat",
     name = "Clusters_no_batch_correction",
 )
 
#Cluster counts no batch correction IterativeLSI
install.packages("pheatmap")
library(pheatmap)
cluster_counts <- as.data.frame(t(as.data.frame(as.vector((table(project2$Clusters_no_batch_correction))))))
rownames(cluster_counts) <- NULL
colnames(cluster_counts) <- names(table(project2$Clusters_no_batch_correction)) 

#Batch correction
project2 <- addHarmony(
     ArchRProj = project2,
     reducedDims = "IterativeLSI",
     name = "Harmony",
     groupBy = "Sample"
 )
 
project2 <- addClusters(
     input = project2,
     reducedDims = "Harmony",
     method = "Seurat",
     name = "Clusters_Harmony",
 )

#Cluster counts batch corrected 
cluster_counts_harmony <- as.data.frame(t(as.data.frame(as.vector((table(project2$Clusters_Harmony))))))
rownames(cluster_counts) <- NULL
colnames(cluster_counts) <- names(table(project2$Clusters_Harmony)) 
 
cM_nobatchcorrection <- confusionMatrix(paste0(project2$Clusters_no_batch_correction), paste0(project2$Sample))
cM <- confusionMatrix(paste0(project2$Clusters_Harmony), paste0(project2$Sample))
 
#RNAseq data read in
RNAseq <- readRDS("seurat.pfc.final.rds")
RNAseq$cellIDs <- gsub('FC-', '', RNAseq$cellIDs)

project2 <- addGeneIntegrationMatrix(
     ArchRProj = project2,
     useMatrix = "GeneScoreMatrix",
     matrixName = "GeneIntegrationMatrix",
     reducedDims = "Harmony",
     seRNA = RNAseq,
     addToArrow = FALSE,
     groupRNA = "cellIDs",
     nameCell = "predictedCell_Un",
     nameGroup = "predictedGroup_Un",
     nameScore = "predictedScore_Un"
 )
 
cM <- as.matrix(confusionMatrix(project2$Clusters_Harmony, project2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1, which.max)]
cbind(preClust, rownames(cM))

unique(unique(project2$predictedGroup_Un))
ExN <- paste0(c(1,2,3,4,8,11), collapse = "|")
InN <- paste0(c(7,10,13,14), collapse = "|")
RG <- paste0(c(5,6), collapse = "|")
MG <- paste0(9, collapse = "|")
Else <- paste0(c(12,15), collapse="|")

clustExN <- rownames(cM)[grep(ExN, preClust)]
clustInN <- rownames(cM)[grep(InN, preClust)]
clustRG <- rownames(cM)[grep(RG, preClust)]
clustMG <- rownames(cM)[grep(MG, preClust)]
clustElse <- rownames(cM)[grep(Else, preClust)]

RNA <- RNAseq[grep(clust, RNAseq$cellIDs)]
