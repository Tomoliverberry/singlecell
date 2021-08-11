getwd()
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")

setwd("C:/Users/Tom/Desktop/single_cell/")

library(ArchR)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force = TRUE)
BiocManager::install("Matrix", force = TRUE)
addArchRGenome("hg38")

SAMPLES <- c("510", "611", "993")

Arrowfiles <- createArrowFiles(
    inputFiles = c(paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[1], "_FC_fragments.tsv.gz"), 
                  paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[2], "_FC_fragments.tsv.gz"),
                  paste0("C:/Users/Tom/Desktop/single_cell/FastFile-3aEnrFqAw78jBfpr/", SAMPLES[3], "_FC_fragments.tsv.gz")),
    sampleNames = SAMPLES, 
    minTSS = 4, 
    minFrags = 1000)



project1 <- ArchRProject(
    ArrowFiles = Arrowfiles, outputDirectory = "C:/Users/Tom/Desktop/single_cell", copyArrows = TRUE
)


saveArchRProject(project1, outputDirectory = "C:/Users/Tom/Desktop/single_cell", load = FALSE)

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
#cluster_counts_harmony <- as.data.frame(t(as.data.frame(as.vector((table(project2$Clusters_Harmony))))))
#rownames(cluster_counts) <- NULL
#colnames(cluster_counts) <- names(table(project2$Clusters_Harmony)) 
 
cM_nobatchcorrection <- confusionMatrix(paste0(project2$Clusters_no_batch_correction), paste0(project2$Sample))
cM <- confusionMatrix(paste0(project2$Clusters_Harmony), paste0(project2$Sample))
 
#RNAseq data read in

RNAseq <- readRDS("seurat.pfc.final.rds")
RNAseq$cellIDs <- gsub('FC-', '', RNAseq$cellIDs)

#unconstrained integration

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
ExN <- paste0(c("ExN-1","ExN-3","ExN-5","ExN-2","ExN-4","ExN-6"), collapse = "|")
InN <- paste0(c("InN-1","InN-2","InN-3","InN-4"), collapse = "|")
RG <- paste0(c("RG-1","RG-2"), collapse = "|")
MG <- paste0("MG", collapse = "|")
Else <- paste0(c("IP","CycPro"), collapse="|")

clustExN <- rownames(cM)[grep(ExN, preClust)]
clustInN <- rownames(cM)[grep(InN, preClust)]
clustRG <- rownames(cM)[grep(RG, preClust)]
clustMG <- rownames(cM)[grep(MG, preClust)]

RNAcells_ExN <- colnames(RNAseq)[grep(pattern = ExN, x = RNAseq$cellIDs)]
RNAcells_InN <- colnames(RNAseq)[grep(pattern = InN, x = RNAseq$cellIDs)]
RNAcells_RG <- colnames(RNAseq)[grep(pattern = RG, x = RNAseq$cellIDs)]

groupList <- SimpleList(
     ExNlist = SimpleList(
         ATAC = project2$cellNames[project2$Clusters_Harmony %in% clustExN],
         RNA = RNAcells_ExN
     ),
     InNlist = SimpleList(
         ATAC = project2$cellNames[project2$Clusters_Harmony %in% clustInN],
         RNA = RNAcells_InN
     ),
     RGlist =  SimpleList(
         ATAC = project2$cellNames[project2$Clusters_Harmony %in% clustRG],
         RNA = RNAcells_RG	
     )
 )

project2 <- addGeneIntegrationMatrix(
    ArchRProj = project2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = RNAseq,
    addToArrow = FALSE,
    groupList = groupList,
    groupRNA = "cellIDs",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co",
	dimsToUse = 1:26,
	k.score = 26,
	npcs = 26,
	dims = 1:26
)

pal <- paletteDiscrete(values = RNAseq$cellIDs)

atac_and_rna_p1 <- plotEmbedding(
    project2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)

atac_and_rna_p2 <- plotEmbedding(
    project2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal
)

#Pseudo-bulk Replicates in ArchR

project3 <- addGeneIntegrationMatrix(
    ArchRProj = project2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = RNAseq,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "cellIDs",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

#changing cluster labels to cell types present in RNA seq data

cM2 <- confusionMatrix(project3$Clusters_Harmony, project3$predictedGroup)
labelOld <- rownames(cM2)
labelNew <- colnames(cM2)[apply(cM2, 1, which.max)]

project3$Clusters_RNAmapped <- mapLabels(project3$Clusters_Harmony, newLabels = labelNew, oldLabels = labelOld)

p1 <- plotEmbedding(project3, colorBy = "cellColData", name = "Clusters_RNAmapped")

#making pseudo bulk replicates

project4 <- addGroupCoverages(ArchRProj = project3, groupBy = "Clusters_RNAmapped")

#Calling peaks with tile matrix

project4@projectMetadata$outputDirectory <- "C:/Users/Tom/Desktop/single_cell"

project4 <- addReproduciblePeakSet(
    ArchRProj = project4, 
    groupBy = "Clusters_RNAmapped",
    peakMethod = "Tiles",
    method = "p",
    cutOff = 0.01,
    extendSummits = 500
)

getPeakSet(project4)
