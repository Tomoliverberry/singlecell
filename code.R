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
cluster_counts_harmony <- as.data.frame(t(as.data.frame(as.vector((table(project2$Clusters_Harmony))))))
rownames(cluster_counts_harmony) <- NULL
colnames(cluster_counts_harmony) <- names(table(project2$Clusters_Harmony)) 
 
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

clust_cM <- pheatmap::pheatmap(
     mat = as.matrix(cM),
     color = paletteContinuous("whiteBlue"),
     border_color = "black", display_numbers = TRUE, number_format = "%.0f"
 ) 

unique(unique(project2$predictedGroup_Un))
ExN <- paste0(c("ExN-1","ExN-3","ExN-5","ExN-2","ExN-4","ExN-6"), collapse = "|")
InN <- paste0(c("InN-1","InN-2","InN-3","InN-4"), collapse = "|")
RG <- paste0(c("RG-1","RG-2"), collapse = "|")
MG <- paste0("MG")


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
    method = "p"
)

getPeakSet(project4)

------------------------------------------------------------------------------------------------------------------------------
TABLES AND FIGURES
------------------------------------------------------------------------------------------------------------------------------

#PRE-CLUSTERING

#Exported plots as PNGs through RStudio

#Read in pre-filter quality control data generated by ArchR

`510-Pre-Filter-Metadata` <- readRDS("~/QualityControl/510/510-Pre-Filter-Metadata.rds")
`611-Pre-Filter-Metadata` <- readRDS("~/QualityControl/611/611-Pre-Filter-Metadata.rds")
`993-Pre-Filter-Metadata` <- readRDS("~/QualityControl/993/993-Pre-Filter-Metadata.rds")

#Table 1: Preclustering QC summary for donors---------------------------------------------------------------------------------

table_precluster_A <- data.frame(
"Sample" = "510",
"Cells_Passed_QC" = sum(`510-Pre-Filter-Metadata`$Keep),
"Cells_Failed_QC" = sum(`510-Pre-Filter-Metadata`$Keep == 0),
"Total_Frags" = sum(`510-Pre-Filter-Metadata`$nFrags),
"Median_Frags" = median(`510-Pre-Filter-Metadata`$nFrags[`510-Pre-Filter-Metadata`$Keep ==1]),
"Median_TSS_Enrichment" = median(`510-Pre-Filter-Metadata`$TSSEnrichment[`510-Pre-Filter-Metadata`$Keep == 1])
)

newrow_precluster_A <- data.frame(
"Sample" = "611",
"Cells_Passed_QC" = sum(`611-Pre-Filter-Metadata`$Keep),
"Cells_Failed_QC" = sum(`611-Pre-Filter-Metadata`$Keep == 0),
"Total_Frags" = sum(`611-Pre-Filter-Metadata`$nFrags),
"Median_Frags" = median(`611-Pre-Filter-Metadata`$nFrags[`611-Pre-Filter-Metadata`$Keep ==1]),
"Median_TSS_Enrichment" = median(`611-Pre-Filter-Metadata`$TSSEnrichment[`611-Pre-Filter-Metadata`$Keep == 1])
)
table_precluster_A <- rbind(table_precluster_A, newrow_precluster_A)

newrow_precluster_A <- data.frame(
"Sample" = "993",
"Cells_Passed_QC" = sum(`993-Pre-Filter-Metadata`$Keep),
"Cells_Failed_QC" = sum(`993-Pre-Filter-Metadata`$Keep == 0),
"Total_Frags" = sum(`993-Pre-Filter-Metadata`$nFrags),
"Median_Frags" = median(`993-Pre-Filter-Metadata`$nFrags[`993-Pre-Filter-Metadata`$Keep ==1]),
"Median_TSS_Enrichment" = median(`993-Pre-Filter-Metadata`$TSSEnrichment[`993-Pre-Filter-Metadata`$Keep == 1])
)
table_precluster_A <- rbind(table_precluster_A, newrow_precluster_A)

#Table 2: Count table for doublet removal-------------------------------------------------------------------------------------

`510-Doublet-Summary` <- readRDS("C:/Users/Tom/Desktop/single_cell/QualityControl/510/510-Doublet-Summary.rds")
`611-Doublet-Summary` <- readRDS("C:/Users/Tom/Desktop/single_cell/QualityControl/510/510-Doublet-Summary.rds")
`993-Doublet-Summary` <- readRDS("C:/Users/Tom/Desktop/single_cell/QualityControl/510/510-Doublet-Summary.rds")

table_doubletremoval_B <- data.frame(
"Sample" = c("510", "611", "993"),
"Cell_Count_Pre_Doublet_Removal" = c("8321","5029","2287"),
"Cell_Count_Post_Doublet_Removal" = c("7629","5029","2235"),
"Percentage_of_Cells_Removed_per_Donor" = c("8.3%","0%","2.3%")
)

#Plot 1: TSS Enrichment vs log10(number of unique fragments)------------------------------------------------------------------

df <- getCellColData(project2, select = c("log10(nFrags)", "TSSEnrichment"))

df611 <- df[grep("611", rownames(df)), ]
df510 <- df[grep("510", rownames(df)), ]
df993 <- df[grep("993", rownames(df)), ]

plot611_1 <- ggPoint(
     x = df611[,1], 
     y = df611[,2], 
     colorDensity = TRUE,
     continuousSet = "sambaNight",
     xlabel = "Log10 Unique Fragments",
     ylabel = "TSS Enrichment",
     xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
     ylim = c(0, quantile(df[,2], probs = 0.99))
 ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
 
 plot510_1 <- ggPoint(
    x = df510[,1], 
    y = df510[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plot993_1 <- ggPoint(
    x = df993[,1], 
    y = df993[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

#Ridge plot for TSS enrichment per donor--------------------------------------------------------------------------------------

#Batch corrected data used for pre clustering with filtering

Ridge_plot <- plotGroups(
     ArchRProj = project2,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "TSSEnrichment",
     plotAs = "ridges"
 )
 
#Fragment size plot-----------------------------------------------------------------------------------------------------------

plot_fragmentsizes <- plotFragmentSizes(project2, groupBy = "Sample")

#POST-CLUSTERING--------------------------------------------------------------------------------------------------------------

#Cell counts per cluster per donor--------------------------------------------------------------------------------------------

table_cell_per_cluster <- as.matrix(confusionMatrix(project2$Sample, project2$Clusters_Harmony))

#Cell counts per cluster------------------------------------------------------------------------------------------------------

#Janitor package will be used for data frame editing

install.packages("janitor")
library(janitor)

#Make data frame copy of table containing number of cells per cluster per donor

table_cell_cluster_all <- as.data.frame(table_cell_per_cluster)

#Move rownames to column as use of adorn_totals creates a "Total" cell which shifts data

setDT(table_cell_cluster_all, keep.rownames = TRUE)[]

#Use janitor feature to apply new row to data frame containing totals

table_cell_cluster_all %<>% adorn_totals(dat = table_cell_cluster_all, where = "row")

#Remove rows containing donor specific data and column containing "Total" cell

table_cell_cluster_all <- table_cell_cluster_all[-c(1,2,3), ]
table_cell_cluster_all <- table_cell_cluster_all[,-1]

#UMAPs------------------------------------------------------------------------------------------------------------------------
#IterativeLSI UMAP------------------------------------------------------------------------------------------------------------
project2 <- addUMAP(
    ArchRProj = project2,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)
plot_UMAP_by_cluster <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Clusters_no_batch_correction", embedding = "UMAP")
plot_UMAP_by_donor <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#Harmony UMAP-----------------------------------------------------------------------------------------------------------------

project2 <- addUMAP(
    ArchRProj = project2,
    reducedDims = "Harmony",
    name = "UMAP_Harmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)
plot_UMAP_by_cluster_harmony <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Clusters_Harmony", embedding = "UMAP_Harmony")
plot_UMAP_by_donor_harmony <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony")

 
 
 
#Final RNA cell type mapping to atac-seq data
 
 atac_and_rna_p3 <- plotEmbedding(
     project4, 
     colorBy = "cellColData",
     embedding = "UMAP_Harmony",
     name = "predictedGroup", 
     pal = pal
 )
 
table_cell_cluster_RNA_ATAC <- as.matrix(confusionMatrix(project4$Sample, project4$Clusters_RNAmapped))


------------------------------------------------------------------------------------------------------------------------------
FINEMAPPING
------------------------------------------------------------------------------------------------------------------------------

#Load in rsnps package, to obtain hg38 information from hg19 data by accessing OpenSNP and NCBI's dbSNP SNP database
install.packages("rsnps")
library(rsnps)

#Load readxl package for PGC3 SNP data
install.packages("readxl")
library(readxl)

#Read .xlsx file containing PGC3 data
SNPs_PGC3 <- read_excel("Copy of PGC3 revised 11-8-21 finemapped SNPs with at least 0.1 posterior probability(1335).xlsx")

#rsnps requires the RSID's to be in a character string. Extract RSID's to new vector
RSID <- as.list(SNPs_PGC3$rsid)

#Perform NCBI search of all RSID's. Package limits to one SNP per second.
snp_query <- ncbi_snp_query(RSID)

#Read in all peakcalling data from ArchR for each RNA-seq cell type cluster
peak_ExN_2 <- readRDS(file = "~/PeakCalls/ExN.2-reproduciblePeaks.gr.rds")
peak_ExN_3 <- readRDS(file = "~/PeakCalls/ExN.3-reproduciblePeaks.gr.rds")
peak_ExN_4 <- readRDS(file = "~/PeakCalls/ExN.4-reproduciblePeaks.gr.rds")
peak_ExN_5 <- readRDS(file = "~/PeakCalls/ExN.5-reproduciblePeaks.gr.rds")
peak_ExN_6 <- readRDS(file = "~/PeakCalls/ExN.6-reproduciblePeaks.gr.rds")
peak_InN_1 <- readRDS(file = "~/PeakCalls/InN.1-reproduciblePeaks.gr.rds")
peak_InN_2 <- readRDS(file = "~/PeakCalls/InN.2-reproduciblePeaks.gr.rds")
peak_InN_3 <- readRDS(file = "~/PeakCalls/InN.3-reproduciblePeaks.gr.rds")
peak_RG_1 <- readRDS(file = "~/PeakCalls/RG.1-reproduciblePeaks.gr.rds")
peak_RG_2 <- readRDS(file = "~/PeakCalls/RG.2-reproduciblePeaks.gr.rds")


SNPs_h19_h38 <- SNPs_PGC3
SNPs_h19_h38 <- cbind(SNPs_h19_h38, hg38_position = snp_query$bp)
SNPs_h19_h38 <- cbind(SNPs_h19_h38, gene = snp_query$gene)

rsID_in_peak_range <- cbind(SNPs_h19_h38[FALSE,], "Cell_type" = character(), "Peak_Start" = integer(), "Peak_Stop" = integer())
rsID_not_in_peak_range <- cbind(SNPs_h19_h38[FALSE,], "Cell_type" = character(), "Peak_Start" = integer(), "Peak_Stop" = integer())

data_names <- c("ExN_2", "ExN_3", "ExN_4", "ExN_5", "ExN_6", "InN_1", "InN_2", "InN_3", "RG_1", "RG_2")

 for (k in data_names) {
     peak_CELLTYPE <- eval(parse(text = paste0("peak_",k,"@ranges@start")))
     for(i in 1:nrow(SNPs_h19_h38)) {
         break_check = 0
         for(j in peak_CELLTYPE) {
             if(j > SNPs_h19_h38[i, 21]){
                 #print(SNPs_h19_h38[i, ])
                 rsID_not_in_peak_range <- add_row(rsID_not_in_peak_range, SNPs_h19_h38[i, ], "Cell_type" = k, "Peak_Start" = j, "Peak_Stop" = (j+500))	
                 break_check = 1
                 break
             } else if (SNPs_h19_h38[i, 21] <= j+500){
                 rsID_in_peak_range <- add_row(rsID_in_peak_range, SNPs_h19_h38[i, ], "Cell_type" = k, "Peak_Start" = j, "Peak_Stop" = (j+500))	
                 break_check = 1
                 break
             }
         }#end of j
         if(break_check == 0){
             rsID_not_in_peak_range <- add_row(rsID_not_in_peak_range, SNPs_h19_h38[i, ], "Cell_type" = k, "Peak_Start" = j, "Peak_Stop" = (j+500))	
         }
     } #end of i
 }#end k
-------------------------------------------------------------------------------------------------------------------------------
