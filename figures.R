`510-Pre-Filter-Metadata` <- readRDS("~/QualityControl/510/510-Pre-Filter-Metadata.rds")
`611-Pre-Filter-Metadata` <- readRDS("~/QualityControl/611/611-Pre-Filter-Metadata.rds")
`993-Pre-Filter-Metadata` <- readRDS("~/QualityControl/993/993-Pre-Filter-Metadata.rds")
table_precluster_A <- data.frame(matrix(ncol = 5, nrow = 3))
rownames(table_precluster_A) <- c('510', '611', '993')
colnames(table_precluster_A) <- c('Passed QC', 'Failed QC', 'Total Fragments', 'Median Fragments', 'Median TSS Enrichment')

#Trying to filter Keep data by 1's and 0's to determine passed QC and failed QC
df[1 %in% `510-Pre-Filter-Metadata`$Keep, ]
`510-Pre-Filter-Metadata`["Keep" ==1]
which(`510-Pre-Filter-Metadata`["Keep"] == 1)
`510-Pre-Filter-Metadata`(Keep)[names(`510-Pre-Filter-Metadata`(Keep)) == 1]

## Solution - note that I removed line 53 from the code sent in the email as it was superfluous

##  Inital ArchR QC -------------------------------------------------------------------
# ArchR does some QC when loading the files in so need to load the pre-QC info
# Pre-filter
for (SAMPLE in 1:length(SAMPLES)) {
  
  # Subset IDs
  sampleID <- substr(SAMPLES[SAMPLE], 1, 7)
  donorID <- substr(SAMPLES[SAMPLE], 1, 3)
  
  # Load Pre-filtered data
  preQC_df <- readRDS(paste0(WK_DIR, "QualityControl/", SAMPLES[SAMPLE], "/", 
                             SAMPLES[SAMPLE], "-Pre-Filter-Metadata.rds"))
  
  # Assign frag plots and counts_df
  assign(paste0("counts_df_", donorID), 
         data.frame("Sample" = sampleID,
                    "Cells_Pass_Filter" = sum(preQC_df$Keep),
                    "Cells_dropped" = sum(preQC_df$Keep == 0),
                    "Total_Frags" = sum(preQC_df$nFrags),
                    "Median_Frags" = median(preQC_df$nFrags[preQC_df$Keep == 1]),
                    "Median_TSS_Enrichment" = median(preQC_df$TSSEnrichment[preQC_df$Keep == 1])))
  
}

## Initial QC reporting  --------------------------------------------------------------
# Counts df
counts_df <- rbind(counts_df_510, counts_df_611, counts_df_993)
