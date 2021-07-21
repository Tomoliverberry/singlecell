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
