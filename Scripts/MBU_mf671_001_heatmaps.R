#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()

library("ggplot2")
library("ggrepel")
library("cowplot")
library("DESeq2")
library("biomaRt")
library("dplyr")
library("RColorBrewer")
library("KEGGREST")
library("ComplexHeatmap")
library("circlize")


Project        <- "MBU_mf671_001_"
baseDir        <- "/Users/mn367/Documents/MBU-Projects/MBU_Michele_Frison/MBU_mf671_001"
setwd(baseDir)

group_cols <-  c("lowHeteroplasmy"="#107A78", "highHeteroplasmy"="#A7740E" )

l2fc_cutoff <- 0.6


message("+-------------------------------------------------------------------------------")
message("+                       load in critical objects                                ")
message("+-------------------------------------------------------------------------------")

sampleTable   <- read.csv("Input/sampleTable_high_glucose.csv", row.names = 1)
sampleTable$Heteroplasmy_group <- factor(sampleTable$Heteroplasmy_group, levels = c("lowHeteroplasmy",  "highHeteroplasmy"))
sampleTable_coll <- sampleTable[!duplicated(sampleTable$Barcode) , -c(1,9)]

countData_filt <- read.csv("Input/raw_counts.csv", row.names = 1)



message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)


message("+-------------------------------------------------------------------------------")
message("+                 DESeq2                                ")
message("+-------------------------------------------------------------------------------")

dds = DESeqDataSetFromMatrix(countData = countData_filt, colData = sampleTable, design = ~ Heteroplasmy_group )
dds <- collapseReplicates(dds, groupby = sampleTable$Barcode)
dds = DESeq(dds, minReplicatesForReplace=5)


cbind(resultsNames(dds))

res <- results(dds, name = "Heteroplasmy_group_highHeteroplasmy_vs_lowHeteroplasmy", alpha = 0.05)
head(res)
summary(res)


rld = rlog(dds)
vsd = vst(dds)


res_highGlucose.ann <- as.data.frame(res[order(res$padj),])
res_highGlucose.ann$external_gene_name <- ensEMBL2id[match( rownames(res_highGlucose.ann), ensEMBL2id$ensembl_gene_id),]$external_gene_name
res_highGlucose.ann <- na.omit(res_highGlucose.ann)
resSig_highGlucose <- res_highGlucose.ann[res_highGlucose.ann$padj < 0.05,]

rld_df <- as.data.frame(assay(rld))

#write.csv(res_highGlucose.ann, paste(Project, "results_highGlucose_m5024_MEFs.csv", sep = "_"), quote = FALSE)
#write.csv(resSig_highGlucose, paste(Project, "resSig_highGlucose_m5024_MEFs_DEGs.csv", sep = "_"), quote = FALSE)
#saveRDS(rld, paste(Project, "rld.Rds", sep = "_"))



message("+-------------------------------------------------------------------------------")
message("+               ubiquitin-mediated proteolysis genes                            ")
message("+-------------------------------------------------------------------------------")


df <- data.frame(genes=keggGet("mmu04120")[[1]]$GENE)
df %>% filter(row_number() %% 2 == 0) -> df ## Select even rows
genes <- gsub( ";.*", "", df$genes)

UbiProt_genes <- unique(genes)

resSig_highGlucose[resSig_highGlucose$external_gene_name %in% UbiProt_genes ,]


### proteasome genes

df <- data.frame(genes=keggGet("mmu03050")[[1]]$GENE)
df %>% filter(row_number() %% 2 == 0) -> df ## Select even rows
genes <- gsub( ";.*", "", df$genes)

Proteosome_genes <- unique(genes)

resSig_highGlucose[resSig_highGlucose$external_gene_name %in% Proteosome_genes ,]




message("+-------------------------------------------------------------------------------")
message("+             PLOT HEATMAP with  ubiquitin-mediated proteolysis                 ")
message("+-------------------------------------------------------------------------------")

rld_df_meanCentered <- rld_df - rowMeans(rld_df)   



# choose selected genes :
selected_genes <- UbiProt_genes
#selected_genes <- UbiProt_genes[UbiProt_genes %in% resSig_highGlucose$external_gene_name]
#selected_genes <- c("Ubb", "Ubc", "Uba52", "Rps27a")
selected_genes <- Proteosome_genes
selected_genes <- Proteosome_genes[Proteosome_genes %in% resSig_highGlucose$external_gene_name]



mat <- rld_df_meanCentered[rownames(rld_df_meanCentered) %in% ensEMBL2id[ ensEMBL2id$external_gene_name %in% selected_genes,]$ensembl_gene_id ,]
mat$gene <- ensEMBL2id[match(rownames(mat) , ensEMBL2id$ensembl_gene_id),]$external_gene_name
rownames(mat) <- mat$gene
mat <- mat[,-ncol(mat)]

sampleTable_coll <- sampleTable_coll[order(sampleTable_coll$Heteroplasmy2),]
mat2 <- mat[, match(sampleTable_coll$Barcode, colnames(mat))]


ht_anno <- data.frame(sample= colnames(mat2))
ht_anno$Heteroplasmy <- sampleTable_coll[match(ht_anno$sample, sampleTable_coll$Barcode),]$Heteroplasmy2
ht_anno$group <- sampleTable_coll[match(ht_anno$sample, sampleTable_coll$Barcode),]$Heteroplasmy_group
ht_anno <- ht_anno[order(ht_anno$Heteroplasmy),]

ha_top = HeatmapAnnotation(Heteroplasmy=anno_points(ht_anno$Heteroplasmy, pch = 16, size = unit(0.8, "mm")),  show_annotation_name = TRUE)
ha_bottom = HeatmapAnnotation(Group = ht_anno$group, col = list(Group = group_cols),  show_annotation_name = TRUE)

mat2 <- mat2[,match( ht_anno$sample, colnames(mat2))]
dim(mat2)

col_SIG_pseudo <- ifelse( rownames(mat2) %in%  resSig_highGlucose$external_gene_name, "green4", "black")


ht1 = Heatmap(as.matrix(mat2),  name = "Bulk RNA-seq",  row_title = "", column_title = "Ubiquitin mediated proteolysis genes", 
              show_row_names = TRUE, cluster_columns = FALSE, cluster_rows = TRUE , 
              heatmap_legend_param = list(title = "Expression in MEFs", legend_height = unit(3, "cm"), title_position = "topleft"), 
              row_names_side ="right", column_split = ht_anno$condition, 
              bottom_annotation = ha_bottom, top_annotation = ha_top, 
              row_names_gp = gpar(fontface = "italic", col = col_SIG_pseudo))


ht1

pdf(paste(Project, "ComplexHeatmap",  "5024_MEFs", "", "ubiquitin_mediated_proteolysis_genes", "all", ".pdf", sep="_"), onefile=FALSE, width=10, height=20)
#pdf(paste(Project, "ComplexHeatmap",  "5024_MEFs", "", "ubiquitin_mediated_proteolysis_genes", "DEGs", ".pdf", sep="_"), onefile=FALSE, width=7, height=4.5)
#pdf(paste(Project, "ComplexHeatmap",  "5024_MEFs", "", "ubiquitin_genes", "all", ".pdf", sep="_"), onefile=FALSE, width=7, height=2.5)
#pdf(paste(Project, "ComplexHeatmap",  "5024_MEFs", "", "proteosome_genes", "all", ".pdf", sep="_"), onefile=FALSE, width=7, height=9)
#pdf(paste(Project, "ComplexHeatmap",  "5024_MEFs", "", "proteosome_genes", "DEGs", ".pdf", sep="_"), onefile=FALSE, width=7, height=2)
par(bg=NA)
draw(ht1, ht_gap = unit(1.5, "cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()






