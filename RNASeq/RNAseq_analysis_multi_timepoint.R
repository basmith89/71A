library("DESeq2")
library("vsn")
library("pheatmap")
library("DEGreport")
library("dplyr")
library("tibble")

#if bioconductor packages haven't been installed already use
# BiocManager::install("package_name")

wdir <- "setwd"
cts <- as.matrix(read.csv("DESeqCountData.csv", sep=",", row.names = "Geneid"))
coldata <- read.csv("condition.csv", row.names = 1)

coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

#i think ideally it should look something like this and using both new and old data if it is somewhat similar
#dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ batch + condition + time + condition:time)
#dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ batch + treatment + time)
#obviously the prefiltering and releveling steps would be prior to running the analysis step

#can call dds to inspect summarized info

#prefiltering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#relevel the factor to specify control group -- will become critical if comparing both data sets together
dds$condition <- relevel(dds$condition, ref = "control")

#difexp analysis for time series
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds)
#it's important to note that when using the LRT model each comparison is stored in your DESEQ analysis object (dds)
#so, when you head() your results you will only see one comparison ie. T1vsT3
#but you can call each comparison by name
#first get names by resultsNames(dds) then specify the name when calling results results(dds, name = "name_of_interest")



#If doing pairwise comparisons (Wald's test) I would shrink the data here
#resLFC <- lfcShrink(dds......

#plot counts by gene. this takes the gene with the smallest p value
#This code can be commented out generally unless interested
d <- plotCounts(dds, gene=which.min(res$padj), intgroup = "condition",
                returnData = TRUE)
ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1, h=0)) +
  scale_y_log10(breaks=c(25,100,400))

#transfrom data
vsd <- vst(dds)
rld <- rlog(dds)

#SD rows over means
#Note this figure may not be beneficial for experiments where there are many expected differences
meanSdPlot(assay(vsd))

#store pca
pca_plot <- plotPCA(vsd)

#plot MA
plotMA(res, ylim=c(-2,2))

##Joshua's filtering step
##I tend to disagree with running this with the time series experiment because we are running the LRT model while this reverts back to the Wald's test 
#resGA <- results(dds, lfcThreshold = 2, altHypothesis = "greaterAbs")

#order results by padj
resOrdered <- res[order(res$padj, res$log2FoldChange),]

#output data and images


#DEG LRT Gene Clustering Analysis
#reformat data to tibble for clustering analysis
res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

#filter padj significance
sigLRT_gene <- res_tb %>%
  dplyr::filter(padj < 0.05)

#transform counts for data vis
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

#clustered hierarchical heatmap (QC)
pheatmap(rld_cor, annotation = coldata)

#adjust data by padj --- could subset data if necessary here
clustering_sig_genes <- sigLRT_gene %>%
  arrange(padj)

#grab rlog values for sig genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene,]

#plot gene clusters across all sample groups
clusters <- degPatterns(cluster_rlog, metadata = coldata, time = "condition", col=NULL)


#inspecting the list
# names(clusters)
# head(clusters$df)

#extract group of interest
group1 <- clusters$df %>%
  dplyr::filter(cluster == 1)

#remove row name column
rownames(group1) <- NULL

cluster_count <- 1

#for loop of above
for (i in unique(clusters$df$cluster)) {
  assign(paste0("group", i), 
         clusters$df %>% dplyr::filter(cluster == i) %>% 
           `rownames<-`(NULL))
}

###Data manipulation for downstream analysis
#for loops to pull the data frame for each comparison in results (control vs T3) and filter out gene cluster groups.
result_names <- resultsNames(dds)
result_names <- setdiff(result_names, "Intercept")

results_list <- list()

for (name in result_names) {
  results_list[[name]] <- results(dds, name = name)
}

#extract time series lfc for funage pro
log2fc_df <- data.frame(row.names = Reduce(intersect, lapply(results_list, rownames)))

for (name in names(results_list)) {
  log2fc_df[[name]] <- results_list[[name]][rownames(log2fc_df), "log2FoldChange"]
}

#run data through funage or any pathway analysis software
#write lfc df to csv for downstream pathway analysis
write.csv(log2fc_df, "funagepro_input.csv")

#match cluster results and extract by designated group
cluster_match_results <- list()

for (name in names(results_list)) {
  cluster_match_results[[name]] <- results_list[[name]][rownames(results_list[[name]]) %in% group5$genes, ] 
}

for (n in names(cluster_match_results)) {
  write.csv(cluster_match_results[[n]], paste0("pwd",n, "_group#.csv"))
}