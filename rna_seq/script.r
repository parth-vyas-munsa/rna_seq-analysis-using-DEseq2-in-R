library(DESeq2)
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
write.csv(counts , "counts.csv" , row.names = FALSE )
head(counts)

metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
write.csv(metadata , "metadata.csv" , row.names = FALSE )
head(metadata)



head(counts)
rownames(counts) = counts$Gene.ID
head(counts)


genes = counts[, c("Gene.ID" , "Gene.Name")]
counts = counts[ ,-c(1, 2)]
head(counts)

head(metadata)
rownames(metadata) = metadata$Run
head(metadata)


res_df = as.data.frame(res)
head(res_df)

head(genes)

res_df = merge(res_df , genes , by = "row.names")
head(res_df)



metadata = metadata[,c("Sample.Characteristic.genotype.") , drop = FALSE]
metadata

colnames(metadata) = c("genotype")
metadata

metadata$genotype[metadata$genotype == "wild type genotype"] = "wildtype"
metadata$genotype[metadata$genotype == "Snai1 knockout"] = 'knockout'
metadata

genes_id = genes$Gene.ID[genes$Gene.Name == 'SNAI1']
genes_count = counts[genes_id,]
genes_count

gene_data = cbind(metadata,counts = as.numeric(genes_count))
gene_data

library(ggplot2)
ggplot(gene_data,aes(x = genotype , y = counts , fill = genotype)) +geom_boxplot()






dds <- DESeqDataSetFromMatrix(countData = counts , colData = metadata , design=~ genotype)
dds

dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)
dds

res = results(dds,contrast = c("genotype" , "knockout" , "wildtype") , alpha = 1e-5)
res
check_genes = c("THY1" , "SFMBT2" , "PASD1")
res_df[res_df$Gene.Name %in% check_genes,]

sig <- subset(res,
              padj < 0.05 &
                abs(log2FoldChange) > 1)
head(sig)
upregulated <- subset(res,
                      padj < 0.05 &
                        log2FoldChange > 1)
downregulated <- subset(res,
                        padj < 0.05 &
                          log2FoldChange < -1)

nrow(upregulated)

nrow(downregulated)


plotMA(res)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res,lab = rownames(res) , x = "log2FoldChange" , y = "pvalue")

BiocManager::install("biomaRt")
library(biomaRt)
ensemble  <- useEnsembl(biomart = "genes")
datasets = listDatasets(ensemble)
head(datasets)


dataset_nb = grep("human" , datasets$description,ignore.case = TRUE )
dataset_nb

dataset = datasets$dataset[dataset_nb]
dataset

ensemble = useDataset(dataset = dataset , mart = ensemble)

attributes <- c("ensembl_gene_id" , "chromosome_name" , "start_position" , "end_position")
values <- list(ensembl_gene_id = c())
all.genes = getBM(attributes = attributes , values = values , mart = ensemble)
head(all.genes)


head(res_df)
colnames(all.genes)[1] = "Gene.ID"
head(all.genes)

merge_data = merge(all.genes , res_df , by = "Gene.ID")
head(merge_data)


merge_data$chromosome_name = paste("chr" , merge_data$chromosome_name , sep = "")
head(merge_data)
write.csv(merge_data, "Book2.csv" , row.names = FALSE)

merge_data_subset = merge_data[merge_data$Gene.Name %in% check_genes,]
head(merge_data_subset)

write.csv(merge_data_subset, "Book1.csv" , row.names = FALSE)





