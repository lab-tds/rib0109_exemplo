setwd('~/Documents/rib0109/2024/projeto_exemplo/rib0109_exemplo/')

library(DESeq2)

raw_counts <- readRDS('data/processed_data/TCGA_UCS_RNASEQ_MATCHED.rds')
maf <- readRDS('data/processed_data/TCGA_UCS_MAF_MATCHED.rds')

maf <- as.data.frame(maf)
table(maf$TP53)

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = maf,
                              design = ~ TP53)

dds <- DESeq(dds)
res <- results(dds)

res
which(res$padj < 0.05)

top_de <- as.data.frame(res[which(res$padj < 0.05), ])

top_de[order(top_de$log2FoldChange),]

boxplot(log(assay(dds)['ENSG00000136750.13', ]+1) ~ maf$TP53)

top_de


rnaseq.dirs <- list.files('data/rnaseq')[1]
rnaseq.dirs

annot <- read.delim(paste0('./data/rnaseq/',  rnaseq.dirs[1], '/', list.files(paste0('data/rnaseq/', rnaseq.dirs[1]), pattern = '.tsv')), skip = 1)[-c(1:4), ]
annot <- annot[, c(1,2)]
head(annot)

top_de$gene_id <- rownames(top_de)
top_de <- dplyr::left_join(top_de, annot, by = 'gene_id')
head(top_de)

top_de <- top_de[order(abs(top_de$log2FoldChange), decreasing = TRUE),]
head(top_de)

saveRDS(top_de, file = 'data/results/DE_TP53_mut_vs_wt.rds')

png('data/results/TP53.png')
boxplot(log(assay(dds)['ENSG00000111049.4', ]+1) ~ maf$TP53, ylab = 'log(MYF5)', xlab = 'TP53 mut')
dev.off()
