rnaseq.dirs <- list.files('data/rnaseq')
print(rnaseq.dirs)

rnaseq.dirs[1]

raw_counts <- read.delim(paste0('./data/rnaseq/',  rnaseq.dirs[1], '/', list.files(paste0('data/rnaseq/', rnaseq.dirs[1]), pattern = '.tsv')), skip = 1)[-c(1:4), ]
raw_counts <- raw_counts[, c(1,4)]
head(raw_counts)
colnames(raw_counts)[2] <- rnaseq.dirs[1]

c = 3
for (i in rnaseq.dirs[-1]) {
    tmp <- read.delim(paste0('./data/rnaseq/',  i, '/', list.files(paste0('data/rnaseq/', i), pattern = '.tsv')), skip = 1)[-c(1:4), ]
    tmp <- tmp[, 4]
    raw_counts <- cbind(raw_counts, tmp)
    colnames(raw_counts)[c] <- i
    c = c+1
    print(c)
}
rm(tmp)
head(raw_counts)
dim(raw_counts)
rownames(raw_counts) <- raw_counts$gene_id
raw_counts <- as.matrix(raw_counts[, -1])
class(raw_counts)
is.numeric(raw_counts)
raw_counts[1:4, 1:4]

saveRDS(raw_counts, file = 'data/processed_data/TCGA_UCS_RNASEQ.rds')
