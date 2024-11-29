
### RNAseq

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


### MAF

maf.dirs <- list.files('data/maf')
print(maf.dirs)


maf <- read.delim(paste0('./data/maf/',  maf.dirs[1], '/', list.files(paste0('data/maf/', maf.dirs[1]), pattern = 'maf.gz')), skip = 7)
head(maf)
colnames(maf)
maf <- maf[, c(1,9,16)]
head(maf)

for (i in maf.dirs[-1]) {
    tmp <- read.delim(paste0('./data/maf/',  i, '/', list.files(paste0('data/maf/', i), pattern = 'maf.gz')), skip = 7)
    tmp <- tmp[, c(1,9,16)]
    maf <- rbind(maf, tmp)
    print(i)
}
head(maf)
barplot(table(maf$Tumor_Sample_Barcode)[order(table(maf$Tumor_Sample_Barcode))])
maf$Tumor_Sample_Barcode[grep('POLE', maf$Hugo_Symbol)]

table(maf$Variant_Classification)
table(maf$Tumor_Sample_Barcode, maf$Hugo_Symbol)

muts <- c('Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation')

saveRDS(maf, file = './data/processed_data/TCGA_UCS_MAF.rds')

maf <- maf[which(maf$Variant_Classification %in% muts), ]
table(maf$Variant_Classification)

saveRDS(maf, file = './data/processed_data/TCGA_UCS_MAF_filterd.rds')
raw_counts[1:4,1:4]
head(maf)


#### Integrar dados

manifest <- read.delim('manifests/gdc_manifest.2024-11-29.084914.txt.map2submitterID')
head(manifest)
manifest$cases.0.samples.0.submitter_id[which(manifest$id == '0f638551-a5b7-4765-841e-23ae1424f14c')]

rownames(manifest) <- manifest$id
manifest <- manifest[colnames(raw_counts), ]
all(colnames(raw_counts) == manifest$id)

colnames(raw_counts) <- manifest$cases.0.samples.0.submitter_id

raw_counts[1:4, 1:4]

head(maf)
maf$barcode <- stringr::str_sub(maf$Tumor_Sample_Barcode, 1,16)

maf$Hugo_Symbol <- as.character(as.factor(maf$Hugo_Symbol))
names(table(maf$Hugo_Symbol)[order(table(maf$Hugo_Symbol), decreasing = TRUE)][1:10])

maf_f <- maf[which(maf$Hugo_Symbol %in% names(table(maf$Hugo_Symbol)[order(table(maf$Hugo_Symbol), decreasing = TRUE)][1:10])), ]
maf_table <- table(maf_f$barcode, maf_f$Hugo_Symbol)

head(maf_table)

dim(maf_table)
dim(raw_counts)


intersect(rownames(maf_table), colnames(raw_counts))
mis <- colnames(raw_counts)[which(!colnames(raw_counts) %in% rownames(maf_table))]

dim(maf_table)
maf_table <- rbind(maf_table, rep(0, ncol(maf_table)))
rownames(maf_table)[c(56,57)] <- mis
maf_table

maf_table <- maf_table[colnames(raw_counts), ]

all(rownames(maf_table) == colnames(raw()))
maf_table[maf_table > 1] <- 1


saveRDS(raw_counts, file = 'data/processed_data/TCGA_UCS_RNASEQ_MATCHED.rds')
saveRDS(maf_table, file = 'data/processed_data/TCGA_UCS_MAF_MATCHED.rds')
