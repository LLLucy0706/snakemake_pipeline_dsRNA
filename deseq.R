library("tximport")
library("readr")
setwd("~/Desktop")

#locate the directory containing salmon files
dir <- "~/Desktop"
files <- file.path(dir, "salmon",list.files("salmon"),"quant.sf")
names(files) <- paste0("sample", 1:9)

#load transcript-to-gene file
#names(mydata) <- c("TXNAME", "GENEID")
#write.csv(mydata, file = "transcript_to_gene.csv", row.names=FALSE)
tx2gene <- read_csv(file.path(dir, "transcript_to_gene.csv"))

#import transcript-level estimates and summarized to the gene-level
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)



library("DESeq2")
#experiment design
sampleTable <- data.frame(condition = factor(rep(c("control", "sprayed_petiole","unsprayed_petiole"),
                                                 each=3)))
rownames(sampleTable) <- colnames(txi$counts)

#pass txi to DESeqDataSetFromTximport to import gene-level estimated counts
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

#filtering
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

#Run the DESeq function
dds <- DESeq(dds)
res <- results(dds)

#control vs sprayed_petiole
res1 <- results(dds, contrast = c("condition","sprayed_petiole","control"))

#control vs unsprayed_petiole
res2 <- results(dds, contrast = c("condition","unsprayed_petiole","control"))

#get differential expression results and order by adjusted p-value
res1Sig <- subset(res1, padj < 0.05 & abs(log2FoldChange) > 1)
res1Sig <- res1Sig[order(res1Sig$padj), ]
res2Sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 1)
res2Sig <- res2Sig[order(res2Sig$padj), ]

#sigGene <- intersect(c(rownames(res1Sig)),c(rownames(res2Sig)))

#write results
write.csv(res1Sig, file="control_vs_sprayed_petiole.csv")
write.csv(res2Sig, file="control_vs_unsprayed_petiole.csv")


#volcano plot control vs sprayed_petiole
png("volcanoplot1.png", 1200, 1000, pointsize=20)
with(res1, plot(log2FoldChange, -log10(padj), pch=20, main="control vs sprayed_petiole", xlim=c(-3,3)))
with(subset(res1, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(v=c(-1,1))
abline(h=-log10(0.05))
dev.off()

#volcano plot control vs unsprayed_petiole
png("volcanoplot2.png", 1200, 1000, pointsize=20)
with(res2, plot(log2FoldChange, -log10(padj), pch=20, main="control vs unsprayed_petiole", xlim=c(-3,3)))
with(subset(res2, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(v=c(-1,1))
abline(h=-log10(0.05))
dev.off()
