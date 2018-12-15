library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(rtracklayer)
library(biomaRt)
library(rtracklayer)
library("org.Mm.eg.db")
library(GenomicRanges)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

# Read bed files
genes_00h <- import("../results/peaks/fold10_30_p1e6/00h_fold10_30_p1e6_peaks.bed",format = "BED")
genes_03h <- import("../results/peaks/fold10_30_p1e6/03h_fold10_30_p1e6_peaks.bed",format = "BED")
genes_06h <- import("../results/peaks/fold10_30_p1e6/06h_fold10_30_p1e6_peaks.bed",format = "BED")
genes_12h <- import("../results/peaks/fold10_30_p1e6/12h_fold10_30_p1e6_peaks.bed",format = "BED")

genes_00h_anno <- annotatePeak(genes_00h,tssRegion = c(-3000,3000),TxDb = txdb, level = "gene", annoDb = "org.Mm.eg.db", sameStrand = FALSE,ignoreOverlap = FALSE,overlap = "TSS")
genes_03h_anno <- annotatePeak(genes_03h,tssRegion = c(-3000,3000),TxDb = txdb, level = "gene", annoDb = "org.Mm.eg.db", sameStrand = FALSE,ignoreOverlap = FALSE,overlap = "TSS")
genes_06h_anno <- annotatePeak(genes_06h,tssRegion = c(-3000,3000),TxDb = txdb, level = "gene", annoDb = "org.Mm.eg.db", sameStrand = FALSE,ignoreOverlap = FALSE,overlap = "TSS")
genes_12h_anno <- annotatePeak(genes_12h,tssRegion = c(-3000,3000),TxDb = txdb, level = "gene", annoDb = "org.Mm.eg.db", sameStrand = FALSE,ignoreOverlap = FALSE,overlap = "TSS")

#Visualizations
plotAnnoPie(genes_00h_anno)
plotAnnoPie(genes_03h_anno)
plotAnnoPie(genes_06h_anno)
plotAnnoPie(genes_12h_anno)

upsetplot(genes_00h_anno,vennpie = TRUE)
upsetplot(genes_03h_anno,vennpie = TRUE)
upsetplot(genes_06h_anno,vennpie = TRUE)
upsetplot(genes_12h_anno,vennpie = TRUE)

head(as.data.frame(genes_00h_anno))
write.table(as.data.frame(genes_00h_anno), "../results/annotation/00h_peaks_anno.txt",row.names = F, col.names = T, sep = "\t", quote = F)
write.table(as.data.frame(genes_03h_anno), "../results/annotation/03h_peaks_anno.txt",row.names = F, col.names = T, sep = "\t", quote = F)
write.table(as.data.frame(genes_06h_anno), "../results/annotation/06h_peaks_anno.txt",row.names = F, col.names = T, sep = "\t", quote = F)
write.table(as.data.frame(genes_12h_anno), "../results/annotation/12h_peaks_anno.txt",row.names = F, col.names = T, sep = "\t", quote = F)

## PATHWAY ENRICHMENT
library(clusterProfiler)
ego_00h <- enrichGO(gene = as.data.frame(genes_00h_anno)$geneId,
                    "org.Mm.eg.db",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05
                    )
dotplot(ego_00h, showCategory = 20)

ego_03h <- enrichGO(gene = as.data.frame(genes_03h_anno)$geneId,
                    "org.Mm.eg.db",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05
)
dotplot(ego_03h, showCategory = 20)

ego_06h <- enrichGO(gene = as.data.frame(genes_06h_anno)$geneId,
                    "org.Mm.eg.db",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05
)
dotplot(ego_06h, showCategory = 20)

ego_12h <- enrichGO(gene = as.data.frame(genes_12h_anno)$geneId,
                    "org.Mm.eg.db",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05
)
dotplot(ego_12h, showCategory = 20)

# kegg pathway
kk_00h <- enrichKEGG(gene = as.data.frame(genes_03h_anno)$geneId,
                     organism = "mmu",
                     pvalueCutoff = 0.05)
dotplot(kk_00h, showCategory = 20)

kk_03h <- enrichKEGG(gene = as.data.frame(genes_03h_anno)$geneId,
                     organism = "mmu",
                     pvalueCutoff = 0.05)
dotplot(kk_03h, showCategory = 20)

kk_06h <- enrichKEGG(gene = as.data.frame(genes_06h_anno)$geneId,
                     organism = "mmu",
                     pvalueCutoff = 0.05)
dotplot(kk_06h, showCategory = 20)

kk_12h <- enrichKEGG(gene = as.data.frame(genes_12h_anno)$geneId,
                     organism = "mmu",
                     pvalueCutoff = 0.05)
dotplot(kk_12h, showCategory = 20)

