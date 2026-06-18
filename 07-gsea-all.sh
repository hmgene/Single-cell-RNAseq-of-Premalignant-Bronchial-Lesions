library(gprofiler2)
library(clusterProfiler)
library(fgsea)
library(data.table)
library(org.Hs.eg.db)

j=c(0, 1, 2, 14, 15, 11, 17)
edger <- fread("data/s13_selfrenewal_DE.csv")
edger = edger[ cluster %in% j]
sig <- edger[padj< 0.05 & abs(logFC) > 1]
genes <- unique(sig$gene)   # adjust column name if needed
gene_df <- bitr( genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_genes <- unique(gene_df$ENTREZID)
ego_bp<-enrichGO(gene=entrez_genes,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.25,readable=TRUE)
dotplot(ego_bp, showCategory = 20) + ggplot2::ggtitle("GO Biological Process Enrichment")

p=barplot(ego_bp, showCategory = 20)
ggsave(p,file="figures/gobp-s3s4epi-selfrenewal-de-barplot.svg",height=15)

p=cnetplot(ego_bp, showCategory = 5)
ggsave(p,file="figures/gobp-s3s4epi-selfrenewal-de-network.svg")

heatplot(ego_bp, showCategory = 10)

## Enricher

input_gmt="data/h.all.v2026.1.Hs.symbols.gmt" 
gmt_files <- Sys.glob(input_gmt)
gmt_df <- rbindlist(lapply(lapply(gmt_files, fgsea::gmtPathways),
    \(x) data.table(term = rep(names(x), lengths(x)), gene = unlist(x, use.names = FALSE)
)))

edger <- fread("data/s3s4_selfrenewal_DE_core_epi.csv")
sig <- edger[padj< 0.05 & abs(logFC) > 1]
genes <- unique(sig$gene)   # adjust column name if needed
res <- enricher( gene = genes, TERM2GENE = gmt_df[, .(term, gene)])

df <- as.data.table(res_df)
df <- df[order(p.adjust)]

ggplot(df, aes( x = reorder(Description, -p.adjust), y = -log10(p.adjust), size = Count, color = FoldEnrichment)) +
  geom_point() + coord_flip() + theme_bw() +
  labs( title = "Hallmark Enrichment (edgeR DE genes)",
    x = "Pathway", y = "-log10 adjusted p-value")

## GSEA

tt=fread("data/gsea_s3s4epi_h.all.tsv")
s3s4_up=tt[ES > 0][head(order(pval), n=10), pathway]
s3s4_dn=tt[ES < 0][head(order(pval), n=10), pathway]
s3s4_pathways = c(s3s4_up,s3s4_dn)


edger=fread("data/s13_selfrenewal_DE.csv") 
edger=fread("data/s3s4_selfrenewal_DE.csv") 
j=c(0, 1, 2, 14, 15, 11, 17)
edger=edger[cluster%in%j]

gmt_files <- "data/h.all.v2026.1.Hs.symbols.gmt"
gmt_list <- fgsea::gmtPathways(gmt_files)

edger <- edger[!is.na(logFC)]

ranks_df <- edger[, .(logFC = mean(logFC, na.rm = TRUE)), by = gene]
ranks <- ranks_df$logFC
names(ranks) <- ranks_df$gene
ranks <- sort(ranks, decreasing = TRUE)

fgseaRes <- fgseaMultilevel( pathways = gmt_list, stats = ranks, minSize = 15, maxSize = 500,nPermSimple = 1000 )
fgseaRes = fgseaRes[ !is.na(padj) ] 
fgseaRes <- fgseaRes[order(padj)]
fwrite(fgseaRes,file="data/gsea_s13epi_h.all.tsv")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

p=plotGseaTable(gmt_list[topPathways], ranks, fgseaRes, gseaParam=0.5)
ggsave(p,file="figures/gsea_s13epi.svg")

p=plotGseaTable(gmt_list[s3s4_pathways], ranks, fgseaRes, gseaParam=0.5)
ggsave(p,file="figures/gsea_s13epi_s3s4path.svg")

