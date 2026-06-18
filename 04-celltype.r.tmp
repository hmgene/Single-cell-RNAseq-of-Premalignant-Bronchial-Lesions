input="bigdata/s13_clustered.rds"
input_info="data/sample_input.csv"
input_marker="data/Marker Genes for Annotation.xlsx";
input_pathology="data/Pathology and sample ID.xlsx";
output="bigdata/s13_celltyped.rds"

library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(readxl)

d=readRDS(input);
#devtools::install_github("satijalab/AzimuthAPI")
#d <- CloudAzimuth(d);
#saveRDS(d,file=output)

## dotplot moumita markers
mk = as.data.table(read_excel(input_marker, sheet = "Pan Marker Gene"))
#p <- DotPlot( d, features = split(mk$Gene,mk$Annotation),group.by = "seurat_clusters") +
#  scale_color_viridis_c() + theme_classic() + theme( axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
p <- DotPlot( d, features = split(mk$Gene, mk$Annotation), group.by = "seurat_clusters") +
  scale_color_viridis_c() +
  theme_classic() + theme( axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 90, face = "bold"))
ggsave(p,file="figures/dotplot_moumita_markers.svg",width=21,height=7)
ggsave(p,file="figures/dotplot_moumita_markers.pdf",width=21,height=7)


## pathology
pt=as.data.table(read_excel(input_pathology))
d@meta.data$Description <- pt$Description[match(d@meta.data$orig.ident, pt$id)]
d$DysRep = paste0(d$Description,d$selfrenewal)



mk=list();
mk[["barry"]]=fread("data/barry2020_long.txt")
mk[["moumita"]]=fread("data/moumita2022_long.txt")
mk[["moumita"]]$feature <- toupper(mk[["moumita"]]$feature)

## Azimuth
#https://satijalab.org/pan_human_azimuth/
#devtools::install_github("satijalab/AzimuthAPI")
#d <- CloudAzimuth(d);
#meta <- data.frame(cell_id = colnames(d), d@meta.data[, c(
#    "full_hierarchical_labels","azimuth_fine","selfrenewal", "final_level_labels",
#    "azimuth_label", "final_level_confidence", "full_consistent_hierarchy"
#)])
#write.csv(meta, "data/meta_azimuth.csv", row.names = FALSE)



x=table( d$full_hierarchical_labels, d$seurat_clusters)
x = x[apply(x, 1, function(r) sum(r)>200 & max(r) > 0),]
x=scale(x)

library(ComplexHeatmap)
library(circlize)  # for colorRamp2
p=Heatmap( x, name = "Proportion", row_names_side = "right",  # row labels on right
  column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10),
  col = colorRamp2(c(0, max(x)), c("white", "red")), row_names_max_width = unit(13, "cm") 
)
svg("figures/heatmap_cluster_vs_celltype_azimuth.svg", width = 11, height = 7)
draw(p)
dev.off()


markers=unique(mk[["moumita"]]$feature)
p <- DotPlot( d, features = markers, group.by = "azimuth_fine") +
  scale_color_viridis_c() + theme_classic() + theme( axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
ggsave(p,file="figures/dotplot_moumita_markers.svg",width=21,height=7)

d=readRDS(input)
mk=list();
mk[["barry"]]=fread("data/barry2020_long.txt")
mk[["moumita"]]=fread("data/moumita2022_long.txt")
mk[["moumita"]]$feature <- toupper(mk[["moumita"]]$feature)






markers <- mk[["moumita"]]
feature_list <- split(markers$feature, markers$type)

p <- DotPlot( d, features = split(mk[["moumita"]]$feature, mk[["moumita"]]$type), group.by = "seurat_clusters") +
  scale_color_viridis_c() + theme_classic() + theme( axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
ggsave(p,file="figures/dotplot_moumita_markers.svg",width=21,height=7)

barry <- mk[["barry"]]
barry <- barry[!duplicated(barry$feature)]
feature_list <- split(barry$feature, barry$type)
feature_list <- feature_list[order(names(feature_list))]
p <- DotPlot( d, features = feature_list, group.by = "seurat_clusters") +
  scale_color_viridis_c() + theme_classic() + theme( axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
ggsave(p,file="figures/dotplot_barry_markers.svg",width=28,height=7)


source("lib.r")
for( x in names(mk)){
    d=add_modules_named(d, mk[[x]],x)
}

output="bigdata/seurat5/s13_celltyped.rds"

library(SingleR)
library(celldex)

ref <- celldex::HumanPrimaryCellAtlasData()
expr <- GetAssayData(d, layer = "data")  # instead of slot="data"
pred <- SingleR( test = expr, ref = ref, labels = ref$label.main)
d$SingleR_label <- pred$labels


options(future.globals.maxSize = 20 * 1024^3)  # 20 GB
hlca.ref=readRDS("bigdata/azimuth/human_lung.Rds")

DefaultAssay(d) <- "RNA"
d <- SCTransform(d, verbose = TRUE)  # converts query to SCT

anchors <- FindTransferAnchors(
  reference = hlca_ref,
  query = d,
  reference.assay = "refAssay",    # SCT assay in reference
  query.assay = "RNA",             # Log-normalized
  normalization.method = "SCT",    # Use SCT mode to handle reference
  reduction = "pcaproject",
  dims = 1:50
)

anchors <- FindTransferAnchors(
  reference = hlca_ref, query = d,
  reference.assay = "refAssay",   # use the assay used in the reference
  query.assay = "RNA",
  reduction = "pcaproject",
  dims = 1:50
)



hlca.ref <- NormalizeData(hlca.ref)
hlca.ref <- FindVariableFeatures(hlca.ref)
hlca.ref <- ScaleData(hlca.ref)
hlca.ref <- RunPCA(hlca.ref)


DefaultAssay(hlca.ref) <- "SCT"
DefaultAssay(hlca.ref) <- "SCT"

# Recompute variable features
hlca.ref <- FindVariableFeatures( hlca.ref, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

# Recompute PCA
hlca.ref <- RunPCA( hlca.ref, features = VariableFeatures(hlca.ref), verbose = FALSE)

hlca.ref <- FindVariableFeatures( hlca.ref, selection.method = "vst", nfeatures = 3000)
anchors <- FindTransferAnchors( reference = hlca.ref, query = d, dims = 1:30,  normalization.method = "SCT")

predictions <- TransferData( anchorset = anchors, refdata = hlca.ref$cell_type, dims = 1:30)

d$HLCA_label <- predictions$predicted.id



#
# manual cell identification 
#

j<- grep("^moumita_", colnames(d@meta.data), value = TRUE) 
svg("figures/dotplot_moumita_marker_s13_clustered.svg",width=7,height=10)
DotPlot(d, features = j) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off();

## use cell similarity 
avg_mat <- aggregate( d@meta.data[, j], by = list(cluster = Idents(d)), FUN = mean)
rownames(avg_mat) <- avg_mat$cluster
avg_mat$cluster <- NULL
sim_mat <- cor(t(avg_mat))
#sim_mat 

d@meta.data$celltype="unknown"
d$celltype[ d$seurat_clusters %in% c(0,6,8,10,11,13,14,15,16,24) ] = "Basal-cells";
d$celltype[ d$seurat_clusters %in% c(1,5,19,20,25,29) ] = "Fibroblasts";
d$celltype[ d$seurat_clusters %in% c(2,21,22) ] = "B-cells";
d$celltype[ d$seurat_clusters %in% c(3,27) ] = "T-cells";
d$celltype[ d$seurat_clusters %in% c(4,7,12,28) ] = "Endothelial-cells";
d$celltype[ d$seurat_clusters %in% c(9) ] = "Ciliated-I";
d$celltype[ d$seurat_clusters %in% c(18) ] = "Ciliated-II";
d$celltype[ d$seurat_clusters %in% c(14) ] = "Club-I";
d$celltype[ d$seurat_clusters %in% c(17) ] = "Macrophages";
d$celltype[ d$seurat_clusters %in% c(19) ] = "Club-II";
d$celltype[ d$seurat_clusters %in% c(23) ] = "Goblet";
d$celltype[ d$seurat_clusters %in% c(26) ] = "NK-cells";
fwrite(d@meta.data,file="data/s13-meta.data.csv",row.names=T)

# barry annotation
j<- grep("^barry", colnames(d@meta.data), value = TRUE) 
svg("figures/dotplot_barry_marker_clustered.svg",width=7,height=10)
DotPlot(d, features = j) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off();

## sviolin
library(Seurat)
library(ggplot2)
library(cowplot)

a <- VlnPlot(d, mk[[1]]$feature, stack = TRUE, sort = TRUE) + theme(legend.position = "none") + ggtitle("Identity on y-axis")



d=readRDS("s13-harmony.rds")

#
# read inventory
#
input = read_sheet("https://docs.google.com/spreadsheets/d/1Rc5sBJCwspZsu9GYABqwmR_gQmumvzD0QFQut5TwDQE/edit#gid=0") %>%
        mutate( cellranger_dir=paste0("../cellranger/",cellranger_dir))


#
# read markers
#
urls="
https://docs.google.com/spreadsheets/d/10lqT7BiU2sitaUYLTP9Sg24l7m5HrqEy9x8hQCyXC1A/edit#gid=0
https://docs.google.com/spreadsheets/d/1wiR9ktw7VdsIP83vWlpbAR-J29tsrpC6JRuGrqyscYc/edit#gid=1034243033
"
urls=strsplit(trimws(urls),"\\s+",perl=T)[[1]]
library(googlesheets4)

# mk[[1]] : barry marker PMID: 32692579 
# mk[[2]] : ghosh markers
mk = lapply( urls, read_sheet )
names(mk) = c("barry","moumita")


mk[[2]] = mk[[2]][ mk[[2]]$type != "TARGET",]


#
# assign clusters to known markers
# 
#d=readRDS("s13-harmony.rds")

d[["H0"]] = input[match(d$orig.ident, input$id),]$H0
d[["H1"]] = input[match(d$orig.ident, input$id),]$H1
features=c("KRT5","MKI67","FOXJ1","SCGB1A1","MUC5B","PTPRC");
FeaturePlot(d, features=mk,ncol=6,by.col=F,split.by="H0")
FeaturePlot(d, features=mk,ncol=6,split.by="H1")
FeaturePlot(d, features=mk,ncol=7,by.col=F,split.by="orig.ident")

DimPlot(d,label=T,group.by="H0",split.by="H0")
DimPlot(d,group.by="H0",cols=c("green","red"))
DimPlot(d,group.by="H1",cols=c("purple","red","green"))
DimPlot(subset(d, subset= H0=="green"),ncol=7, split.by="orig.ident")
DimPlot(subset(d, subset= H0=="red"),ncol=7, split.by="orig.ident")
DimPlot(subset(d, subset= H0=="red"),ncol=7, split.by="H0")
DimPlot(subset(d, subset= H0=="green"),ncol=7, split.by="H0")



library(tibble)
library(Matrix.utils)
d1 =  aggregate.Matrix( t( d$RNA@data ), groupings = d@meta.data$seurat_clusters, fun = "sum") %>%
        t() %>% data.frame(check.names=F) %>% rownames_to_column(var="feature")
i=apply(d1[,-1], 1, var) > 0
d1=d1[i,]



#tmp=merge(mk[[1]],mk[[2]],all=T,by="feature") 
tmp=mk[[1]]
d2 = merge(d1,tmp,by="feature") %>% distinct( feature, .keep_all=T)
m=d2[,2:ncol(d1)];
row.names(m)=d2$feature


library(ComplexHeatmap)

y1=data.frame(rbind(table(Idents(d),d$H0)))
#y1 = y1/apply(y1,1,sum)

y2=data.frame(rbind(table(Idents(d),d$H1)))
#y2 = y2/apply(y2,1,sum)
ha=do.call(HeatmapAnnotation, apply(cbind(y1,y2),2,function(x) anno_barplot(x)))
#ha=HeatmapAnnotation( green_vs_red= anno_barplot( log2( y[,1]/sum(y[,1])) -log2( y[,2]/sum(y[,2])))) 

Heatmap(t(scale(t(m))), 
	top_annotation=ha, row_title_rot=0, column_names_rot=0,
	row_split=d2$type,  
	row_names_gp = gpar(fontsize = 7));

d2 = merge(d1,y,by="feature") %>% arrange( p_val)


