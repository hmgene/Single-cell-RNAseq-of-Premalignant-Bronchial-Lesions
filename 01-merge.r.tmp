library(data.table)
input=fread("data/sample_input.csv")
input[,cellranger_dir := paste0("bigdata/cellranger/",cellranger_dir)]
input_info="data/sample_input.csv"
output="bigdata/seurat5/s13-merged.rds"

m=fread(input_info)


d.list = list()
library(Seurat)
for ( i in 1:nrow(input)){
    id=input[i,]$id
    f= input[i,]$cellranger_dir
    if(id !="" && file.exists(f)){
        counts <- Read10X(data.dir = f)
        d.list[[id]] <- CreateSeuratObject(counts = counts, project = id)
        #d.list[[id]]= Read10X( data.dir = f) %>% CreateSeuratObject(project = id)
    }else{
                warning(paste0(f," doesnot exist!"));
    }
}

 
d = merge(d.list[[1]],d.list[2:length(d.list)])
d$selfrenewal <- m[match(d$orig.ident, m$id), H2]
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")
d[["percent.rb"]] <- PercentageFeatureSet(d, pattern = "^RP[SL]")

#> dim(d)
#[1] 45068 73773
# table(d$orig.oident)
#   s1   s10   s11   s12   s13    s2    s3    s4    s5    s6    s7    s8    s9 
#  495  6358 13170  2152  2761  1065  1439  2742  1793  1233 18960 10290 11315 


library(viridis)
library(ggplot2)
library(dplyr)
df <- d@meta.data
df$mt_status <- ifelse(df$percent.mt > 25, "MT > 25%", "MT ≤ 25%")

# MT
p=ggplot(df, aes(x = orig.ident, y = percent.mt, fill = selfrenewal)) +
  geom_violin(alpha = 0.7, color = NA) +
  geom_hline(yintercept = 25, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = c( "high_selfrenewal" = "#2C7BB6", "low_selfrenewal"  = "#D7191C")) +
  theme_classic() + theme( axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(p,file="figures/qc_mt.svg",width=7,height=4)

# nCount vs nFeature

df$qc_mt <- ifelse(df$percent.mt < 25, "MT < 25%", "MT ≥ 25%")
p <- ggplot(df) +
  geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = qc_mt), size = 0.6, alpha = 0.35) +
  geom_vline(xintercept = 500, linetype = "dotted", linewidth = 0.8, color = "black") +
  geom_hline(yintercept = 200, linetype = "dotted", linewidth = 0.8, color = "black") +
  scale_color_manual( name = "Mitochondrial %", values = c("MT < 25%" = "#2C7BB6", "MT ≥ 25%" = "#E15759")) +
  scale_x_log10() + scale_y_log10() + facet_grid(. ~ orig.ident) +
  theme_classic() + theme( legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold")) +
  labs(x = "UMI counts (log10)", y = "Detected genes (log10)")
ggsave("figures/qc_ncount_vs_nfeature.svg", plot = p, width = 14, height = 4)

saveRDS(d,file=output)


