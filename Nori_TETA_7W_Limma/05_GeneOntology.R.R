suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(UpSetR)
library(RColorBrewer)
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(clusterProfiler)
})

# Gene Onto

d2 <- read.table("dge/TETA_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "TETA_All")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_TETA_All.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "dge/GeneOnto_TETA_All.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

# Split by signature
d2 <- read.table("dge/TETA_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "TETA_Upreg")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_TETA_Upreg.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "dge/GeneOnto_TETA_Upreg.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)



d2 <- read.table("dge/TETA_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "TETA_Downreg")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_TETA_Downreg.pdf",width=8,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "dge/GeneOnto_TETA_Downreg.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)
