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

d2 <- read.table("dge/Both_Dge_MouseID.txt",header=T)

d2 <- d2 %>% filter(Direction == "Both")

GOI <- bitr(as.character(d2$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_TETA_Both_BP.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "dge/GeneOnto_TETA_Both_BP.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)



D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "MF", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_TETA_Both_MF.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "dge/GeneOnto_TETA_Both_MF.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)



D2_GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "CC", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_TETA_Both_CC.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(D2_GeneOnto)
dev.off()

openxlsx::write.xlsx(as.data.frame(D2_GeneOnto), 
                     file = "dge/GeneOnto_TETA_Both_CC.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)
