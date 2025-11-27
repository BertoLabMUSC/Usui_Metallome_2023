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
library(RRHO2)
})

# Gene Ontology
dge <- read.table("dge/Metallome_Dge_MouseID.txt",header=T,sep="\t")

all <- dge %>%
        filter(Direction == "All")


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_all.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(GeneOnto)
dev.off()

pdf("dge/GeneOnto_all_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(GeneOnto, categorySize="pvalue")
dev.off()

openxlsx::write.xlsx(as.data.frame(GeneOnto), 
                     file = "dge/GeneOnto_all.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

# UpReg
all <- dge %>%
        filter(Direction == "Upreg")


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

pdf("dge/GeneOnto_UpReg.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(GeneOnto)
dev.off()

pdf("dge/GeneOnto_UpReg_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(GeneOnto, categorySize="pvalue")
dev.off()

openxlsx::write.xlsx(as.data.frame(GeneOnto), 
                     file = "dge/GeneOnto_UpReg.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


# DownReg
all <- dge %>%
        filter(Direction == "Downreg")


GOI <- bitr(as.character(all$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

GeneOnto <- enrichGO(gene = unique(GOI$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

pdf("dge/GeneOnto_DownReg.pdf",width=6,height=5,useDingbats=FALSE)
dotplot(GeneOnto)
dev.off()

pdf("dge/GeneOnto_DownReg_Network.pdf",width=8,height=8,useDingbats=FALSE)
cnetplot(GeneOnto, categorySize="pvalue")
dev.off()

openxlsx::write.xlsx(as.data.frame(GeneOnto), 
                     file = "dge/GeneOnto_DownReg.xlsx", 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)