# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(future.apply)
library(DESeq2)
library(pheatmap)
library(sva)
})


dir.create("dge")

load("futcounts/Expression_Input.RData")

pd <- data.frame(row.names=colnames(exp), Genotype = as.factor(c(rep("CTL",5),rep("TETA",4))))

# Filter the expression by condition
filter=apply(rpkm, 1, function(x) (all(x[1:5] >= 0.5) | all(x[6:9] >= 0.5)))
count_filt <- exp[filter,]
rpkm_filt <- rpkm[filter,]

logCPM <- log2(rpkm_filt+1)
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)

pdf("dge/PCA.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette=c("black","red"), 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() 
dev.off()

# DGE 
design <- model.matrix(~.,pd) # describe model to be fit

fitLM = lmFit(p,design)#,method="robust");
fitEb = eBayes(fitLM)##robust = TRUE);

TETA7W_FullTab = topTable(fitEb, coef = "GenotypeTETA",number=nrow(p)) %>%
             rownames_to_column("Gene")


TETA7W_DGE <- TETA7W_FullTab %>%
                  mutate(Abs = abs(logFC)) %>%
                  filter(adj.P.Val < 0.05) %>%
                  arrange(desc(Abs))

save(TETA7W_FullTab,TETA7W_DGE,pd, file = "dge/TETA_Dge_Data.RData")


openxlsx::write.xlsx(TETA7W_FullTab, 
                     file = "dge/TETA_DGE_FullStats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(TETA7W_DGE, 
                     file = "dge/TETA_DGE.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


# Input for viz Plot
df <- TETA7W_FullTab %>% 
        mutate(LOG = -log10(adj.P.Val), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(adj.P.Val < 0.05, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0 & adj.P.Val < 0.05 ~ "UpReg", logFC < 0 & adj.P.Val < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(adj.P.Val) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
mat <- p[rownames(p)%in% top_labelled$Gene,] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("dge/Boxplots_TopGenes_TETA.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("black", "red")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/Vulcano_Plot_TETA.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-4,+4)
dev.off()

# heatmap
mat <- p[rownames(p)%in% TETA7W_DGE$Gene,]
anno <- pd
Genotype        <- c("red", "black")
names(Genotype) <- c("TETA", "CTL")
anno_colors <- list(Genotype = Genotype)
pdf("dge/Heatmap_DGE.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()










