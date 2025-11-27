suppressPackageStartupMessages({
library(biomaRt)
library(tidyverse)
})

# Load Data
load("dge/Metallome_Dge_Data.RData")

# Convert to human 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = Metallome_DGE$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(Metallome_DGE,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df,"dge/Metallome_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df <- Metallome_DGE %>%
                mutate(Direction = case_when(logFC > 0 ~ "Upreg", logFC < 0  ~ "Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp <- data.frame(Gene = df$Gene, Direction = rep("All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df,"dge/Metallome_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)

