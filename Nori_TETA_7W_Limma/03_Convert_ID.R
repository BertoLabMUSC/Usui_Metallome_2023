suppressPackageStartupMessages({
library(biomaRt)
library(tidyverse)
})

# Load Data
load("dge/TETA_Dge_Data.RData")

# Convert to human 
human = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")


MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = TETA7W_DGE$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(TETA7W_DGE,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(logFC > 0 ~ "TETA_Upreg", logFC < 0  ~ "TETA_Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("TETA_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene)) %>%
        arrange(Direction)

write.table(df,"dge/TETA_Dge_HumanID.txt",sep="\t",quote=F,row.names=F)


df <- TETA7W_DGE %>%
                mutate(Direction = case_when(logFC > 0 ~ "TETA_Upreg", logFC < 0  ~ "TETA_Downreg")) %>%
                dplyr::select(Gene,Direction)

tmp <- data.frame(Gene = df$Gene, Direction = rep("TETA_All",nrow(df)))

df <- rbind(df,tmp) %>%
        mutate(Gene = as.character(Gene))%>%
        arrange(Direction)

write.table(df,"dge/TETA_Dge_MouseID.txt",sep="\t",quote=F,row.names=F)
