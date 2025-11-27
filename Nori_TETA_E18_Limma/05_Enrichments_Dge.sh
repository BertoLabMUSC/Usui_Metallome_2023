# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir dge/enrichments_ctx/

cp utils/geneset/*.RData dge/enrichments_ctx/
cp utils/Enrichment.r dge/enrichments_ctx/
cp dge/Metallome_Dge_HumanID.txt dge/enrichments_ctx/
cp dge/Metallome_Dge_MouseID.txt dge/enrichments_ctx/

cd dge/enrichments_ctx/

mkdir STATS/

# scRNA healty
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 15585 -o STATS/Allen_Markers -W 4 -H 4
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l GeneSets_BBLake_2018.RData -p -b 15585 -o STATS/BBLake_Markers -W 4 -H 4
Rscript Enrichment.r -g Metallome_Dge_MouseID.txt -l GeneSets_scMouse.RData -p -b 15585 -o STATS/scMouse -W 4 -H 5
Rscript Enrichment.r -g Metallome_Dge_MouseID.txt -l GeneSets_Cowan_PFC.RData -p -b 19982 -o STATS/scMouse_Cowan_PFC -W 4 -H 5
Rscript Enrichment.r -g Metallome_Dge_MouseID.txt -l GeneSets_AdultMouse.RData -p -b 19982 -o STATS/scMouse_AdultMouse_Brain -W 4 -H 6
Rscript Enrichment.r -g Metallome_Dge_MouseID.txt -l GeneSets_YoungMouse.RData -p -b 19982 -o STATS/scMouse_YoungMouse_Brain -W 4 -H 6
Rscript Enrichment.r -g Metallome_Dge_MouseID.txt -l GeneSets_OligoMouse.RData -p -b 19982 -o STATS/scMouse_OligoMouse_Brain -W 4 -H 4

# scRNA Disorders
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell -W 4 -H 5
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell -W 4 -H 5

# Neuropsy
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 4 -H 3
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 4 -H 5
Rscript Enrichment.r -g Metallome_Dge_HumanID.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari -W 4 -H 2

rm *.RData

