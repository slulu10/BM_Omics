library(dplyr)

#prefix = "Seurat.MetsGBM.CombinedAnalysis.Apr2022.Lymphoid"
prefix = "Seurat.MetsGBM.CombinedAnalysis.Apr2022.Myeloid"

marker = "CD14"

outfile = paste(prefix, "CelltypeMarker.gmt", sep=".")


infile = "Seurat.MetsGBM.CombinedAnalysis.Apr2022.Myeloid.MinD0.1.sprd5.NN20.seed122.ep1000_CelltypeDiffGene.csv"
infile <- read.csv(file = infile)

celltypes <- levels(infile$cluster)

list <- c(1:length(celltypes))
for(i in 1:length(celltypes)){
  celltype = celltypes[i]
  subset = infile[infile$cluster==celltype&infile$p_val_adj<=0.05,]
  
  subset <- subset[order(-subset$DeltaPerLogFC),]
  
  genelist = as.character(subset$gene[1:100])
  
  line = do.call(paste, c(as.list(genelist), sep = "\t"))
  
  line = paste(celltype, prefix, marker, line, sep = "\t")
  
  list[i] <- line
  
}

writeLines(as.character(list), file(outfile), sep="\n")

