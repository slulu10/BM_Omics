require(ggplot2)
library(reshape)
library(ggpubr)
library(rgl)
library(car)
library(dplyr)
library(RColorBrewer)

prefix = "Seurat.MetsGBM.CombinedAnalysis.Apr2022.Lymphoid.CD8.DiffusionMapEmbedding"
infile = "Seurat.MetsGBM.CombinedAnalysis.Apr2022.Lymphoid.CD8.DiffusionMapEmbedding.txt"

df <- read.delim(infile, sep = "\t", header = T)
df$rownames = rownames(df)

df$celltype = factor(df$celltype, levels = c("CD8-IL7R-CD69-Tm",
                                                "CD8-ITGA1-ITGAL-Trm",
                                                "CD8-TCF7-CD226-Tprog.exh",
                                                "CD8-ISG High",
                                                "CD8-GZMK-CCL4-CCL5-Teff",
                                                "CD8-GZMH-GZMA-CD52-Teff",
                                                "CD8-VCAM1-IFNG-Tex",
                                                "CD8-CXCL13-LAG3-Tex"))


######assign the color if annotate the DM with celltypes##############
######################################################################
values = levels(df$celltype)
cols = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","maroon","#FDBF6F","#FF7F00")
#cols = c("#FFCDF3","#CC6677","#FDBF6F","#FF7F00","#FF1493","#C71585","#F0E68C")
#cols = c( "#88CCEE","#44AA99","#117733")

df <- df %>%
  mutate(Col = case_when(
    celltype == values[1] ~ cols[1],
    celltype == values[2] ~ cols[2],
    celltype == values[3] ~ cols[3],
    celltype == values[4] ~ cols[4],
    celltype == values[5] ~ cols[5],
    celltype == values[6] ~ cols[6],
    celltype == values[7] ~ cols[7],
    celltype == values[8] ~ cols[8]
  ))

rownames(df) <- df$rownames

xval = df$DC1
yval = df$DC2
zval = df$DC3

plot3d(x=xval, y=yval, z=zval,xlab="DC1",ylab="DC2",zlab="DC3", bty="n", lit=FALSE, col = df$Col,size=4)
#rgl.viewpoint(userMatrix = pp$userMatrix, fov = pp$FOV, zoom = pp$zoom)
pp <- par3d(no.readonly=TRUE)

## Save the list to a text file
dput(pp, file="MDM.irisView.R", control = "all")

outfile = paste(prefix,"3D","celltype","png", sep=".")
rgl.snapshot(outfile, fmt='png') # Check your default directory using command getwd() for the location of the image.
dev.off()


### 3d/2d scatter plot of gene expression################################
#########################################################################
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
pp <- dget("MDM.irisView.R")

gene = "JAML"
expr = df[[gene]] 

cols <- myColorRamp(c("grey", "blue"), expr) 

plot3d(x=xval, y=yval, z=zval,xlab="DC1",ylab="DC2",zlab="DC3",size=4, col=cols, lit=TRUE)
rgl.viewpoint(userMatrix = pp$userMatrix, fov = pp$FOV, zoom = pp$zoom)

outfile = paste(prefix,"3D","Expr",gene,"png", sep=".")
rgl.snapshot(outfile, fmt='png') # Check your default directory using command getwd() for the location of the image.


outfile = paste(prefix,"2D","Expr",gene,"png", sep=".")
png(filename = outfile,width = 600, height = 600)
plot(x=xval, y=yval,xlab="DC1",ylab="DC2" ,col = cols, pch = 16, cex=0.4)
title(main = gene)
dev.off()

##assign the color if annotate the DM with conditions########################
#############################################################################

xlim = range(df$DC1)
ylim = range(df$DC2)
zlim = range(df$DC3)

values = c("GBM.rec","GBM.pembro","met","met.pembro")
cols = c(alpha("cyan",alpha=0.9),
         alpha("blue",alpha = 0.9),
         alpha("pink",alpha=0.9),
         alpha("red",alpha = 0.9))

for(i in 1:length(values)){
  df.sub = subset(df, df$condition %in% c(values[i]))
  
  df.sub <- df.sub %>%
    mutate(Col = case_when(
      condition == values[1] ~ cols[1],
      condition == values[2] ~ cols[2],
      condition == values[3] ~ cols[3],
      condition == values[4] ~ cols[4]    
    ))
  
  rownames(df.sub) <- df.sub$rownames
  
  xval = df.sub$DC1
  yval = df.sub$DC2
  zval = df.sub$DC3
  
  
  plot3d(x=xval, y=yval, z=zval,xlab="DC1",ylab="DC2",zlab="DC3", type="p", lit=FALSE, col = df.sub$Col,size=4,
         xlim=xlim,ylim=ylim,zlim=zlim)
  rgl.viewpoint(userMatrix = pp$userMatrix, fov = pp$FOV, zoom = pp$zoom)
  
  outfile = paste(prefix,"3D","condition",values[i],"png", sep=".")
  rgl.snapshot(outfile, fmt='png') # Check your default directory using command getwd() for the location of the image.
  
}


