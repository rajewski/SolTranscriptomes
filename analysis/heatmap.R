#simple script to offload heatmap construction onto a batch job
library(RColorBrewer)
library(gplots)
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)
setwd("~/bigdata/Slycopersicum/analysis/")
mats <- readRDS("mats.RDS")
clust.div <- readRDS("clust.div.RDS")


pdf("DEGs_filtered_pam15heatmap.pdf")
heatmap.2(mats,
          main = "Tomato Pericarp RNA-Seq",
          notecol="black",
          col=my_palette,            # this sets the color pallet made above
          Colv=FALSE,
          Rowv=FALSE,
          labRow=F,    
          cexRow=1,
          #labCol=fignames[col.use],                  # override column names for figure
          dendrogram="none",
          keysize=1,                                            #key fontsize
          cexCol=1,                                             #changes font size of the column labels
          trace="none",                                         # removes trace lines from heatmap
          density.info="none",                             # "histogram", "density" or "none" for the key
          #colsep=c(3,6,11,15,18,21,24,29),                     # column indices of vertical white lines 
          rowsep=clust.div,                                     # row indices of horizontal white lines
          sepcolor="black",
          sepwidth=c(0.001,0.001),                                  # sets thickness of white lines
          margins=c(10,10),                                       # adjust space for 1-column labels, 2-row labels 
          #lwid=c(0.1,1),                                       # 2x2 relat. widths - Left-key & Rdendo, Right-Cdendo & heatmap
          #lhei=c(0.1,1),                                       # 2x2 relat. heights - Top-key & Cdendo, Bott.-Rdendo & heatmap
          breaks=seq(-3,3,by=6/length(my_palette))             # heat range and bin size - must match no. of colors
)
dev.off()

