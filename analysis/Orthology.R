#compare nobt input with nobt output from diamond through Slyc
NobtUptoSlyc <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/NobtUptoSlycHits.tsv")
NobtUptoSlyc$V2 <- gsub(".{2}$", "", NobtUptoSlyc$V2)
NobtUptoNobt <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/NobtUptoNobtHits.tsv")
NobtUptoNobt$V1 <- gsub(".{2}$", "", NobtUptoNobt$V1)


#merge the tables of the Nobt-Slyc and Slyc-Nobt hits
NobtUp <- merge(NobtUptoSlyc, NobtUptoNobt, by.x="V2", by.y="V1")
#compare if the Nobt genes match
NobtUp$Ortho <- as.character(NobtUp$V1) == as.character(NobtUp$V2.y)
#Flag any Slyc gene where 1) input and output differ or 2) there are multiple matches
#aggregate the NobtUp set by the Slyc gene they map to and flag them as false (0) if any Slyc gene is a "bridge" between a mismatched input and output gene
NobtUpAggr <- aggregate(NobtUp$Ortho~NobtUp$V2, FUN=min)
#filter this aggregate to remove Slyc genes that are bad bridges
NobtUpAggr <- NobtUpAggr[NobtUpAggr$`NobtUp$Ortho`==0,]
#filter the orthoset to remove the multiple Slyc genes
NobtUpFilt <- NobtUp[!(NobtUp$V2 %in% NobtUpAggr$`NobtUp$V2`),]
#Final Set of Slyc Orthologs for Nobt DEGs
NobtUpFilt <- NobtUpFilt[,c(2,1)]
colnames(NobtUpFilt) <- c("Nobt Hit", "Slyc Ortholog")

#repeat for Nobt Down
NobtDowntoSlyc <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/NobtDowntoSlycHits.tsv")
NobtDowntoSlyc$V2 <- gsub(".{2}$", "", NobtDowntoSlyc$V2)
NobtDowntoNobt <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/NobtDowntoNobtHits.tsv")
NobtDowntoNobt$V1 <- gsub(".{2}$", "", NobtDowntoNobt$V1)
NobtDown <- merge(NobtDowntoSlyc, NobtDowntoNobt, by.x="V2", by.y="V1")
NobtDown$Ortho <- as.character(NobtDown$V1) == as.character(NobtDown$V2.y)
NobtDownAggr <- aggregate(NobtDown$Ortho~NobtDown$V2, FUN=min)
NobtDownAggr <- NobtDownAggr[NobtDownAggr$`NobtDown$Ortho`==0,]
NobtDownFilt <- NobtDown[!(NobtDown$V2 %in% NobtDownAggr$`NobtDown$V2`),]
NobtDownFilt <- NobtDownFilt[,c(2,1)]
colnames(NobtDownFilt) <- c("Nobt Hit", "Slyc Ortholog")

#repeat for Slyc up
SlycUptoNobt <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/SlycUptoNobtHits.tsv")
SlycUptoSlyc <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/SlycUptoSlycHits.tsv")
SlycUptoSlyc$V2 <- gsub(".{2}$", "", SlycUptoSlyc$V2)
SlycUp <- merge(SlycUptoNobt, SlycUptoSlyc, by.x="V2", by.y="V1")
SlycUp$Ortho <- as.character(SlycUp$V1) == as.character(SlycUp$V2.y)
SlycUpAggr <- aggregate(SlycUp$Ortho~SlycUp$V2, FUN=min)
SlycUpAggr <- SlycUpAggr[SlycUpAggr$`SlycUp$Ortho`==0,]
SlycUpFilt <- SlycUp[!(SlycUp$V2 %in% SlycUpAggr$`SlycUp$V2`),]
SlycUpFilt <- SlycUpFilt[,c(2,1)]
colnames(SlycUpFilt) <- c("Slyc Hit", "Nobt Ortholog")

#repeat for Slyc Down
SlycDowntoNobt <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/SlycDowntoNobtHits.tsv")
SlycDowntoSlyc <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/SlycDowntoSlycHits.tsv")
SlycDowntoSlyc$V2 <- gsub(".{2}$", "", SlycDowntoSlyc$V2)
SlycDown <- merge(SlycDowntoNobt, SlycDowntoSlyc, by.x="V2", by.y="V1")
SlycDown$Ortho <- as.character(SlycDown$V1) == as.character(SlycDown$V2.y)
SlycDownAggr <- aggregate(SlycDown$Ortho~SlycDown$V2, FUN=min)
SlycDownAggr <- SlycDownAggr[SlycDownAggr$`SlycDown$Ortho`==0,]
SlycDownFilt <- SlycDown[!(SlycDown$V2 %in% SlycDownAggr$`SlycDown$V2`),]
SlycDownFilt <- SlycDownFilt[,c(2,1)]
colnames(SlycDownFilt) <- c("Slyc Hit", "Nobt Ortholog")


#Now begin to compare the Up lists for both species to see which genes are shared
SlycUpDEG <- read.csv("~/bigdata/Slycopersicum/slyc-WT/Slyc_1v15_DESeq_Up.csv")
NobtUpFilt$InSlycUp <- NobtUpFilt$`Slyc Ortholog` %in% SlycUpDEG$X
NobtUpDEG <- read.csv("~/bigdata/Slycopersicum/slyc-WT/Nobt_Prev16_DESeq_Up.csv")
SlycUpFilt$InNobtUp <- SlycUpFilt$`Nobt Ortholog` %in% NobtUpDEG$X
ConservedUp <- merge(SlycUpFilt[SlycUpFilt$InNobtUp==TRUE,], NobtUpFilt[NobtUpFilt$InSlycUp==TRUE,],
      by.x="Slyc Hit", by.y="Slyc Ortholog")
ConservedUp <- ConservedUp[,c(1,4)] #should do a go enrichment with this set of genes
write.csv(ConservedUp, file="ConservedUp.csv",quote = FALSE, row.names = FALSE)

#compare the down-regulated genes for both species
SlycDownDEG <- read.csv("~/bigdata/Slycopersicum/slyc-WT/Slyc_1v15_DESeq_Down.csv")
NobtDownFilt$InSlycDown <- NobtDownFilt$`Slyc Ortholog` %in% SlycDownDEG$X
NobtDownDEG <- read.csv("~/bigdata/Slycopersicum/slyc-WT/Nobt_Prev6_DESeq_Down.csv")
SlycDownFilt$InNobtDown <- SlycDownFilt$`Nobt Ortholog` %in% NobtDownDEG$X
ConservedDown <- merge(SlycDownFilt[SlycDownFilt$InNobtDown==TRUE,], NobtDownFilt[NobtDownFilt$InSlycDown==TRUE,],
                     by.x="Slyc Hit", by.y="Slyc Ortholog")
library(dplyr) #idk why but I'm getting duplicate rows
ConservedDown <- distinct(ConservedDown)
ConservedDown <- ConservedDown[,c(1,4)]
write.csv(ConservedDown, file="ConservedDown.csv", quote=FALSE, row.names=FALSE)



