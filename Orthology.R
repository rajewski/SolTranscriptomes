#compare nobt input with nobt output from diamond through Slyc
NobtUptoSlyc <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/NobtUptoSlycHits.tsv")
NobtUptoNobt <- read.table("~/bigdata/Slycopersicum/slyc-WT/analysis/NobtUptoNobtHits.tsv")

#merge the tables of the Nobt-Slyc and Slyc-Nobt hits
NobtUp <- merge(NobtUptoSlyc, NobtUptoNobt, by.x="V2", by.y="V1")
#compare if the Nobt genes match
NobtUp$Ortho <- as.character(NobtUp$V1) == as.character(NobtUp$V2.y)
#Flag any Slyc gene where 1) input and output differ or 2) there are multiple matches
NobtUpAggr <- aggregate(NobtUp$Ortho~NobtUp$V2, FUN=min)
NobtUpAggr <- NobtUpAggr[NobtUpAggr$`NobtUp$Ortho`==0,]
#filter the orthoset to remove the multiple Slyc genes
NobtUpFilt <- NobtUp[!(NobtUp$V2 %in% NobtUpAggr$`NobtUp$V2`),]
#Final Set of Slyc Orthologs for Nobt DEGs
NobtUpFilt <- NobtUpFilt[,c(2,1)]
