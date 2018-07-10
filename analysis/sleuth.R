setwd("~/bigdata/Slycopersicum/analysis/")
#Install all the things required for Sleuth
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

library(sleuth)
library(dplyr)
library(ggplot2)
library(cowplot) #for pretty graphs
#make experiemental information matrix from scratch rather than an import
s2c <- read.table("../replist.txt", header=F, stringsAsFactors = F)
colnames(s2c) <- "sample"
s2c$rep <- rep(c(1,2,3),5)
s2c$stage <- substr(s2c[,1],1,nchar(s2c[,1])-1)
s2c$apprDPA <- c(1,1,1,3,3,3,15,15,15,32,32,32,40,40,40) #only for spline model
resultspath <- "~/bigdata/Slycopersicum/Results_"
s2c$path <- paste(resultspath,s2c$sample,sep="")
s2c

#make a design matrix for spline based modeling
#admittedly Im not quite certain how this would differ from just using ~day as the model
#I am going to try to do both  and see if it differs
library(splines)
day <- s2c$apprDPA
full_design <- model.matrix(formula(~ ns(day, df = 3)))

#construct sleuth object or read it if you have already done this
#so <- readRDS("Sleuth_so.rds")
so <- sleuth_prep(s2c, full_model = full_design, max_bootstrap = 100, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so)
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")
saveRDS(so, "Sleuth_so.rds") #save this so we dont have to repeat it

#lets do some plots
plot_qq(so, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)

#extract LRT test results
lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
#"true" means these transcripts' expression is better explained by a spline than noise
table(lrt_results[,"qval"] < 0.05) 
head(lrt_results[lrt_results$qval<0.05,])
#Who's on top?
head(lrt_results, n=20)[,1:3]

#plot a heatmap
plot_transcript_heatmap(so, head(lrt_results[lrt_results$qval<0.05,], n=100)$target_id, cluster_transcripts = TRUE, 'est_counts')


#plot a single DE transcript
tmp <- so$obs_raw %>% dplyr::filter(target_id == 'Solyc06g069430.3.1') #FUL1 (FUL2=Solyc03g114830.3.1)
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp #its expression
#make actual graph
tmp <- transform(tmp, day = as.numeric(day))
ggplot(tmp, aes(x=day, y=est_counts)) + geom_point(shape=1) + geom_smooth(method = loess)

#try shiny
sleuth_live(so)


###You must run the glimma script before you can do this analysis. It depends on files made by that script
#try to make a plot of a cluster of genes
#get a list of expression for each gene
library(tidyr)
data <- so$obs_raw[,c(1,2,6)] %>% spread(sample,est_counts)
rownames(data) <- data[,1] #give row names
data <- data[,-1] #drop names
data <- data[,c(4:9,1:3,10:15)] #reorder columns by time

#subset this by cluster
clust1 <- read.delim("DEGs_filtered15_Clust1GENE",header=F)
clust1 <- data[clust1$V1,]
tmp <- so$obs_raw[so$obs_raw$target_id %in% clust1$V1,] 
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp #its expression
tmp <- tmp[,c("apprDPA","target_id","est_counts")]
tmp <- aggregate(est_counts~apprDPA+target_id,tmp,mean)
#make actual graph
ggplot(tmp, aes(x=apprDPA, y=est_counts, group=target_id)) + geom_line()

#try to get them all at once
clusters <- read.delim("DEGs_filtered_pam15reordered.txt")
clusters <- clusters[,c("X","clust.sort")]
colnames(clusters) <- c("target_id","clust.sort")
expr <- dplyr::full_join(so$sample_to_covariates, so$obs_raw, by = 'sample')
expr <- expr[,c("apprDPA","target_id","est_counts")]
expr <- aggregate(est_counts~apprDPA+target_id,expr,mean)
expr <- dplyr::full_join(clusters,expr, by = "target_id")
exprsumm <- aggregate(est_counts~apprDPA+clust.sort, expr, mean)

#Make Graph itself
pdf(file = "ClusterProfiles.pdf", width=14, height=10)
ggplot(expr %>% drop_na(), aes(x=apprDPA, y=est_counts, group=target_id)) + 
  ylab("Estimated Counts") +
  xlab("Days Post-Anthesis") +
  geom_line() +
  facet_wrap(~clust.sort, ncol=5, scales="free")
dev.off()

pdf(file="ClusterProfilesBasic.pdf", width=14, height=10)
ggplot(data=exprsumm %>% drop_na(), aes(x=apprDPA, y=est_counts)) + 
  ylab("Estimated Counts") +
  xlab("Days Post-Anthesis") +
  ggtitle("Clustered Expression Profiles") +
  geom_line() +
  facet_wrap(~clust.sort, ncol=5, scales="free")
dev.off()




