setwd("~/bigdata/Slycopersicum/analysis/")
#Install all the things required for Sleuth
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

library(sleuth)
library(dplyr)
library(ggplot2)
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
full_design <- model.matrix(formula(~ ns(day, df = 4)))

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

#Who's on top?
head(lrt_results, n=20)[,1:3]

#plot a heatmap
plot_transcript_heatmap(so, head(lrt_results, n = 20)$target_id, 'est_counts')

#plot a single DE transcript
tmp <- so$obs_raw %>% dplyr::filter(target_id == 'Solyc06g069430.3.1') #FUL1 (FUL2=Solyc03g114830.3.1)
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp #its expression
#make actual graph
tmp <- transform(tmp, day = as.numeric(day))
ggplot(tmp, aes(x=day, y=est_counts)) + geom_point(shape=1) + geom_smooth(method = loess)

#try shiny
sleuth_live(so)
