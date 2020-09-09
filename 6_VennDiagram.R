# Orthogroup Venn Diagram -------------------------------------------------
library("VennDiagram")
library("dplyr")
source("X_Functions.R")

Orthos <- read.table("Orthofinder/OrthoFinder/Results_Aug31/Orthogroups/Orthogroups.tsv",
                     sep="\t",
                     stringsAsFactors = F,
                     header=T)
#All Orthogroups
# Remove empties for each species
Set1 <- Orthos[,c(1:2)] %>% filter_all(all_vars(!grepl("^$",.)))
Set2 <- Orthos[,c(1,3)] %>% filter_all(all_vars(!grepl("^$",.)))
Set3 <- Orthos[,c(1,4)] %>% filter_all(all_vars(!grepl("^$",.)))
Set4 <- Orthos[,c(1,5)] %>% filter_all(all_vars(!grepl("^$",.)))
venn.diagram(
  x = list(Set1[,1], Set2[,1], Set3[,1], Set4[,1]),
  category.names = c("Arabidopsis","Cucumis", "Nicotiana", "Solanum"),
  filename = "Figures/AllOrthogroups.png",
  imagetype = "png",
  main="Orthogenes",
  main.fontfamily = "sans",
  main.fontface = "bold",
  cat.fontface="italic",
  cat.fontfamily="sans",
  fontfamily="sans",
  scaled=T,
  euler.d=T,
  output=F,
  lwd=2,
  lty=1)

# All Species, single copy
#SingleOrthos <- Orthos %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
SingleSet1 <- Set1 %>% filter_all(all_vars(!grepl(',',.)))
SingleSet2 <- Set2 %>% filter_all(all_vars(!grepl(',',.)))
SingleSet3 <- Set3 %>% filter_all(all_vars(!grepl(',',.)))
SingleSet4 <- Set4 %>% filter_all(all_vars(!grepl(',',.)))
venn.diagram(
  x = list(SingleSet1[,1], SingleSet2[,1], SingleSet3[,1], SingleSet4[,1]),
  category.names = c("Arabidopsis","Cucumis", "Nicotiana", "Solanum"),
  filename = "Figures/SingleCopyOrthogroups.png",
  imagetype = "png",
  main="Shared Single-Copy Orthogenes",
  main.fontfamily = "sans",
  main.fontface = "bold",
  cat.fontface="italic",
  cat.fontfamily="sans",
  fontfamily="sans",
  scaled=T,
  euler.d=T,
  output=F,
  lwd=2,
  lty=1)


