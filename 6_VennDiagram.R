# Orthogroup Venn Diagram -------------------------------------------------
library("VennDiagram")
library("wesanderson")
Orthos <- read.table("Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                     sep="\t",
                     stringsAsFactors = F,
                     header=T)
#All Orthogroups
# Remove empties for each species
Set1 <- Orthos[,c(1:2)] %>% filter_all(all_vars(!grepl("^$",.)))
Set2 <- Orthos[,c(1,3)] %>% filter_all(all_vars(!grepl("^$",.)))
Set3 <- Orthos[,c(1,4)] %>% filter_all(all_vars(!grepl("^$",.)))
venn.diagram(
  x = list(Set1[,1], Set2[,1], Set3[,1]),
  category.names = c("Arabidopsis", "Nicotiana", "Solanum"),
  filename = "Figures/AllyOrthogroups.png",
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
  lty=1,
  fill=wes_palette("Zissou1", 3, "continuous"))

# All Species, single copy
SingleOrthos <- Orthos %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
SingleSet1 <- SingleOrthos[,c(1:2)] %>% filter_all(all_vars(!grepl("^$",.)))
SingleSet2 <- SingleOrthos[,c(1,3)] %>% filter_all(all_vars(!grepl("^$",.)))
SingleSet3 <- SingleOrthos[,c(1,4)] %>% filter_all(all_vars(!grepl("^$",.)))
venn.diagram(
  x = list(SingleSet1[,1], SingleSet2[,1], SingleSet3[,1]),
  category.names = c("Arabidopsis", "Nicotiana", "Solanum"),
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
  lty=1,
  fill=wes_palette("Zissou1", 3, "continuous"))

# Dry Species Only
DryOrthos <- SingleOrthos[,1:3] #For just Dry Species
DrySet1 <- DryOrthos[,c(1:2)] %>% filter_all(all_vars(!grepl("^$",.)))
DrySet2 <- DryOrthos[,c(1,3)] %>% filter_all(all_vars(!grepl("^$",.)))

venn.diagram(
  x = list(DrySet1[,1], DrySet2[,1]),
  category.names = c("  Arabidopsis", "Nicotiana  "),
  filename = "Figures/DryOrthogroups.png",
  imagetype = "png",
  main="Shared Single-Copy Orthogenes",
  main.fontfamily = "sans",
  main.fontface = "bold",
  sub="Dry-Fruited Species",
  sub.fontfamily = "sans",
  sub.fontface = "plain",
  cat.fontface="italic",
  cat.fontfamily="sans",
  fontfamily="sans",
  scaled=T,
  euler.d=T,
  output=F,
  lwd=2,
  lty=1,
  fill=wes_palette("Zissou1", 3, "continuous")[1:2])


