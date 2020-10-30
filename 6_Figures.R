source("X_Functions.R")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("lemon")


# Read Data In ------------------------------------------------------------
Cluster_Nobt <- readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
Nobt_Labs <- ClusterLabs(Cluster_Nobt)

# Tobacco Clusters -----------------------------------------------------------
C1 <- ggplot(Cluster_Nobt$normalized, 
                      aes(x=DAP, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression",
       x="Stage") +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=Nobt_Labs),
                 nrow = 2) +
  scale_fill_manual(values=palfill[2]) +
  scale_color_manual(values=palline[2]) +
  scale_x_discrete(labels=c("0"="1", "3"="2", "6"="3", "11"="Br")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none",
        legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))


# Tobacco GO Plots --------------------------------------------------------
# Do the GO Enrichment
Tables_Nobt <- list()
for (i in levels(as.factor(Cluster_Nobt$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% Cluster_Nobt$normalized$genes[Cluster_Nobt$normalized$cluster==i]))
  names(tmpList) <- Cluster_Nobt$normalized$genes
  Tables_Nobt[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Nobt.gene2go.tsv",
                                GOIs=tmpList)
}
names(Tables_Nobt) <- Nobt_Labs[names(Tables_Nobt)]
# Create a list of the GO Plots
GO_Nobt <- list()
GO_Nobt <- lapply(seq_along(Tables_Nobt), 
                  function(i) GOPlot(Tables_Nobt[[i]], Title = paste("Cluster ",i, " GO Enrichment"))+
                    theme(legend.position = "none"))

# Solanum Comparison ------------------------------------------------------


# Five-Species Venn Diagram -----------------------------------------------


# Five-Species GO Plots ---------------------------------------------------


# Five-Species PCA --------------------------------------------------------


# Five-Species Clusters ---------------------------------------------------


