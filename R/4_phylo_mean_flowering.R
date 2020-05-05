library(phytools)
library(geiger)

phy <- read.tree("Poaceae_hyp1_allC4_maxcred.phy")
mean <- read.csv("data/slopes/mean_flowering.csv")

rownames(mean) <- mean$binomial
td <- treedata(phy, mean)
x <- mean$sp_mean
names(x) <- mean$binomial
se <- mean$sp_sterr
names(se) <- mean$binomial

# estimate phylogenetic signal
mean_psigk <- phylosig(td$phy, x = x, method="K", test=TRUE, nsim=1000, se = se, start=NULL,control=list())
mean_psigl <- phylosig(td$phy, x = x, method="lambda", test=TRUE, nsim=1000, se = se, start=NULL,control=list())

## plot phylogeny depicting mean flowering date ####
  # colours genereated from scales::show_col(viridis(11, begin = 0.3, end = 0.9))
ygre = c("#bbdf27", "#93d741","#6ece58", "#4dc36b", "#33b67a", "#22a884", "#1f9a8a", "#228b8d", "#287d8e", "#2e6f8e", "#35608d")
trait.color = c()
ind = c()
  # pick colour based on mean flowering date of species
for(i in 1:length(td$phy$tip.label)){
  ind <- td$data[,2][i] < c(160,170,180,190,200,210,220,230,240,250,260)
  trait.color[i] = ygre[which(ind == TRUE)[1]]
  }
names(trait.color) <- td$data[,1]
trait.color <- trait.color[order(factor(names(trait.color), levels = td$phy$tip.label))]
  # save phy plot with species names (coloured) to depict trait value
png("figures/fig_1.png", units="in", width=8, height=10, res=300)
plot.phylo(td$phy,tip.color = trait.color, align.tip.label = TRUE,  cex =0.5, label.offset = 1, font = 1)
  # manually add plot color legend
text(x =1.1,y =55,"260", cex = 0.7)
text(x =1.1,y =71,"160", cex = 0.7)
text(x =6.7,y =71,"June", cex = 0.7)
text(x =6.5,y =65.3,"July", cex = 0.7)
text(x =7.3,y =60.3,"August", cex = 0.7)
text(x =8.4,y =55,"September", cex = 0.7)
par(fig=c(0.1415,0.2115, 0.68, 0.908), mar=c(1,1,1,1), new = T)
image(y = 160:260,x = 1, z = as.matrix(t.default(-160:-260)),col = viridis(100, begin = 0.9, end = 0.3), xlab = "", ylab = "", xaxt = "n" ,las = 1, axes = F)
dev.off()
