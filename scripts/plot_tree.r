#tree
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("YuLab-SMU/treedataverse")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")

library(treedataverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)

#plot
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(reshape2)
library(ggrepel)

tree <- read.tree("./data/all62.nwk")
anno <- read.table("./data/tree_ano",header=T,sep ="\t")
p=ggtree(tree)+geom_tiplab(align = T)
#gheatmap(p, anno["typeABS"], offset=0, width=0.2)
library(aplot)
p=ggtree(tree,layout = "fan")
library(ggstar)
p1=p + geom_fruit(
          data=anno,
          geom=geom_star,
          mapping=aes(y=`node`,fill=typeABS, size=year, starshape=typeABS),
          position="identity",
          starstroke=0.1
      )+ 
    # Adjust star size and its legend
      scale_size_continuous(
          range=c(1, 3), # the range of size.
          guide=guide_legend(
                    keywidth=0.5, 
                    keyheight=0.5,
                    override.aes=list(starshape=15),
                    order=2 # order of legend
                )
      ) +
    # Adjust star color and its legend
      scale_fill_manual(
          values=c("#8dd3c7", "#ffed6f"),
          guide="none" 
      ) 
p1 "#8dd3c7", "#ffed6f", "#bebada"
library(ggnewscale)
p2=p1 + 
      ggnewscale::new_scale_fill() + 
      geom_fruit(data=anno,geom=geom_tile,mapping=aes(y=node, fill=typeRELA),
          offset=0.08,   # Adjust
          pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
      ) + 
      scale_fill_manual(
          values=c("#fdbf6f", "#fb9a99", "#a6bce3"))

p2 
p3=p2+ ggnewscale::new_scale_fill() +
      geom_fruit(
          data=anno,
          geom=geom_col,
          mapping=aes(y=node, x=year, fill=typeABS), 
          pwidth=0.4,offset = 0.1,
          axis.params=list(
                          axis="x", # Add x-axis text
                          text.size=2, # text size
                          text.angle=-45, # angle
                          hjust=0  # Adjust
                      ),
          grid.params=list() # Add grid lines
      ) + 
      scale_fill_manual(
          values=c("#bebada", "#fdbf6f", "#fb9a99"),
          guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)
      )
#bebada
ggsave(p3, file='tree_anno.pdf', width=10, height=10)


p4 = p3+ggnewscale::new_scale_fill() + 
      geom_fruit(data=anno,geom=geom_tile,mapping=aes(y=node, fill=order),
          offset=0.08,   # Adjust
          pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
      ) + 
      scale_fill_manual(
          values=c("#947A6D", "#D7B98E", "#a6bce3","#808F7C","#019a99","#0077b0","#ffba4d","#282152","#caa59a","#fdbf6f","#69B0AC","#EEA9A9","#854836","#E1A679","#2B5F75","#D19826","#aac8eb","#caa59a","#DC9FB4"))

ggsave(p4, file='tree_anno_order.pdf', width=10, height=16)

