library(ggplot2)
x <- read.table("Lcancer_Scancer",col.names = c("class","rate"))
p <- ggplot(x, aes(x =rate,fill = class)) +
     geom_density(alpha = 0.8)+
     theme(panel.background = element_blank(), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line =element_line(colour = "black"))+
     scale_x_continuous(limits = c(0,1))+
     xlab("Root-to-Tip distance of cancer")+ 
     ylab("Density of evolutionary distance")+
     #scale_fill_wsj()+
     #scale_fill_d3()+
     theme(legend.position = c(0.8,0.8),legend.text=element_text(size=16),
           axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
     scale_fill_manual("Label title",breaks = c("Lcancer", "Scancer"),
                       values = c("Lcancer"="#C2D9EB", "Scancer"="#FED498"))
ggsave(p, file='Root-to-Tip.cancer.long.short.Density.pdf', width=5, height=4.5)

library(ggplot2)
library(ggprism)
x <- read.table("raxml.cancerlong.cancershort",header=F,sep="\t",col.names = c("class","rate"))
p <- ggplot(x, aes(x = class, y = rate, fill = class)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c("#C2D9EB", "#FED498")) +
  labs(x = " ", y = "root to tip distance", cex.lab = 20, color = "Order") +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 1, linetype = 1),
    axis.ticks = element_line(colour = "black", size = 1, linetype = 1),
    axis.ticks.length = unit(3, "mm"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) 
ggsave(p, file='Root-to-Tip.cancer.long.short.box.pdf', width=5, height=4.5)
