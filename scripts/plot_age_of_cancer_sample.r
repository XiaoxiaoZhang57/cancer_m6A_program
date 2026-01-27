x <- read.table("./data/alltissue.age1",col.names = c("tissue","age"),sep="\t")
library(ggplot2)
library(ggridges)
p <- ggplot(data = x, mapping = aes(
  x = age, 
  y = tissue,
  fill = tissue))
p1 <- p + geom_density_ridges(alpha = 0.4) +
guides(fill = FALSE) +
labs(x = "Age",
y = "Density")+
theme(panel.background = element_blank(), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line =element_line(colour = "black"))+
xlab("The age of the cancer sample")+
ylab("Cumulative distribution density for each age")+
theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
ggsave(p1, file='The age of the cancer sample.Density.pdf', width=5, height=6)
