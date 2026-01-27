#install.packages("gcookbook")
library(gcookbook)
library(ggplot2)
library(dplyr)
a <- read.table("./data/polyseldraw",header = T,stringsAsFactors = F,sep ="\t")
head(a)
p <- ggplot(a,mapping = aes(x=log10p,y=Pathway, colour = setScore))+
  geom_point(size = a$setSize*0.2) +
  scale_colour_continuous(low = 'lightblue', high = 'darkblue', breaks = c(10,20,30,40,50,60,70))+
  theme_gray()+
  theme(axis.text.x = element_text(size=12,angle=0,face="plain"),
        axis.text.y = element_text(size=12,angle=0,face="plain"),  
        axis.title.x = element_text(size=12,face="plain"),
        axis.title.y = element_text(size=12,face="plain"))+
        xlab("log(p)")+
        ylab("Pathway of the cancer mutations are more conserved genes in long-lived species")
  
  #+
 # theme(panel.background = element_blank(), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line =element_line(colour = "black"))
 
 ggsave(p, file='Pathway of the cancer mutations are more conserved genes in long-lived species.pdf', width=10, height=6)
