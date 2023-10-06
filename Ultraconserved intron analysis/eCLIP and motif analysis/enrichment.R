library(ggplot2)
df2 <- read.table("enriched_motifs_dif_1+.txt", sep="\t", header=TRUE)
input = c('plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'bold')
plot1 = ggplot(df2, aes(x=1, y=RBP, label=RBP, fill=Enrichment))+
  geom_tile()+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.text = element_text(size=5), plot.title = element_text(size=10, hjust = 0.5),
                     panel.grid.minor = element_blank(), axis.text = element_text(size=5, face = input), legend.title = element_text(size=10), legend.key.size = unit(0.75, 'cm'))+
  coord_equal()+
  theme( axis.title.x=element_blank())+
  ggtitle('Conserved NMD')+
  labs(fill = 'predicted motifs/region difference from constitutive exon controls')+
  theme( axis.text.x=element_blank())+
  theme( axis.ticks.x=element_blank())+
  theme( axis.ticks.y=element_blank())+
  theme( axis.title.y=element_blank())+
  scale_y_discrete(limit=rev(df2$RBP))+
  scale_fill_gradient2(midpoint=0, low="dodgerblue4", mid="white",
                        high="red", limits = c(-5.5,10))#min(df2$Enrichment),max(df2$Enrichment)))
plot1
ggsave(plot = plot1, width = 5, height = 15, dpi = 1700, filename = "top eclip cons 5%+.png")
