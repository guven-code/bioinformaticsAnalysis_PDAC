library(ggplot2)



terms<-top10$Term
term2=gsub("GO:[0-9]{7}~","",terms)
top10$Term<-term2






ggsave("Goplots.pdf")#width=15,height = 10,units="cm")

options(repr.plot.width=5, repr.plot.height=3)
ggplot(top10[1:40,], aes(x = Term, y = Count, main="GO and KEGG terms")) +
  geom_bar(stat = "identity") +
  coord_flip() + scale_y_continuous(name="Gene Count") +
  scale_x_discrete(name="") +
  theme(axis.text.x = element_text(face="bold", color="black",
                                   size=8, angle=0),
        axis.text.y = element_text(face="bold", color="black",
                                   size=8, angle=0))
dev.off()
#ggplot(bp_fat15, aes(x = Term,
#                     y = -log10(PValue),
 #                    fill = gene_Count)) + geom_bar(stat = "identity",position = "dodge")+ 
 # scale_color_gradient(low="blue",
      #                  high="red")+coord_flip()