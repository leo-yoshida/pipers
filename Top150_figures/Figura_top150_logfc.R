library(ggplot2)
library(ggsci) 
library("scales")
library(RColorBrewer)

parb_barplot <- read.csv('arquivo_R_parb.txt', sep='\t')
pgau_barplot <- read.csv('arquivo_R_pgau.txt', sep='\t')
pmal_barplot <- read.csv('arquivo_R_pmal.txt', sep='\t')
ptub_barplot <- read.csv('arquivo_R_ptub.txt', sep='\t')

parb_barplot$category2 <- gsub("_", " ", parb_barplot$category)
pgau_barplot$category2 <- gsub("_", " ", pgau_barplot$category)
pmal_barplot$category2 <- gsub("_", " ", pmal_barplot$category)
ptub_barplot$category2 <- gsub("_", " ", ptub_barplot$category)

cinza_claro <- "#D9D9D9"
cinza_escuro <- "#666666"
teste <- "#CCECE6"

simp <- pal_simpsons("springfield")(16) #as cores em codigo
simp[[5]] <- "#F39B7FFF" #substituindo laranja claro por escuro
simp[[2]] <- "#9970AB"
simp[[3]] <- "#BDBDBD"
cinza_escuro <- "#666666"
verde_normal <- "#7FC97F"
rosa_estranho <- "#F0027F"
roxo_claro <- "#BEAED4"
amarelo <- "#FFFFCC"
azulao <- "#084594"

palete <- c(simp, cinza_escuro, verde_normal, rosa_estranho, roxo_claro,  amarelo, azulao)
nature <- pal_npg("nrc")(5)
simp_nature <- c(simp,nature)

set3 <- brewer.pal(n = 12, name = "Set3")
sum_set <- c(set2,set3)
sum_set

#nb.cols <- 21
#mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

level_order_full <- c( '(-16,-14)', '(-14,-12)', '(-12,-10)', '(-10,-8)', '(-8,-6)', '(-6,-4)', '(-4,-2)', '(-2,0)', '(0,2)', '(2,4)', '(4,6)', '(6,8)', '(8,10)', '(10,12)', '(12,14)', '(14,16)', '(16,18)' )

level_order_parb <- c( '(-14,-12)', '(-12,-10)', '(-10,-8)', '(-8,-6)', '(-6,-4)', '(-4,-2)', '(-2,0)', '(0,2)', '(2,4)', '(4,6)', '(6,8)', '(8,10)', '(10,12)', '(12,14)', '(14,16)', '(16,18)' )

ggplot(parb_barplot, aes(fill=category2, x = factor(range, level = level_order_parb), y=value)) + 
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = palete) +
  labs(title="Piper arboreum" , fill = "Category") +
  ylab("Number of Genes\n") +
  xlab("\nlogFC bin\n\n")+
  annotate("rect", xmin=0, xmax=7.5, ymin=0, ymax=Inf, alpha=0.1, fill="#F8766D")+
 # annotate("rect", xmin=1.5, xmax=2.5, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  annotate("rect", xmin=7.5, xmax=Inf, ymin=0, ymax=Inf, alpha=0.1, fill="#00BFC4")+
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = palete) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)) , breaks = seq(0, 50, by = 5))+
  scale_x_discrete(drop=F) +
  #coord_cartesian(xlim = c(1, 14)) +
  theme(plot.title = element_text(face="italic",size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom",
        legend.key.size = unit(1.8,"line"),
        #legend.title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=18))

library('Cairo')
ggsave("barplot_parb2.png", width = 6, height = 3, type = 'cairo', dpi= 300)

#PGAU
level_order_pgau <- c( '(-14,-12)', '(-12,-10)', '(-10,-8)', '(-8,-6)', '(-6,-4)', '(-4,-2)', '(-2,0)', '(0,2)', '(2,4)', '(4,6)', '(6,8)', '(8,10)', '(10,12)', '(12,14)' )

mycolors <- c("#F8766D","#00BFC4","#69b3a2")
ggplot(pgau_barplot, aes(fill=category2, x = factor(range, level = level_order_pgau), y=value)) + 
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = palete) +
  #scale_fill_manual("Legend", values = c("normal" = "red", "stress" = "green", "Nitrogen" = "blue"))+
  labs(title="Piper gaudichaudianum" , fill = "Category") +
  ylab("Number of Genes\n") +
  xlab("\nlogFC bin\n\n")+
  annotate("rect", xmin=0, xmax=7.5, ymin=0, ymax=Inf, alpha=0.1, fill="#F8766D")+
  # annotate("rect", xmin=1.5, xmax=2.5, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  annotate("rect", xmin=7.5, xmax=Inf, ymin=0, ymax=Inf, alpha=0.1, fill="#00BFC4")+
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = palete) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)) , breaks = seq(0, 50, by = 5))+
  scale_x_discrete(drop=F) +
  #coord_cartesian(xlim = c(1, 14)) +
  theme(plot.title = element_text(face="italic",size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom",
        legend.key.size = unit(1.8,"line"),
        #legend.title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=18))


#PMAL
level_order_pmal <- c( '(-14,-12)', '(-12,-10)', '(-10,-8)', '(-8,-6)', '(-6,-4)', '(-4,-2)', '(-2,0)', '(0,2)', '(2,4)', '(4,6)', '(6,8)', '(8,10)', '(10,12)', '(12,14)' )

ggplot(pmal_barplot, aes(fill=category2, x = factor(range, level = level_order_pmal), y=value)) + 
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = simp_nature) +
  labs(title="Piper malacophyllum", fill = "Category") +
  ylab("Number of Genes\n") +
  xlab("\nlogFC bin\n\n")+
  annotate("rect", xmin=0, xmax=7.5, ymin=0, ymax=Inf, alpha=0.1, fill="#F8766D")+
  annotate("rect", xmin=7.5, xmax=Inf, ymin=0, ymax=Inf, alpha=0.1, fill="#00BFC4")+
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = palete) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), , breaks = seq(0, 50, by = 5))+
  scale_x_discrete(drop=F) +
  coord_cartesian(xlim = c(1, 14)) +
  theme(plot.title = element_text(face="italic",size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom",
        legend.key.size = unit(1.8,"line"),
        #legend.title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=18))


#PTUB
level_order_ptub <- c( '(-14,-12)', '(-12,-10)', '(-10,-8)', '(-8,-6)', '(-6,-4)', '(-4,-2)', '(-2,0)', '(0,2)', '(2,4)', '(4,6)', '(6,8)', '(8,10)', '(10,12)', '(12,14)' )

ggplot(ptub_barplot, aes(fill=category2, x = factor(range, level = level_order_ptub), y=value)) + 
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = simp_nature) +
  #scale_fill_manual("Legend", values = c("normal" = "red", "stress" = "green", "Nitrogen" = "blue"))+
  labs(title="Piper tuberculatum", fill = "Category") +
  ylab("Number of Genes\n") +
  xlab("\nlogFC bin\n\n")+
  annotate("rect", xmin=0, xmax=7.5, ymin=0, ymax=Inf, alpha=0.1, fill="#F8766D")+
  # annotate("rect", xmin=1.5, xmax=2.5, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  annotate("rect", xmin=7.5, xmax=Inf, ymin=0, ymax=Inf, alpha=0.1, fill="#00BFC4")+
  geom_bar(position="stack", stat="identity", width=0.6, alpha=0.6) +
  scale_fill_manual(values = palete) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 50, by = 5))+
  scale_x_discrete(drop=F) +
  coord_cartesian(xlim = c(1, 14)) +
  theme(plot.title = element_text(face="italic",size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom",
        legend.key.size = unit(1.8,"line"),
        #legend.title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=18))
