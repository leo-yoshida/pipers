library(ggplot2)
library('ggthemes')

level_order_parb <- c( '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '>14' )

ptub_tmhmm <- read.csv('ptub_r_tmhmm.txt', sep='\t')
parb_tmhmm <- read.csv('parb_r_tmhmm.txt', sep='\t')
pgau_tmhmm <- read.csv('pgau_r_tmhmm.txt', sep='\t')
pmal_tmhmm <- read.csv('pmal_r_tmhmm.txt', sep='\t')

roxo_claro <- "#BEAED4"
azulao <- "#084594"
teste <- "#CCECE6"
azul1 <- "#74A9CF"
azul2 <- "#3690C0"
azul3 <- "#0570B0"
azul4 <- "#034E7B"
ggplot(parb_tmhmm, aes(x=factor(tmhmm, level=level_order_parb), y=number)) +  
  geom_bar(stat="identity", width=0.7, fill="#3690C0", alpha=0.7) +
  geom_text(aes(label=number), vjust=-0.3, size=4.5,)+
  # coord_flip()+
  ggtitle("Piper arboreum") +
  ylab("Number of Transcripts\n") +
  xlab("\nNumber of Predicted Helices") +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000,6000), limits=c(0,6000), expand = expansion(mult = c(0, 0.05)))+
#  ylim(c(0, 6000))+
  #theme_minimal()+
  theme(plot.title = element_text(face="italic",vjust=-7,size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.key.size = unit(1.5,"line"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=15))

ggplot(pgau_tmhmm, aes(x=factor(tmhmm, level=level_order_parb), y=number)) +  
  geom_bar(stat="identity", width=0.7, fill="#3690C0", alpha=0.7) +
  geom_text(aes(label=number), vjust=-0.3, size=4.5)+
  # coord_flip()+
  ggtitle("Piper gaudichaudianum") +
  ylab("Number of Transcripts\n") +
  xlab("\nNumber of Predicted Helices") +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000,6000), limits=c(0,6000), expand = expansion(mult = c(0, 0.05)))+
#  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  #theme_minimal()+
  theme(plot.title = element_text(face="italic",vjust=-7,size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.key.size = unit(1.5,"line"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=15))

ggplot(pmal_tmhmm, aes(x=factor(tmhmm, level=level_order_parb), y=number)) +  
  geom_bar(stat="identity", width=0.7, fill="#3690C0", alpha=0.7) +
  geom_text(aes(label=number), vjust=-0.3, size=4.5)+
  # coord_flip()+
  ggtitle("Piper malacophyllum") +
  ylab("Number of Transcripts\n") +
  xlab("\nNumber of Predicted Helices") +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000,6000), limits=c(0,6000), expand = expansion(mult = c(0, 0.05)))+
#  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  #theme_minimal()+
  theme(plot.title = element_text(face="italic",vjust=-7,size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.key.size = unit(1.5,"line"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=15))

ggplot(ptub_tmhmm, aes(x=factor(tmhmm, level=level_order_parb), y=number)) +  
  geom_bar(stat="identity", width=0.7, fill="#3690C0", alpha=0.7) +
  geom_text(aes(label=number), vjust=-0.3, size=4.5)+
  # coord_flip()+
  ggtitle("Piper tuberculatum") +
  ylab("Number of Transcripts\n") +
  xlab("\nNumber of Predicted Helices") +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000), limits=c(0,4000), expand = expansion(mult = c(0, 0.05)))+
#  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  #theme_minimal()+
  theme(plot.title = element_text(face="italic",vjust=-7,size=22, hjust=0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.key.size = unit(1.5,"line"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=15))

signal_pipers <- read.csv('pipers_signal_table.txt', sep='\t')
signal_pipers

ggplot(signal_pipers, aes(x=species, y=number, fill=signalpeptide)) + 
  ggtitle("Predicted N-terminal peptide sequences") +
  ylab("Number of Transcripts\n") +
  xlab("\nSpecies\n\n") +
  labs(fill='Target') +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000), limits=c(0,5200), expand = expansion(mult = c(0, 0.05)))+
  geom_bar(stat="identity", width=0.7,alpha=0.7, position=position_dodge(width=0.8)) +
  geom_text(aes(label=number), position=position_dodge(width=0.8), vjust=-0.5, size=4)+
  theme(plot.title = element_text(vjust=1,size=20, hjust=0.5),
        axis.text.x = element_text(face='italic',size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.key.size = unit(1.5,"line"),
        legend.position="bottom",
        legend.title=element_text(size=18), 
        legend.text=element_text(size=15))

  
#geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+ #dodge deve ser maior que o width escolhido
 # geom_text(aes(label=number), position=position_dodge(width=0.8), hjust=1.5, size=4.5)+

