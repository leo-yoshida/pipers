#Nome das amostras e condição experimental
samples_parb <- read.table("samples_parb.txt", header=TRUE)
rownames(samples_parb) <- samples_parb$sample

samples_pgau <- read.table("samples_pgau.txt", header=TRUE)
rownames(samples_pgau) <- samples_pgau$sample

samples_pmal <- read.table("samples_pmal.txt", header=TRUE)
rownames(samples_pmal) <- samples_pmal$sample

samples_ptub <- read.table("samples_ptub.txt", header=TRUE)
rownames(samples_ptub) <- samples_ptub$sample

#Diretórios onde a quantificação do salmon se encontra
dir_parb <- "salmon_results_parb"
dir_pgau <- "salmon_results_pgau"
dir_pmal <- "salmon_results_pmal"
dir_ptub <- "salmon_results_ptub"

files_parb <- file.path(dir_parb, samples_parb$sample, "quant.sf")
files_pgau <- file.path(dir_pgau, samples_pgau$sample, "quant.sf")
files_pmal <- file.path(dir_pmal, samples_pmal$sample, "quant.sf")
files_ptub <- file.path(dir_ptub, samples_ptub$sample, "quant.sf")

names(files_parb) <- samples_parb$sample
names(files_pgau) <- samples_pgau$sample
names(files_pmal) <- samples_pmal$sample
names(files_ptub) <- samples_ptub$sample

#Arquivos tx2gene, lido com 'readr'
library("readr")
tx2gene_parb <- read_csv("E:/Doutorado_Destilado_2024/Transcriptomas_Pipers_2025/tx2gene/cat_parb_tx2gene.txt")
tx2gene_pgau <- read_csv("E:/Doutorado_Destilado_2024/Transcriptomas_Pipers_2025/tx2gene/cat_pgau_tx2gene.txt")
tx2gene_pmal <- read_csv("E:/Doutorado_Destilado_2024/Transcriptomas_Pipers_2025/tx2gene/cat_pmal_tx2gene.txt")
tx2gene_ptub <- read_csv("E:/Doutorado_Destilado_2024/Transcriptomas_Pipers_2025/tx2gene/cat_ptub_tx2gene.txt")

#tximport, importa os counts (números quebrados do salmon) de transcritos do salmon, e sumariza em termos de genes
library("tximport")
txi_parb <- tximport(files_parb, type="salmon", tx2gene=tx2gene_parb)
txi_pgau <- tximport(files_pgau, type="salmon", tx2gene=tx2gene_pgau)
txi_pmal <- tximport(files_pmal, type="salmon", tx2gene=tx2gene_pmal)
txi_ptub <- tximport(files_ptub, type="salmon", tx2gene=tx2gene_ptub)

txi_parb_isoforms <- tximport(files_parb, type="salmon", txOut = TRUE)
txi_pgau_isoforms <- tximport(files_pgau, type="salmon", txOut = TRUE)
txi_pmal_isoforms <- tximport(files_pmal, type="salmon", txOut = TRUE)
txi_ptub_isoforms <- tximport(files_ptub, type="salmon", txOut = TRUE)

#Modelagem e quantificação dos counts com DESeq2
library("DESeq2")

#Função DESeqDataSetFromTximport() converte counts continuos em counts discretos, cria um objeto dds
ddsTxi_parb <- DESeqDataSetFromTximport(txi_parb,
                                        colData = samples_parb,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_parb
ddsTxi_pgau <- DESeqDataSetFromTximport(txi_pgau,
                                        colData = samples_pgau,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_pgau
ddsTxi_pmal <- DESeqDataSetFromTximport(txi_pmal,
                                        colData = samples_pmal,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_pmal
ddsTxi_ptub <- DESeqDataSetFromTximport(txi_ptub,
                                        colData = samples_ptub,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_ptub

#Isoforms
ddsTxi_parb <- DESeqDataSetFromTximport(txi_parb_isoforms,
                                        colData = samples_parb,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_parb
ddsTxi_pgau <- DESeqDataSetFromTximport(txi_pgau_isoforms,
                                        colData = samples_pgau,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_pgau
ddsTxi_pmal <- DESeqDataSetFromTximport(txi_pmal_isoforms,
                                        colData = samples_pmal,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_pmal
ddsTxi_ptub <- DESeqDataSetFromTximport(txi_ptub_isoforms,
                                        colData = samples_ptub,
                                        design = ~ condition) #design identifica a coluna "condition" no samples_ptub

#Com os counts discretos, primeiramente será realizado um PCA com os counts estabilizados em relação à variância
#Na escala log, genes com baixos níveis de expressão exibem maior variancia, por isso o rlog é aplicado para estabilizar essa variancia
#Torna os dados de counts homocedásticos

rlog_dds_parb <- rlog(ddsTxi_parb_isoform)
rlog_dds_pgau <- rlog(ddsTxi_pgau_isoform)
rlog_dds_pmal <- rlog(ddsTxi_pmal_isoform)
rlog_dds_ptub <- rlog(ddsTxi_ptub_isoform)

library(ggplot2)
#PCA para checar similaridade das amostras
plotPCA(rlog_dds_parb, intgroup="condition") + ggtitle('Piper arboreum') + theme(plot.title = element_text(face='italic', hjust=0.5, size=20)) + theme(aspect.ratio = 1) + theme(axis.text.y = element_text(colour = "black", size = 12)) +  theme(axis.text.x = element_text(colour = "black", size = 12)) + theme(axis.title.x = element_text(size=16)) + theme(axis.title.y = element_text(size=16)) + theme(legend.title = element_text(size=14), legend.text = element_text(size = 14))
plotPCA(rlog_dds_pgau, intgroup="condition") + ggtitle('Piper gaudichaudianum') + theme(plot.title = element_text(face='italic', hjust=0.5, size=20)) + theme(aspect.ratio = 1) + theme(axis.text.y = element_text(colour = "black", size = 12)) +  theme(axis.text.x = element_text(colour = "black", size = 12)) + theme(axis.title.x = element_text(size=16)) + theme(axis.title.y = element_text(size=16)) + theme(legend.title = element_text(size=14), legend.text = element_text(size = 14))
plotPCA(rlog_dds_pmal, intgroup="condition") + ggtitle('Piper malacophyllum') + theme(plot.title = element_text(face='italic', hjust=0.5, size=20)) + theme(aspect.ratio = 1) + theme(axis.text.y = element_text(colour = "black", size = 12)) +  theme(axis.text.x = element_text(colour = "black", size = 12)) + theme(axis.title.x = element_text(size=16)) + theme(axis.title.y = element_text(size=16)) + theme(legend.title = element_text(size=14), legend.text = element_text(size = 14))
plotPCA(rlog_dds_ptub, intgroup="condition") + ggtitle('Piper tuberculatum') + theme(plot.title = element_text(face='italic', hjust=0.5, size=20)) + theme(aspect.ratio = 1) + theme(axis.text.y = element_text(colour = "black", size = 12)) +  theme(axis.text.x = element_text(colour = "black", size = 12)) + theme(axis.title.x = element_text(size=16)) + theme(axis.title.y = element_text(size=16)) + theme(legend.title = element_text(size=14), legend.text = element_text(size = 14))

#Clusterização também para checar similaridade das amostras
#colData() em um dds mostra o sampledata embutido no DESeqDataSetFromTximport()
colData(rlog_dds_parb)
colData(rlog_dds_parb)$condition
plot(hclust(dist(t(assay(rlog_dds_parb)))), main = expression(italic("Piper arboreum:")~"Cluster dendogram") , labels=colData(ddsTxi_parb)$condition, cex=1.4, cex.lab=1.4, cex.axis=1.1, cex.main=1.7,sub=NA,xlab=NA)
plot(hclust(dist(t(assay(rlog_dds_pgau)))), main = expression(italic("Piper gaudichaudianum:")~"Cluster dendogram") , labels=colData(ddsTxi_pgau)$condition, cex=1.4, cex.lab=1.4, cex.axis=1.1, cex.main=1.7,sub=NA,xlab=NA)
plot(hclust(dist(t(assay(rlog_dds_pmal)))), main = expression(italic("Piper malacophyllum:")~"Cluster dendogram") , labels=colData(ddsTxi_pmal)$condition, cex=1.4, cex.lab=1.4, cex.axis=1.1, cex.main=1.7,sub=NA,xlab=NA)
plot(hclust(dist(t(assay(rlog_dds_ptub)))), main = expression(italic("Piper tuberculatum:")~"Cluster dendogram") , labels=colData(ddsTxi_ptub)$condition, cex=1.4, cex.lab=1.4, cex.axis=1.1, cex.main=1.7,sub=NA,xlab=NA)

#Expressão diferencial DEG

#Primeiro, checar o desenho experimental incluido na hora de criar objeto dds
design(ddsTxi_parb)
design(ddsTxi_pgau)
design(ddsTxi_pmal)
design(ddsTxi_ptub)

ddsTxi_parb$condition #Seedling está por último, é necessário trocar
ddsTxi_parb$condition <- relevel(ddsTxi_parb$condition, "Seedling")
ddsTxi_parb$condition #Agora "Adult" está por último

ddsTxi_pgau$condition #Seedling está por último, é necessário trocar
ddsTxi_pgau$condition <- relevel(ddsTxi_pgau$condition, "Seedling")
ddsTxi_pgau$condition #Agora "Adult" está por último

ddsTxi_pmal$condition #Seedling está por último, é necessário trocar
ddsTxi_pmal$condition <- relevel(ddsTxi_pmal$condition, "Seedling")
ddsTxi_pmal$condition #Agora "Adult" está por último

ddsTxi_ptub$condition #Seedling está por último, é necessário trocar
ddsTxi_ptub$condition <- relevel(ddsTxi_ptub$condition, "Seedling")
ddsTxi_ptub$condition #Agora "Adult" está por último

#Função DESeq() nos counts discretos para analise de expressão diferencial
ddsTxi_parb <- DESeq(ddsTxi_parb)
res_parb <- results(ddsTxi_parb)

ddsTxi_pgau <- DESeq(ddsTxi_pgau)
res_pgau <- results(ddsTxi_pgau)

ddsTxi_pmal <- DESeq(ddsTxi_pmal)
res_pmal <- results(ddsTxi_pmal)

ddsTxi_ptub <- DESeq(ddsTxi_ptub)
res_ptub <- results(ddsTxi_ptub)


#Resultados com logFC (lfcThreshold) diferente e valores de p-valor ajustados (alpha) diferentes 
res_parb <- results(ddsTxi_parb, lfcThreshold = 1, alpha=0.05)
res_pgau <- results(ddsTxi_pgau, lfcThreshold = 1, alpha=0.05)
res_pmal <- results(ddsTxi_pmal, lfcThreshold = 1, alpha=0.05)
res_ptub <- results(ddsTxi_ptub, lfcThreshold = 1, alpha=0.05)

#Salvando os resultados em um dataframe
resultados_parb <- data.frame(res_parb)
resultados_pgau <- data.frame(res_pgau)
resultados_pmal <- data.frame(res_pmal)
resultados_ptub <- data.frame(res_ptub)

#Salvando os resultados em .tsv
write.table(resultados_parb, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(resultados_pgau, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(resultados_pmal, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(resultados_ptub, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")

#Tabela counts normalizados
parb_normalized_counts <- counts(ddsTxi_parb, normalized=TRUE)
pgau_normalized_counts <- counts(ddsTxi_pgau, normalized=TRUE)
pmal_normalized_counts <- counts(ddsTxi_pmal, normalized=TRUE)
ptub_normalized_counts <- counts(ddsTxi_ptub, normalized=TRUE)

write.table(parb_normalized_counts, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(pgau_normalized_counts, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(pmal_normalized_counts, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(ptub_normalized_counts, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")

#Tabela counts normalizados rlog

rlog_parb <- assay(rlog_dds_parb)
rlog_pgau <- assay(rlog_dds_pgau)
rlog_pmal <- assay(rlog_dds_pmal)
rlog_ptub <- assay(rlog_dds_ptub)

write.table(rlog_parb, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(rlog_pgau, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(rlog_pmal, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
write.table(rlog_ptub, file = "Path_to_output", row.names=TRUE, col.names=NA, sep="\t")
