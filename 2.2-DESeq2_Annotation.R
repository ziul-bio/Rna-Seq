######################## Análise de RNA-seq com DESeq2 #################################

# Autor: Luiz Carlos Vieira
# 25/09/2021

#------------------------- Import & pre-processamento ----------------------------------

library(DESeq2)

# Anotação gênica
library(AnnotationDbi)
library(org.Mm.eg.db)  #db de camundongo

# Manipulação e visualizaçãoo de dados
library(dplyr)
library(writexl)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Padronizando tamanho de imagens
options(repr.plot.width = 16, repr.plot.height = 12)

## Carregando a tabela com as counts por gene
countData <- read.table("featureCounts.txt", header=TRUE, row.names=1)

## Excluindo as colunas (chr, start, end, strand, length)
countData <- countData[,6:ncol(countData)]

## Renomeando as colunas:
colnames(countData) <- c('SRR10493810', 'SRR10493811', 'SRR10493818', 'SRR10493819')

## Convertendo o dataframe counsts em matriz
countData <- as.matrix(countData)

## Carregando a tabela coldata
coldata <- read.table("coldata.txt", header=TRUE, row.names=1, sep ='\t')

# Transformando a coluna "condition" em "factor" com dois levels (day4 e day7)
coldata$condition <- as.factor(coldata$condition)


#----------------------------- Análise com DESeq2 --------------------------------------


## Passando os dados de countData, coldata e design para a funÃ§Ã£o DESeqDataSetFromMatrix()
dds <- DESeqDataSetFromMatrix(countData=countData, colData=coldata, design= ~condition)

## Realizando um filtro para os genes que tenham a soma das readcounts &gt;= 10.
filtro <- rowSums(counts(dds)) >= 10
dds <- dds[filtro,]


### Definindo o nível de referência. 
dds$condition <- relevel(dds$condition, ref = "day4")

# Executando o deseq2
ddsDE <- DESeq(dds)


### Conferindo se a análise do DESeq foi realizada sobre as condições desejadas:
resultsNames(ddsDE)


# Gerando um objeto da função results() com um p-valeu ajustado < 0.05.
res05 <- results(ddsDE, alpha = 0.05, contrast=c("condition","day7","day4"), pAdjustMethod = "BH")


# -------------- Resumo dos resultados -------------------

# Número total de genes
sum(!is.na(res05$padj))

# Número de genes com um 'p-value ajustado' nemor que 0,05
sum(res05$padj < 0.05, na.rm=TRUE)

# Resumo dos resultados
summary(res05)


# Informaçõeses sobre quais variáveis e testes foram usados
mcols(res05)$description



#----------------------------- Visualização dos dados --------------------------------


# O Histograma abaixo mostra a distribuição das counts não normalizados.
png("No_norm_reads.png", width=800, height = 600, pointsize=20)
hist(countData, xlab = 'ReadsCounts não Normalizada', main= 'Dados não transformados')
dev.off()


# Regularized log transformation para os plots PCAs, clustering e heatmaps, 
rld <- rlog(ddsDE, blind = FALSE)
assay_rld <-assay(rld)

# histplot dos dados transformados
png("norm_reads.png", width=800, height = 600, pointsize=20)
hist(assay_rld, xlab = 'ReadsCounts Normalizada', main= 'Dados transformados' )
dev.off()


## Distribuição dos dados transformados
x <- assay_rld
corRainbow = rainbow(dim(x)[2])

png("density_counts.png", width=1000, height = 800, pointsize=20)
plot(density(x[,1]), col = corRainbow[1], lwd=2,
     xlab="valores de ExpressÃ£o", ylab="Densidade", main= "DistribuiÃ§Ã£o dos dados transformados",
     ylim=c(0, max(density(x[,1])$y)+.02 ) )
for( i in 2:dim(x)[2] )
lines(density(x[,i]), col=corRainbow[i], lwd=2)
legend("topright", cex=1.1, colnames(x), lty=rep(1,dim(x)[2]), col=corRainbow)
dev.off()


## Distribuição das read counts entre as replicatas.
png("dist_counts.png", width=1000, height = 800, pointsize=20)
par(mfrow=c(1,2))
plot(countData[,1:2], xlab='day4_rep1', ylab='day4_rep2', pch=16, cex=0.5)
plot(countData[,3:4], xlab='day7_rep1', ylab='day7_rep2', pch=16, cex=0.5)
mtext("Comparando as distribuição das replicatas", line = -1.5, cex = 1.5, outer = TRUE)
dev.off()


## Estimativa da DispersÃ£o
png("Plot_dispersions.png", width=1000, height = 800, pointsize=20)
plotDispEsts(ddsDE, main="Estimativa de DispersÃ£o")
dev.off()


## PCA
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("PCA.png", width=800, height = 400, pointsize=20)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=2) +
  ggtitle('Observação da variação entre as amostras do dia 4 e dia 7') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste0("PC1: ",percentVar[1],"% variação")) +
  ylab(paste0("PC2: ",percentVar[2],"% variação")) + 
  coord_fixed()
dev.off()


## PCA 2
png("PCA-2.png", width=800, height = 400, pointsize=20)
plotPCA(rld, intgroup = "condition", ntop = 500) +
  theme_bw() + # remove o tema padrÃ£o do ggplot2
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-5, 5)) + # Modifica os limites para corrigir as dimensÃµes da figura
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Variação dos genes em relação das amostras")
dev.off()


## Distância entre as amostras
sampleDists <- dist(t(assay_rld))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(assay_rld), " - ", rld$condition)
colnames(sampleDistMatrix) <- colnames(assay_rld)
corBlues <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("sample_Dist.png", width=900, height=600, pointsize=20)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=corBlues)
dev.off()


## Analisando os valores de p-value ajustado
png("p-values_plot.png", width=900, height=600, pointsize=20)
hist(res05$padj, breaks=50, col="gray", main =('Histrograma'),
     xlab = 'p-value ajustado', 
     ylab = 'Frequência')
dev.off()


## MA-Plot de LFCs vs contagens normalizadas
png("MA-plot.png", width=900, height=600, pointsize=20)
plotMA(res05, ylim=c(-5,5), main='Expressão diferenciada (day4 vs day7)', 
       xlab='Média das counts Normalizada')
dev.off()



#---------------------- shrink dos valores de log2 fold changes --------------------

## Consultando o resultsNames
resultsNames(ddsDE)
resultsNames(ddsDE)[2]

## Executando o lfcShrink()
res05_shrink <- lfcShrink(ddsDE, coef=2, res=res05, type = "apeglm")

## MA-Plot de LFCs dos lfc reduzidos
png("shrink.png", width=700, height=400, pointsize=20)
plotMA(res05_shrink,  ylim=c(-5,5))
dev.off()


## Comparando os plotMA em relação a os valores de LFC-shrink com os valores de LFC normais.
png("shrink_plots.png", width=900, height=600, pointsize=20)
par(mfrow=c(1,2), mar=c(4,4,4,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res05, xlim=xlim, ylim=ylim, main="normal")
plotMA(res05_shrink, xlim=xlim, ylim=ylim, main="apeglm")
dev.off()



## Plot counts
png("count_plot.png", width=900, height=600, pointsize=20)
par(mfrow=c(1,3), mar=c(6,4,4,1))
gene_idx1 = which.min(res05$padj)
plotCounts(ddsDE, gene=gene_idx1, intgroup="condition")
gene_idx2 = which.min(res05$log2FoldChange)
plotCounts(ddsDE, gene=gene_idx2, intgroup="condition")
gene_idx3 = which.max(res05$padj)
plotCounts(ddsDE, gene=gene_idx3, intgroup="condition")
dev.off()


### Agrupamento de genes pelo heatmap
topVarGenes <- head(order(-rowVars(assay_rld)),50) # - significa decresing
cores_heat<- colorRampPalette(brewer.pal(9, "YlOrRd"))(250)
mat <- assay_rld[ topVarGenes, ]
mat <- mat - rowMeans(mat)

png("heat_map.png", width=600, height=1000, pointsize=20)
pheatmap(mat, scale="row", cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, annotation_col=coldata, col=cores_heat)
dev.off()


## Volcano plot
png("volcano-DE.png", height=500, width=700, pointsize=20)
mylabel <- c('Significante', "Nãoo-Significante")
plot(res05_shrink$log2FoldChange, -1*log10(res05_shrink$padj), col=ifelse(res05_shrink$padj<0.05, "green", "black"),
     xlab="log Counts", ylab="log Fold Change", pch=20, main='Volcano Plot de Expressão Diferencial')
legend('topright', mylabel, fill=c('green', 'black'))
dev.off()



#-------------------------------- Anotação Gênica ----------------------------------

# Adicionando uma coluna de significância em res05_shrink
res05_shrink$significant <- ifelse(res05_shrink$padj< 0.05, "True", "False")
res05_shrink[which(abs(res05_shrink$log2FoldChange)< 0.5),'significant'] <- "False"


# Transformando a tabela res05 em df, para o merge com a tabela de dados homÃ³logos de camundongos
df_res05_shrink <- as.data.frame(res05_shrink)

# Ordenando a tabela res05 pelo valor de p-value ajustado
df_res05_shrink <- df_res05_shrink[order(df_res05_shrink$padj), ]
df_res05_shrink$Gene.stable.ID <- row.names(df_res05_shrink)


# Dataset com informaÃ§Ãµes de genes de camundongos homologos Ã¡s cÃ©lulas CHO (chinnese hamister ovary cells)
homoBiomart <- as.data.frame(read.delim('mart_export.txt', header = T, sep = '\t'))

# Unindo os datasets pela coluna em comun gene.stable.ID
results <- as.data.frame(merge(df_res05_shrink, homoBiomart, by='Gene.stable.ID'))

# Removendo dados duplicados da Mouse.gene.stable.ID
results <- distinct(results, Mouse.gene.stable.ID, .keep_all= TRUE)


# Informaçõess sobre os dados de entrada para a query.

### Consultar as colunas do db
columns(org.Mm.eg.db)

### Consultar as keys do db
keytypes(org.Mm.eg.db)


#--------------- Adicionando as anotaçõees na tabela de resultados -------------------

# Adicionando o símbolo do gene
results$GeneSymbol <- mapIds(x = org.Mm.eg.db,
                              keys = results$Mouse.gene.stable.ID,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")


# Adicionando a descrição do gene
results$GeneDescription <- mapIds(org.Mm.eg.db,
                              keys = results$Mouse.gene.stable.ID,
                              column = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Adicionando o ENTREZID
results$ENTREZID <- mapIds(org.Mm.eg.db,
                              keys = results$Mouse.gene.stable.ID,
                              column = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# --------------------------------------------------------------------------------------

# Criando um subset dos genes com expressÃ£o significativa
res_sig = as.data.frame(subset(results, padj< 0.05))
res_sig <- na.omit(res_sig)
res_sig = res_sig[order(res_sig$log2FoldChange, decreasing=TRUE),]


# Volcano Plot
volcano = ggplot(res_sig, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=significant)) +
  scale_color_manual(values=c("red", "green"))

png("volcano-DE2.png", height=500, width=700, pointsize=20)
volcano + geom_text_repel(data=filter(res_sig, abs(log2FoldChange) > 2 & padj <1e-50), aes(label=GeneSymbol)) +
  coord_cartesian(clip = "off")
dev.off()



#--------------------------- Salvando o arquivo resultados ----------------------------

# Criando uma cópia do df
df <- res_sig

# Aplicando a funçãoo as.character as colunas GeneSymbol, GeneDescription e ENTREZID c(10:12)
alter_cols <- apply(df[ , c(9:11)], 2, as.character)  # retorna um df com as colunas modificadas

# Substituindo as colunas
df[ , colnames(df) %in% colnames(alter_cols)] <- alter_cols  

# Salvando os resultados do DESeq2
# write.csv(as.data.frame(df), 'results_DESeq2_shrink.csv')
# writexl::write_xlsx(as.data.frame(df), 'results_DESeq2_shrink.xlsx')
# 
# # Sanvando os resultados DE para e up regulados down regulados 
# res05_shrink_sig_up <- filter(df, significant == "True" & log2FoldChange > 0)
# writexl::write_xlsx(as.data.frame(res05_shrink_sig_up), 'results_DESeq2_shrink_genes_up.xlsx')
# 
# res05_shrink_sig_down <- filter(df, significant == "True" & log2FoldChange < 0)
# writexl::write_xlsx(as.data.frame(res05_shrink_sig_down), 'results_DESeq2_shrink_genes_down.xlsx')



## Informações da sessão
sessionInfo()

# FIM