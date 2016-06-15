setwd("/groups2/joshi_grp/guillaume/otherProject/earMicroarrays/")
library(limma)
library(gplots)
library(ggplot2)
library(dplyr)
library(magrittr)
library(gridExtra)
library(data.table)
library(parallel)
library(readr)
library(cowplot)
load("rma_data_branArrayENST.RData")

metaData <- read.table("metaData.txt", header = TRUE, row.names = 1)
metaData[23, 1] <- 2
metaData$Tissue <- factor(metaData$Tissue)
levels(metaData$Tissue) <- c("blood", "serous", "mucoid")
row.names(metaData) <- sub(".CEL", "", row.names(metaData))
lExprs$data <- lExprs$data[, row.names(metaData)]
rownames(lExprs$data) <- 1:nrow(lExprs$data)
rownames(lExprs$anno) <- 1:nrow(lExprs$anno)

ENSTdico <- read_tsv("ENST2symbol_bioMart_2015-10-20.txt")
ENSTdico <- as.data.frame(ENSTdico, stringAsFactors = FALSE)
colnames(ENSTdico) <- c("geneID", "transcriptID", "description", "geneName", "transcriptName", "transcriptType")
ENSTdico <- ENSTdico[,c(1,2,4,5,6,3)]
ENSTdico <- unique(ENSTdico)

lExprs$anno$transcriptID <- sub("_at", "", lExprs$anno$unitName, fixed = TRUE)

dim(lExprs$anno)
lExprs$anno <- merge(lExprs$anno, ENSTdico, by = "transcriptID", all.x = TRUE, all.y = FALSE)
dim(lExprs$anno )

doMeanOfTranscriptsForGene <- function(myGene) {
    lExprsSubset <- lapply(lExprs, function(x) subset(x, lExprs$anno$geneName == myGene & lExprs$anno$transcriptType == "protein_coding"))
    lExprsSubset$anno <- lExprsSubset$anno[1, c("geneID", "geneName")]
    lExprsSubset$data <- colMeans(lExprsSubset$data)
    finalLine <- c(lExprsSubset$anno, lExprsSubset$data) %>% as.data.frame
    return(finalLine)
}

a <- Sys.time()
mergedTranscripts <- do.call(rbind, mclapply(unique(lExprs$anno$geneName), doMeanOfTranscriptsForGene, mc.cores = 32))
Sys.time() - a # long : 5minutes
dim(mergedTranscripts)
mergedTranscripts <- subset(mergedTranscripts, !is.nan(mergedTranscripts[,3]))
row.names(mergedTranscripts) <- mergedTranscripts$geneName
dim(mergedTranscripts)
tail(mergedTranscripts)

lExprs$anno <- mergedTranscripts[,1:2]
lExprs$data <- mergedTranscripts[,3:ncol(mergedTranscripts)]

# log and logFC transformation -------------
any(rowMeans(lExprs$data) == 0) # FALSE
lExprs$logFC_all <- log2(lExprs$data / rowMeans(lExprs$data))
lExprs$log2 <- log2(lExprs$data)

# limma DGE analysis -----------------
design <- model.matrix(~0+metaData$Tissue)
colnames(design) <- c("blood", "serous", "mucoid")
fit <- lmFit(lExprs$log2, design = design)

contrast.matrix <- makeContrasts(serous - blood, mucoid - blood, mucoid - serous, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH", number = 10)

# venn diagrams ----------------
pdf(file = "Rplot/brainArray_tm/DGE_vennDiagrams.pdf", width = 9, height = 6)
layout(matrix(1:6, nrow = 2, byrow = TRUE))
results <- decideTests(fit2)
vennDiagram(results, main = "Up or Down, p-value < 0.05, |log2FC| > 0", cex = c(1,1,1), mar = c(0,0,1,0))
vennDiagram(results, main = "Up, p-value < 0.05, log2FC > 0", include = "up", cex = c(1,1,1), mar = c(0,0,1,0))
vennDiagram(results, main = "Down, p-value < 0.05, log2FC < 0", include = "down",  cex = c(1,1,1), mar = c(0,0,1,0))
results <- decideTests(fit2, p.value = 0.01, lfc = 1)
vennDiagram(results, main = "Up or Down, p-value < 0.01, |log2FC| > 1", cex = c(1,1,1), mar = c(0,0,1,0))
vennDiagram(results, main = "Up, p-value < 0.01, log2FC > 1", include = "up", cex = c(1,1,1), mar = c(0,0,1,0))
vennDiagram(results, main = "Down, p-value < 0.01, log2FC < -1", include = "down", cex = c(1,1,1), mar = c(0,0,1,0))
dev.off()

system("firefox Rplot/brainArray_tm/DGE_vennDiagrams.pdf &")

# heatmaps -------------------
results <- decideTests(fit2, p.value = 0.01, lfc = 1)
degl <- rowSums(abs(results))
degl <- names(degl)[which(degl > 0)]

colByTissue <- factor(metaData$Tissue)
levels(colByTissue) <- c("green", "red", "purple")
colByTissue <- as.character(colByTissue)

MvsS <- names(results[,3][results[,3] != 0])
colMvsS <- rep("white", length(degl))
names(colMvsS) <- degl
colMvsS[degl %in% MvsS] <- "red"

any(rowMeans(lExprs$log2[,metaData$Tissue == "blood"]) == 0)
lExprs$log2FC_blood <- lExprs$log2 - rowMeans(lExprs$log2[,metaData$Tissue == "blood"])

mybreaks <- seq(-3, 3, length.out = 128)
myGradient <- colorRampPalette(c("green", "black", "magenta"))(127)

pdf(file = "Rplot/brainArray_tm/heatmap_all_DGE_p0.01_fc1.pdf", width = 8.3, height = 11.7)
heatmap.2(as.matrix(lExprs$log2FC_blood[degl,]), Colv = TRUE, dendrogram = "both", breaks = mybreaks, col = myGradient,
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          trace = "none", ColSideColors = colByTissue, RowSideColors = colMvsS, 
          density.info = "none", key.xlab = "log2(FC)", denscol = "white", key.title="",
          labRow = "", labCol = sub("_.HuGene.2_0.st.", "", colnames(lExprs$log2FC_blood), fixed = TRUE),
          lhei = c(1, 5), main = "Supplementary figure 3:\n2682 differentially expressed genes", useRaster = TRUE)
dev.off()    

# serous vs mucoid
degl2 <- names(which(results[,3] !=0))

pdf(file = "Rplot/brainArray_tm/heatmap_sereous_vs_mucoid_DGE_p0.01_lfc1.pdf", width = 8.3, height = 11.7)
heatmap.2(as.matrix(lExprs$log2FC_blood[degl2,]), Colv = TRUE, dendrogram = "both", breaks = mybreaks, col = myGradient,
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          trace = "none", ColSideColors = colByTissue, cexRow = 0.42,
          density.info = "none", key.xlab = "log2FC", denscol = "white", key.title="",
          labCol = sub("_.HuGene.2_0.st.", "", colnames(lExprs$log2FC_blood), fixed = TRUE),
          lhei = c(1, 5), main = "209 differentially expressed genes\nin mucoid vs serous (|log2(FC)| > 1)", useRaster = TRUE)
dev.off()    

lExprs$log2FC_mvss <- lExprs$log2[,metaData$Tissue != "blood"]
lExprs$log2FC_mvss <- lExprs$log2FC_mvss - rowMeans(lExprs$log2FC_mvss)

pdf(file = "Rplot/brainArray_tm/heatmap_sereous_vs_mucoid_DGE_p0.01_lfc1_neutral.pdf", width = 8.3, height = 11.7)
heatmap.2(as.matrix(lExprs$log2FC_mvss[degl2,]), Colv = TRUE, dendrogram = "both", breaks = mybreaks, col = myGradient,
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          trace = "none", ColSideColors = colByTissue[metaData$Tissue != "blood"], cexRow = 0.42,
          density.info = "density", key.xlab = "log2FC", denscol = "white", key.title="",
          labCol = sub("_.HuGene.2_0.st.", "", colnames(lExprs$log2FC_blood)[metaData$Tissue != "blood"], fixed = TRUE),
          lhei = c(1, 5), main = "209 differentially expressed genes\nin mucoid vs serous (|log2(FC)| > 1)", useRaster = TRUE)
dev.off()    



# --------------

fcTables <- list(
    serous_vs_blood = topTable(fit2, coef = 1, adjust = "BH", number = 1000000),
    mucoid_vs_blood = topTable(fit2, coef = 2, adjust = "BH", number = 1000000),
    mucoid_vs_serous = topTable(fit2, coef = 3, adjust = "BH", number = 1000000)
)
lapply(fcTables, dim)

# volcanoes ---------------
lapply(fcTables, function(x) any(x$P.Value == 0)) # FALSE
fcTablesMessedUp <- lapply(fcTables, function(x)
    mutate(x, "log10.p.value" = -log10(adj.P.Val))
)

fcTablesMessedUp <- lapply(fcTablesMessedUp, function(myFcTable, alFCt = 1, pvalt = 0.01) {
    mutate(myFcTable, "significant" = ifelse(abs(logFC) >= alFCt & adj.P.Val <= pvalt, TRUE, FALSE))
})

myVolcanoes <- list(
    serous_vs_blood = ggplot(fcTablesMessedUp$serous_vs_blood, aes(x = logFC, y = log10.p.value, color = significant)) + geom_point(size = 1) +
        ggtitle("Serous - Blood") + xlab("log2(FC)") + ylab("-log10(P-value)") + theme_bw() +
        xlim(-5, 5) + ylim(0, 18) +  theme(legend.position = "none") + scale_color_manual(values = c("#808080", "#000000")) +
        annotate("text", x =  c(-4, 4), y = 18, label = c("858", "1138"), size = 4),
    
    mucoid_vs_blood = ggplot(fcTablesMessedUp$mucoid_vs_blood, aes(x = logFC, y = log10.p.value, color = significant)) + geom_point(size = 1) +
        ggtitle("Mucoid - Blood") + xlab("log2(FC)") + ylab("-log10(P-value)") + theme_bw() +
        xlim(-5, 5) + ylim(0, 18) +  theme(legend.position = "none") + scale_color_manual(values = c("#808080", "#000000")) + 
        annotate("text", x =  c(-4, 4), y = 18, label = c("967", "991"), size = 4),
    
    mucoid_vs_serous = ggplot(fcTablesMessedUp$mucoid_vs_serous, aes(x = logFC, y = log10.p.value, color = significant)) + geom_point(size = 1) +
        ggtitle("Mucoid - Serous") + xlab("log2(FC)") + ylab("-log10(P-value)") + theme_bw() +
        xlim(-5, 5) + ylim(0, 18) +  theme(legend.position = "none") + scale_color_manual(values = c("#808080", "#000000")) +
        annotate("text", x =  c(-4, 4), y = 18, label = c("144", "65"), size = 4)
)

pdf(file = "Rplot/brainArray_tm/volcanoes_adj_coloured.pdf", width = 9, height = 3)
plot_grid(myVolcanoes[[1]], myVolcanoes[[2]], myVolcanoes[[3]], "labels" = LETTERS[1:3], "ncol" = 3, label_size = 20)
dev.off()

system("firefox Rplot/brainArray_tm/volcanoes_adj_coloured.pdf &")


# exporting lists ---------------

for (i in 1:3){
    write.table(fcTables[[i]],
                file = paste0("ProtGenes_ba19_fc_table_", names(fcTables)[i],".txt"),
                sep = "\t", quote = FALSE, row.names = TRUE
    )
}

writeGeneListsFor <- function(myTable, file.prefix, abs.log2.FC = 1.0, p.val = 0.01, ref = "hugo") {
    upTable <- subset(myTable, logFC >= abs.log2.FC & adj.P.Val <= p.val)
    upTable <- row.names(upTable)
    downTable <- subset(myTable, logFC <= -abs.log2.FC & adj.P.Val <= p.val)
    downTable <- row.names(downTable)
    write.table(upTable, file = paste0(file.prefix, "_up.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(downTable, file = paste0(file.prefix, "_down.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

writeGeneListsFor(fcTables$serous_vs_blood, "geneLists/ProtGenes_ba19_serous_vs_blood_adj")
writeGeneListsFor(fcTables$mucoid_vs_blood, "geneLists/ProtGenes_ba19_mucoid_vs_blood_adj")
writeGeneListsFor(fcTables$mucoid_vs_serous, "geneLists/ProtGenes_ba19_mucoid_vs_serous_adj")


write.table(cbind(lExprs$anno[,c(1,2)], lExprs$data), file = "geneLists/ProtGenes_ba19_all_exprs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#plot(fcTables$serous_vs_blood$P.Value, fcTables$serous_vs_blood$adj.P.Val, pch = 20)


# qPCR -------------
qPCR <- read.table("gene_qPCR.txt", sep="\t", header = FALSE, stringsAsFactors = FALSE)
qPCR <- toupper(qPCR[,1])

qPCRvalid <- subset(cbind(lExprs$anno[,c(1,2)], lExprs$data), geneName %in% qPCR)
qPCRvalid <- qPCRvalid[order(qPCRvalid$geneName),]

magicalyFormatData <- function(x) {
    tdt <- data.frame(sample = rownames(metaData),
                      tissue = as.character(metaData$Tissue),
                      expr = t(qPCRvalid[x, 3:ncol(qPCRvalid)]),
                      stringsAsFactors = FALSE)
    row.names(tdt) <- NULL
    colnames(tdt) <- c("sample", "tissue", "expr")
    return(tdt)
}

nm <- 1:nrow(qPCRvalid)
tdtlist <- lapply(nm, magicalyFormatData)

doNicePlotWithggplot2 <- function(x) {
    naziplot <- ggplot(tdtlist[[x]], aes(x = tissue, y = expr, fill = tissue)) + geom_boxplot() +
        ggtitle(qPCRvalid[x,"geneName"]) +
        ylab("expression level (A.U.)") + theme_bw() +
        scale_fill_manual(values = c("green", "red", "purple")) + coord_cartesian(ylim = c(0, ceiling(max(tdtlist[[x]]$expr)/10)*10))
    return(naziplot)
}

myGraphs <- lapply(nm, doNicePlotWithggplot2)

pdf(file = "Rplot/brainArray_tm/valid_qPCR_boxplot.pdf", width = 4, height = 3)
for (i in nm) {
    print(myGraphs[[i]])
}
dev.off()

system("firefox Rplot/brainArray_tm/valid_qPCR_boxplot.pdf &")


doMyRatios <- function(x){
    result <- qPCRvalid[, 2*x+2]/qPCRvalid[, 2*x+1]
}
ratioValid <- do.call(cbind, lapply(1:11, doMyRatios))
ratioValid <- cbind(ratioValid, qPCRvalid[,25]/qPCRvalid[,23])
geneName <- qPCRvalid[order(apply(ratioValid, 1, median)), "geneName"]
ratioValid <- ratioValid[order(apply(ratioValid, 1, median)),]

library(reshape2)

rownames(ratioValid) <- geneName
mratioValid <- melt(ratioValid)
colnames(mratioValid) <- c("geneID", "patientID", "expr")
mratioValid$expr <- log2(mratioValid$expr)

pdf(file = "Rplot/brainArray_tm/valid_qPCR_big_boxplot.pdf", width = 5, height = 10)
ggplot(mratioValid, aes(x = factor(geneID), y = expr)) + geom_boxplot() + theme_bw() + coord_flip() +
    scale_x_discrete(labels = geneName) + geom_hline(yintercept = 0) + labs(y = "log2(FC)", x = "gene") +
    theme(axis.text.y  = element_text(size = 10))
dev.off()

system("firefox Rplot/brainArray_tm/valid_qPCR_big_boxplot.pdf &")

library(readr)
qPCRdata <- read_tsv("qPCR_data.txt") %>% as.data.frame
qPCRdata[,1] <- toupper(qPCRdata[,1])
qPCRl <- list(metadata = qPCRdata[,1:2], data = as.matrix(qPCRdata[,3:ncol(qPCRdata)]))
colnames(qPCRl$data) <- NULL
rownames(qPCRl$data) <- qPCRl$metadata[,1]
qPCRl$data <- melt(qPCRl$data)
colnames(qPCRl$data) <- c("geneID", "patientID", "expr")
qPCRl$data$expr <- log2(qPCRl$data$expr)

qPCRl$data <- cbind(qPCRl$data, "experiment" = "RTqPCR")
mratioValid <- cbind(mratioValid, "experiment" = "microarray")
bxpdata <- rbind(qPCRl$data, mratioValid)

pdf(file = "Rplot/brainArray_tm/valid_qPCR_big_boxplot_compare.pdf", width = 7, height = 11.4)
ggplot(bxpdata, aes(x = factor(geneID, levels = geneName, ordered = TRUE), y = expr, fill = experiment)) + geom_boxplot(outlier.size = 1) + theme_bw() + coord_flip() +
    geom_hline(yintercept = 0) + labs(y = "log2(FC)", x = "gene") + ggtitle("Supplementary figure X: Microarray\nvalidation using RTqPCR") +
    theme(axis.text.y  = element_text(size = 10), legend.position = "top") + scale_fill_manual(values = c("#FFFFFF", "#808080"), breaks = c("microarray","RTqPCR")) +
    scale_y_continuous(breaks = seq(-10, 15, 5))
dev.off()

system("firefox Rplot/brainArray_tm/valid_qPCR_big_boxplot_compare.pdf &")

subsetHV <- c("NRP2",
              "NRP1",
              "JUN",
              "VEGFA",
              "KIT",
              "SPHK1",
              "PGF",
              "MECOM",
              "IL6",
              "FIGF",
              "KDR",
              "CAV1",
              "NOS2",
              "FLT4",
              "SHC2",
              "VEGFC",
              "KRT18",
              "PTGS2",
              "SMAD2",
              "IL1A",
              "SMAD3"
)

bxpdataHV <- subset(bxpdata, bxpdata$geneID %in% subsetHV)
bxpdataLV <- subset(bxpdata, !(bxpdata$geneID %in% subsetHV))

doStatTestsForGene <- function(geneName2) {
    valuesArray <- ratioValid[geneName2, , drop = TRUE] %>% log2
    tArray <- t.test(valuesArray, mu = 0)$p.value
    valuesPCR <- qPCRdata[which(qPCRdata$Assay == geneName2), 3:ncol(qPCRdata)] %>% as.matrix %>% as.vector %>% log2
    tPCR <- t.test(valuesPCR, mu = 0)$p.value
    return(cbind(array = tArray, qPCR = tPCR))
}

statTestResults <- do.call(rbind, lapply(unique(qPCRdata$Assay), doStatTestsForGene))
rownames(statTestResults) <- unique(qPCRdata$Assay)

pdf(file = "Rplot/brainArray_tm/valid_qPCR_big_boxplot_compare_2p.pdf", width = 8, height = 10)
grid.arrange(
    ggplot(bxpdataHV, aes(x = factor(geneID, levels = geneName, ordered = TRUE), y = expr, fill = experiment)) + geom_boxplot(outlier.size = 1) + theme_bw() + coord_flip() +
        geom_hline(yintercept = 0) + labs(y = "log2(FC)", x = "gene") + ggtitle("Supplementary figure 2: Microarray\nvalidation using RTqPCR") +
        theme(axis.text.y  = element_text(size = 10), legend.position = "top") + scale_fill_manual(values = c("#FFFFFF", "#808080"), breaks = c("microarray","RTqPCR")) +
        scale_y_continuous(breaks = seq(-10, 15, 5)) +
        annotate("text", x = c(0.7, 1.1, 3.7, 5.7, 7.7, 8.7, 9.7, 11.7, 12.7, 13.1, 13.7, 14.1, 14.7, 15.1, 15.7, 16.1, 16.7, 17.1, 17.7, 18.1, 18.7, 19.1, 19.7, 20.1),
                 y = -9, label = "*", size = 6, col = "darkgray"),
    ggplot(bxpdataLV, aes(x = factor(geneID, levels = geneName, ordered = TRUE), y = expr, fill = experiment)) + geom_boxplot(outlier.size = 1) + theme_bw() + coord_flip() +
        geom_hline(yintercept = 0) + labs(y = "log2(FC)", x = "gene") +
        theme(axis.text.y  = element_text(size = 10), legend.position = "none") + scale_fill_manual(values = c("#FFFFFF", "#808080"), breaks = c("microarray","RTqPCR")) +
        annotate("text", x = c(0.65, 1.05, 1.65, 2.05, 2.65, 3.05, 3.65, 4.65, 5.65, 6.65, 7.65, 8.65, 9.65, 10.65, 12.65, 16.05, 19.05, 19.65, 20.05, 20.65, 21.05, 22.05, 22.65, 23.05, 24.05,
                               24.65, 25.05, 25.65, 26.05, 27.05, 27.65, 28.05, 28.65, 29.05, 29.65, 30.05, 30.65, 31.05, 31.65, 32.05),
                 y = -4, label = "*", size = 6, col = "darkgray"),
    ncol = 2
)
dev.off()

system("firefox Rplot/brainArray_tm/valid_qPCR_big_boxplot_compare_2p.pdf &")




