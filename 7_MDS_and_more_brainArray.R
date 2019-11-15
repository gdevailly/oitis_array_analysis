setwd("/groups2/joshi_grp/guillaume/otherProject/earMicroarrays/")
load("rma_data_branArrayENST.RData")

metaData <- read.table("metaData.txt", header = TRUE, row.names = 1)
metaData[23, 1] <- 2
metaData$Tissue <- factor(metaData$Tissue)
levels(metaData$Tissue) <- c("blood", "serous", "mucoid")
row.names(metaData) <- sub(".CEL", "", row.names(metaData))
lExprs$data <- lExprs$data[, row.names(metaData)]

# logFC transformation
any(rowMeans(lExprs$data) == 0)
lExprs$logFC_all <- log2(lExprs$data / rowMeans(lExprs$data))

# Classical MDS ------------------------

d <- dist(t(lExprs$logFC_all)) # euclidean distances between the rows
fit <- cmdscale(d, eig = TRUE, k = 2) # k is the number of dim
fit # view results

# by tissue
colByTissue <- factor(metaData$Tissue)
levels(colByTissue) <- c("green", "red", "purple")
colByTissue <- as.character(colByTissue)

pdf(file = "Rplot/brainArray/MDS_all.pdf", width = 5, height = 5)
plot(fit$points[,1], fit$points[,2], xlab="logFC dim 1", ylab="logFC dim 2", main="MDS plot", pch = 20, xlim = c(-200, 240), col = colByTissue)
text(fit$points[,1], fit$points[,2], labels = sub("_(HuGene-2_0-st)", "", row.names(fit$points), fixed = TRUE), cex = 1, pos = 4, col = colByTissue)
legend("topright", col = levels(factor(colByTissue)), legend = c("blood", "mucoid", "serous"), pch = 20)
dev.off()

system("firefox Rplot/brainArray/MDS_all.pdf &")


# by patient
colByPatients <- factor(metaData$Patients)
levels(colByPatients) <- rainbow(11)
colByPatients <- as.character(colByPatients)

pdf(file = "Rplot/brainArray/MDS_all_by_patient.pdf", width = 5, height = 5)
plot(fit$points[,1], fit$points[,2], xlab="logFC dim 1", ylab="logFC dim 2", main="MDS plot", pch = 20, xlim = c(-200, 240), col = colByPatients)
text(fit$points[,1], fit$points[,2], labels = sub("_(HuGene-2_0-st)", "", row.names(fit$points), fixed = TRUE), cex = 1, pos = 4, col = colByPatients)
for(i in levels(factor(metaData$Patients))) {
    temp <- fit$points[metaData$Patients == i,]
    colP <- colByPatients[metaData$Patients == i]
    segments(temp[1,1], temp[1,2], temp[2,1], temp[2,2], col = colP[1], lty = 2)
    if (i ==58) segments(temp[1,1], temp[1,2], temp[3,1], temp[3,2], col = colP[1], lty = 2)
}
dev.off()

pdf(file = "Rplot/brainArray/new_MDS_all_patientLinks.pdf", width = 5, height = 5)
plot(fit$points[,1], fit$points[,2], xlab="dimension 1", ylab="dimension 2", main="Multidimensional\nscaling of the microarray data", pch = 20, xlim = c(-200, 240), col = colByTissue)
text(fit$points[,1], fit$points[,2], labels = sub("_(HuGene-2_0-st)", "", row.names(fit$points), fixed = TRUE), cex = 1, pos = 4, col = colByTissue)
legend("topright", col = levels(factor(colByTissue)), legend = c("blood", "mucoid", "serous"), pch = 20)
for(i in levels(factor(metaData$Patients))) {
    temp <- fit$points[metaData$Patients == i,]
    colP <- colByPatients[metaData$Patients == i]
    segments(temp[1,1], temp[1,2], temp[2,1], temp[2,2], col = "gray", lty = 2)
    if (i ==58) segments(temp[1,1], temp[1,2], temp[3,1], temp[3,2], col = "gray", lty = 2)
}
dev.off()

png(file = "Rplot/brainArray/new_MDS_all_patientLinks.png", width = 5, height = 5, units = "in", res = 300)
plot(fit$points[,1], fit$points[,2], xlab="dimension 1", ylab="dimension 2", main="Multidimensional\nscaling of the microarray data", pch = 20, xlim = c(-200, 240), col = colByTissue)
text(fit$points[,1], fit$points[,2], labels = sub("_(HuGene-2_0-st)", "", row.names(fit$points), fixed = TRUE), cex = 1, pos = 4, col = colByTissue)
legend("topright", col = levels(factor(colByTissue)), legend = c("blood", "mucoid", "serous"), pch = 20)
for(i in levels(factor(metaData$Patients))) {
    temp <- fit$points[metaData$Patients == i,]
    colP <- colByPatients[metaData$Patients == i]
    segments(temp[1,1], temp[1,2], temp[2,1], temp[2,2], col = "gray", lty = 2)
    if (i ==58) segments(temp[1,1], temp[1,2], temp[3,1], temp[3,2], col = "gray", lty = 2)
}
dev.off()


system("firefox Rplot/brainArray/MDS_all_patientLinks.pdf &")

# by batch
colByBatch <- factor(metaData$Batch)
levels(colByBatch) <- c("blue", "red")
colByBatch <- as.character(colByBatch)

pdf(file = "Rplot/brainArray/MDS_all_byBatch.pdf", width = 5, height = 5)
plot(fit$points[,1], fit$points[,2], xlab="logFC dim 1", ylab="logFC dim 2", main="MDS plot", pch = 20, xlim = c(-200, 240), col = colByBatch)
text(fit$points[,1], fit$points[,2], labels = sub("_(HuGene-2_0-st)", "", row.names(fit$points), fixed = TRUE), cex = 1, pos = 4, col = colByBatch)
legend("topright", col = levels(factor(colByBatch)), legend = c("Batch 1", "Batch 2"), pch = 20)
dev.off()

system("firefox Rplot/brainArray/MDS_all_byBatch.pdf &")

# predictive MDS -----------------

d <- dist(t(lExprs$logFC_all[,metaData$Tissue == "blood"])) # euclidean distances between the rows
fit <- cmdscale(d, eig = TRUE, k = 2) # k is the number of dim
fit # view results

colByBatch <- factor(metaData$Batch[metaData$Tissue == "blood"])
levels(colByBatch) <- c("purple", "red")
colByBatch <- as.character(colByBatch)

pdf(file = "Rplot/brainArray/predictiveMDS_all_byBatch.pdf", width = 5, height = 5)
plot(fit$points[,1], fit$points[,2], xlab="logFC dim 1", ylab="logFC dim 2", main="MDS plot", pch = 20, xlim = c(-100, 100), col = colByBatch)
text(fit$points[,1], fit$points[,2], labels = sub("_(HuGene-2_0-st)", "", row.names(fit$points), fixed = TRUE), cex = 1, pos = 4, col = colByBatch)
legend("topright", col = levels(factor(colByBatch)), legend = c("Batch 1, mucoid", "Batch 2, serous"), pch = 20)
dev.off()

system("firefox Rplot/brainArray/predictiveMDS_all_byBatch.pdf &")

