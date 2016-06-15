# running localy for some reasons

# convert text cdf into binary cdf
# cdf download from:
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/19.0.0/enst.asp
# http://mbni.org/customcdf/19.0.0/enst.download/hugene20st_Hs_ENST_19.0.0.zip
library(affxparser)
convertCdf("aromaBrainArray19/annotationData/chipTypes/HuGene-2_0-st/hugene20st_Hs_ENST.cdf", "aromaBrainArray19/annotationData/chipTypes/HuGene-2_0-st/HuGene-2_0-st.cdf")


#setwd("/groups2/joshi_grp/guillaume/otherProject/earMicroarrays/aromaBrainArray19/")
setwd("C:/DevaillyDell/aromaBrainArray19")
library("aroma.affymetrix")

# RMA with aroma --------------------
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

chipType <- "HuGene-2_0-st"
cdf <- AffymetrixCdfFile$byChipType(chipType)
print(cdf)

cs <- AffymetrixCelSet$byName("earArrays", cdf = cdf)
print(cs)

bc <- RmaBackgroundCorrection(cs)
csBC <- process(bc, verbose = verbose)

qn <- QuantileNormalization(csBC, typesToUpdate="pm")
print(qn)
csN <- process(qn, verbose = verbose)

plm <- RmaPlm(csN)
print(plm)
fit(plm, verbose=verbose)

qam <- QualityAssessmentModel(plm)

pdf(file = "qc_NUSE_RLE_ba.pdf", width = 11, height = 8)
par(mar = c(12,4,4,2))
layout(cbind(1,2))
plotNuse(qam)
plotRle(qam)
dev.off()

ces <- getChipEffectSet(plm)
gExprs <- extractDataFrame(ces, units=NULL, addNames=TRUE)

lExprs <- list("anno" = gExprs[,1:5], "data" = gExprs[,6:28])
save(lExprs, file="rma_data_branArrayENST.RData")

