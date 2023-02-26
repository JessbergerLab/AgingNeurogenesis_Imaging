rm(list = ls())
setwd("~/Desktop/in vivo imaging/RSTUDIO/ROI/")

library(HippoLinTools)
library(beeswarm)
require(dae)
require(igraph)
require(lattice)
require(pheatmap)
require(grid)
require(RColorBrewer)
require(ggfortify)
require(ade4)
require(plotrix)
require(ggplot2)
library(knitr)
library(Hmisc)
library(ggbeeswarm)

myConfidenceThreshold <- 0.95
lateDeathThreshold <- 7
onlyRroots <- TRUE
thresholdDivisionWindow <- 6

RoiDirectoryPathGli1 <- "~/Desktop/In vivo imaging/RStudio/Analysis/Aging/ROI/"
CorrespTimePointDPIPath <- "~/Desktop/in vivo imaging/RStudio/Analysis/Aging/Correspondence/"
decodedDataGli1 <- decodeROIData(ROIpath = RoiDirectoryPathGli1, DPITimePointPath = CorrespTimePointDPIPath, fillUnchangedTimePoint = T)

plotLineage2(decodedData = decodedDataGli1, saveGraph = F, myTime = "generation")
plotLineage2(decodedData = decodedDataGli1, saveGraph = F, myTime = "DPI")

#Or
plotLineage2(decodedData = decodedDataGli1, saveGraph = F, myTime = "generation", displayLabels = F, myFileName = "Rlineages.pdf")
plotLineage2(decodedData = decodedDataGli1, saveGraph = T, myTime = "DPI", displayLabels = F, myFileName = "Rlineages.pdf")

decodedDataGli1 <- neuronAssignment(decodedData = decodedDataGli1, confidenceThreshold = myConfidenceThreshold)
write.table(x = decodedDataGli1, file = "Gli1_to_test.txt", quote = FALSE, sep = "\t")
