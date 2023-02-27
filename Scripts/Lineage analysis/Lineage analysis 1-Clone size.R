################################################################################################################################
# Removes everything from the workspace
rm(list = ls())
################################################################################################################################

#############################################
# Parameters
#############################################

confidenceThreshold <- 0.95 # 5% chance to be wrong when assigning an uncertain neuroblast/neuron as certain neuron depending on how long it has been observed
lateDeathThreshold <- 7 # if neuron dies after surviving more than 7 days it is considered late cell death (activity dependent)
onlyRroots <- TRUE # considers only the trees with a R root

#############################################
# Load library
#############################################

library(HippoLinTools)
library(beeswarm)
require(dae)
require(igraph)
require(lattice)
require(pheatmap)
require(grid)
require(RColorBrewer)
require(FactoMineR)
require(ggfortify)
require(ade4)
require(plotrix)
require(Hmisc)
require(ggplot2)

#############################################
# Load data, generates the different kinds of subsets of datatables
#############################################

myDate <- '080920' 
RoiDirectoryPath <- "../Data/ROI/"
CorrespTimePointDPIPath <- "../Data/Correspondence/"



# Convert the ROI code into a nice table
decodedData <- decodeROIData(ROIpath = RoiDirectoryPath, DPITimePointPath = CorrespTimePointDPIPath, fillUnchangedTimePoint = T )

# Make uncertain neurons observed for a long enough time are certain (not necessary if there is no uncertainty)
decodedData <- neuronAssignment(decodedData, confidenceThreshold = confidenceThreshold)

#############################################
## To export this table as a txt file that you can import in excel:
## x = which table to export
## file = name of the created file
#############################################

write.table(x = decodedData, file = paste("DecodedRoiData", myDate,".txt",sep=""), quote = FALSE, sep = "\t", row.names = F)

#############################################
## Transform the tree-structure table into timepoint-structured table, and also datasets with only first and/or last time point for every cell
#############################################

myLegend <-  c("R","NR","N","Cell_Death","Uncertain")

# Dataset with one line per time step and total number of cells of each types (Celltype 2 if myLegend <-  c("R","NR","N","Cell_Death","Uncertain"))
# Tp: time point 
dataTp <- transformTreeIntoTp(data2Transform = decodedData, namesTypes = myLegend)

# Cumulated number of cell death per lineage
dataTp$cumCellDeath <- unlist(with(dataTp, tapply(Cell_Death,INDEX = Clone_uniqueID,FUN = cumsum)))

# One row per cell, we take into account only the last time point the cell is seen
# CellType: 1, 2, 3 which has been named during recording
# CellType2: 6 means 'cell death'
# CellType3: final assignment of the cell type: 3 -> 4
# Rank: at which round of division when it is recorded, mother cell (root cell) is '1', first-genration progeny is '2'
# Rank NR: at which round of division when it is recorded (only consider NR cell), mother cell (root cell) is '1', first-genration progeny is '2'
dataMax <- transformFullTreeIntoFirstLastTree(decodedData,rows2Keep = c("last")) 

dataMax$UncertaintyType2 <- as.character(dataMax$UncertaintyType2)

# One row per cell, we take into account only the first time point the cell is seen
dataMin <- transformFullTreeIntoFirstLastTree(decodedData, rows2Keep = c("first")) 

# Every cell has 2 rows, one when it is first seen and when it is last seen, except if the cell is seen only 1 time point
dataMinMax <- rbind(dataMax, dataMin) 

#############################################
## Selecting only trees starting with R cell
#############################################

if(onlyRroots == T){
  
  treeID <- unique(decodedData$Clone_uniqueID)
  
  dataMaxR <- dataMax[dataMax$CellType2 == 1 & dataMax$UncertaintyType2 == 1,]
  RlineageNames <- unique(dataMaxR$Clone_uniqueID)
  
  dataMax <- dataMax[dataMax$Clone_uniqueID %in% RlineageNames,]
  dataMin <- dataMin[dataMin$Clone_uniqueID %in% RlineageNames,]
  dataMinMax <- dataMinMax[dataMinMax$Clone_uniqueID %in% RlineageNames,]
  
  decodedData <- decodedData[decodedData$Clone_uniqueID %in% RlineageNames,]
  dataTp <- dataTp[dataTp$Clone_uniqueID %in% RlineageNames,]
  
}

# Try
if(onlyRroots == T){
  
  treeID <- unique(decodedData$Clone_uniqueID)
  
  dataMaxR <- dataMax[dataMax$CellType2 == 1,]
  
  RlineageNames <- unique(dataMaxR$Clone_uniqueID)
  
  dataMax <- dataMax[dataMax$Clone_uniqueID %in% RlineageNames,]
  dataMin <- dataMin[dataMin$Clone_uniqueID %in% RlineageNames,]
  dataMinMax <- dataMinMax[dataMinMax$Clone_uniqueID %in% RlineageNames,]
  
  decodedData <- decodedData[decodedData$Clone_uniqueID %in% RlineageNames,]
  dataTp <- dataTp[dataTp$Clone_uniqueID %in% RlineageNames,]
  
}

#############################################
# to export this table as a txt file that you can import in excel:
# x = which table to export
# file = name of the created file
#############################################
# write.table(x = decodedData, file = paste("DecodedRoiData", myDate,'_',myDate2,".txt",sep=""), quote = FALSE, sep = "\t", row.names = F)

#############################################
## Number of successive divisions in each tree
#############################################

numberSuccessiveDivision <- with(dataMax, tapply(X = Rank, INDEX =Clone_uniqueID, FUN = max)) - 1 # Rank 1 = no successive divisions : explains the -1
quiescentClonesID <-  names(numberSuccessiveDivision)[numberSuccessiveDivision == 0]
nonQuiescentID <- as.character(names(numberSuccessiveDivision)[!names(numberSuccessiveDivision) %in% quiescentClonesID ])

# Plot lineages 6 per page
# source('./test_plotlineage_6plotsPerPage.R')
# plotLineage2(decodedData,saveGraph = T,displayLabels = F, myFileName = "Rlineages_6perPage_070717.pdf")
# 
# # number of cells per mouse
# with(dataMax, tapply(Spot_ID, INDEX = Mouse_ID, length))
# # number of clones per mouse
# unique(dataMax$Clone_uniqueID)
# 
# numberLineagePerMouse <- c(7,3,6,3,6,2,3,9,24)
# mean(numberLineagePerMouse)
# sd(numberLineagePerMouse)/sqrt(length(numberLineagePerMouse))


################################################################################################################################
# Cell type uncertainty (don't have to run when there is no uncertainty)
#############################################
# Either considering the cell type uncertainty at every time point (potentially several time steps per cell) or considering only the cell type 
# uncertainty at one time point for every cell, namely the last time point before division or death. Cell death is always certain in the dataset.
################################################################################################################################

#############################################
## Proportion of certain / uncertain types and transitions, here all cells, all type points
#############################################

# Two lines run together (with 'legend')
barplot(table(decodedData$UncertaintyType2, decodedData$CellType2), 
        main = "Proportion uncertain cell types \n all time points", names.arg = c("R", "NR","N","Death"), ylab = "# cells", beside = TRUE, col = c("black","gray","white"))
legend("topleft",legend = c("uncertain","semi-certain", "certain"),fill = c("black","gray","white"), cex =0.8,bty = 'n')

barplot(prop.table(table(decodedData$UncertaintyType2, decodedData$CellType2),2), 
        main = "Frequency uncertain cell types \n all time points", 
        names.arg = c("R", "NR","N","Death"), ylab = "Frequency of cells", beside = TRUE, col = c("black","gray","white"))

#############################################
## Proportion of certain / uncertain types and transitions, here all cells, only last time point per cell
#############################################

barplot(table(dataMax$UncertaintyType3, dataMax$CellType3), 
        main = paste("Proportion uncertain cell types \n only last time point per cell n=", 
        dim(dataMax)[1]," cells"), names.arg = c("R", "NR","N"), ylab = "# cells", beside = TRUE, col = c("black","gray","white"), cex.main = 0.8)
legend("topleft",legend = c("uncertain","semi-certain", "certain"),fill = c("black","gray","white"), cex =0.8,bty = 'n')

barplot(prop.table(table(dataMax$UncertaintyType3, dataMax$CellType3),2), 
        main = "Frequency uncertain cell types \n only last time point per cell", 
        names.arg = c("R", "NR","N"), ylab = "Frequency of cells", beside = TRUE, col = c("black","gray","white"), cex.main = 0.9)

tab <- prop.table(table(dataMax$UncertaintyType3, dataMax$CellType3),2)
colnames(tab) <- c("R", "NR","N")
print("% of certain (1) uncertain (0 or 0-5) types among the different cell types")
print(round(tab,2)*100)


################################################################################################################################
# Lineage uncertainty
################################################################################################################################

#############################################
## Remove the root cell of the lineage for computing, relationship uncertainty, because it is never known by construction
#############################################

dataMinWithoutRoot <- dataMin[dataMin$Cell_ID != "1" & dataMin$Cell_ID != "1-1",]

barplot((table(dataMinWithoutRoot$UncertaintyLineage)), col = "white", 
        main = paste("Proportions of uncertain transitions \n", dim(dataMinWithoutRoot)[1]," total transitions"), 
        ylab = "# transitions", cex.main = 0.9)
barplot(prop.table(table(dataMinWithoutRoot$UncertaintyLineage)), col = "white", 
        main = "Frequencies of uncertain transitions", ylab ="Transition frequency", cex.main = 0.9)

#############################################
## Number of divisions with 2 daughters, with more than 2 daughters
## If bug here it might be because there are no 0-5 unsure transitions (2 factors instead of 3). 
## If the case do not forget to change the colors in the graphs and the legends ("black", "white")
#############################################

numberDaughtersPerMother <- as.vector(table(dataMinWithoutRoot$Mother_uniqueID))
NumberDaughtersUncertainty <- as.matrix(table(dataMinWithoutRoot$Mother_uniqueID, dataMinWithoutRoot$UncertaintyLineage))
NumberDaughtersUncertainty <- as.data.frame(matrix(NumberDaughtersUncertainty, nrow = length(unique(dataMinWithoutRoot$Mother_uniqueID)), ncol = 3))
NumberDaughtersUncertainty$totalTransitions <- rowSums(NumberDaughtersUncertainty)

barplot(table(NumberDaughtersUncertainty$V3,NumberDaughtersUncertainty$totalTransitions), 
        main = paste("Proportions of divisions with \n2 or more daughters n= ",length(numberDaughtersPerMother),"divisions",sep =" "), 
        ylab = "# divisions", xlab = "number of daughters per mother", cex.main = 0.9, col = c("black","gray","white"), beside = F)
legend("topright", legend = c("2 uncertain","2 semi-certain", "2 certain"),fill = c("black","gray","white"), cex = 0.8, bty = "n")

propTab <- prop.table(table(NumberDaughtersUncertainty$V3,NumberDaughtersUncertainty$totalTransitions))
barplot(propTab, main = "Frequency of divisions \n with 2 or more daughters", 
        ylab = "Frequency of divisions", xlab = "number of daughters per mother", cex.main = 0.9, col= c("black","gray","white"))
legend("topright", legend = c("2 uncertain","2 semi-certain", "2 certain"),fill = c("black","gray","white"), cex = 0.8, bty = "n")


################################################################################################################################
# Attempt to cluster lineages (the results may not be very correct)
# Characteristics to cluster lineages: maximum number of cells, cell death rate, proportion of R cell divisions, number of divisions per lineage.
################################################################################################################################

# Build a table with characteristics of each lineage (Clone ID, number successive divisions, max number of cells, final number of cells, cell death rate, self-renewing time, number of R successive divisions):
# !!!!!!!!!!!!!!!!!be careful the order of lineages from dataTp and decodedData is different!!!!!!!!!!!!!!!!!!!!!!!!
# We use the order of decodedData!!

cloneFinalComposition <- cloneOutcome(dataTp) 
treeID <- unique(decodedData$Clone_uniqueID)
reorderDataTp <- match(treeID, cloneFinalComposition$Clone_uniqueID)

#############################################
## Maximum cellcnumber
#############################################

maxCellNumberLineages <- tapply(X = dataTp$TotalCellNumber, as.factor(dataTp$Clone_uniqueID),  max)
maxCellNumberLineages <- maxCellNumberLineages[reorderDataTp]

# Checks that the orders of lineages is the same everywhere:
# match(treeID, names(maxCellNumberLineages)) - c(1:length(cloneFinalComposition$Clone_uniqueID))

#############################################
## Number of successive divisions in each tree
#############################################

numberSuccessiveDivision <- with(dataMax, tapply(X = Rank, INDEX = Clone_uniqueID, FUN = max)) - 1 # rank 1 = no successive divisions : explains the -1
# To check whethe 'numberSuccessiveDivision' is correct

quiescentClonesID <-  names(numberSuccessiveDivision)[numberSuccessiveDivision == 0]

#############################################
## Final cell number
#############################################

finalCellNumberLineages <- cloneFinalComposition$TotalCellNumber
finalCellNumberLineages <- finalCellNumberLineages[reorderDataTp]

# match(treeID, cloneFinalComposition$Clone_uniqueID[reorderDataTp]) - c(1:length(cloneFinalComposition$Clone_uniqueID)) # test lineage order

#############################################
## Death rate per lineage
#############################################

deathRate <- c(1:length(treeID))

for(t in seq_along(treeID)){
  myTree <- dataMax[dataMax$Clone_uniqueID == treeID[t],]
  deathRate[t] <- getDeathRate(myTree)
}

# Generate the deathRate with lineage names: 
# deathRateLineage <- cbind(treeID, deathRate)

#############################################
## Number of R successive divsisions
#############################################

dataMaxProg <- dataMax[dataMax$Cell_uniqueID %in% dataMax$Mother_uniqueID,]
dataMaxRprog <- dataMaxProg[dataMaxProg$CellType2 == 1 & dataMaxProg$UncertaintyType2 == 1,]

maxRRankAll <- c(1:length(treeID)) # Vector containing every lineages
nonRLineageID <- match(treeID[!treeID %in% dataMaxR$Clone_uniqueID],treeID)
RlineageID <- match(treeID[treeID %in% dataMaxR$Clone_uniqueID],treeID)

#############################################
## Maximum rank of a R cell in every lineage which contains at least one R cell
#############################################

maxRRankAll[nonRLineageID] <- 0
maxRRank <- with(dataMaxRprog, tapply(Rank, as.factor(Clone_uniqueID), max))  
maxRRank <- with(dataMax, tapply(Rank, as.factor(Clone_uniqueID), max))  
maxRRank[names(maxRRank) %in% quiescentClonesID ] <- 0 # Quiescent cells do not have any R divisions

maxRRankAll[RlineageID] <- maxRRank

#############################################
## Total number of divisions in the tree
#############################################

numberMothers <- function(lineageMothers){
  mothers <- unique(lineageMothers[2:length(lineageMothers)])
  number <- length(mothers) # starts at 2 because the mother of the root is unknown
  number 
}

numberMotherPerLineage <- with(dataMax, tapply(Mother_uniqueID, INDEX = Clone_uniqueID, FUN = numberMothers))
numberMotherPerLineage[numberSuccessiveDivision == 0] <- 0 # Quiescent clones do not have mothers because there is no division

#############################################
## Proportions of R divisions in the lineage 
#############################################

numberMotherPerLineage <- vector("numeric", length = length(treeID))
motherTypePerLineage <- matrix(nrow = length(treeID), ncol =3)

for (i in seq_along(treeID)){
  myTree <- dataMax[dataMax$Clone_uniqueID == treeID [i],]
  mothers <- unique(myTree$Mother_uniqueID[2:length(myTree$Mother_uniqueID)])
  numberMotherPerLineage[i] <- length(mothers)
  
  motherTypes <- myTree$CellType2[myTree$Cell_uniqueID %in% mothers & myTree$UncertaintyType2 == 1]
  motherTypes <- factor(motherTypes,levels = c("1","2"))
  motherTypes <- relevel(motherTypes, ref = "1")
  motherTypePerLineage[i,] <- c(table(motherTypes), sum(table(motherTypes)))
}

proportionRmother <- motherTypePerLineage[,1]/motherTypePerLineage[,3]
proportionRmother[is.na(proportionRmother)] <- 0 # Cells with uncertain types are omitted

#############################################
## Creates the table of lineage characteristics + makes sure all data are between 0 and 1
#############################################

dataLineages <- data.frame(treeID, maxCellNumberLineages, finalCellNumberLineages, 
                           numberSuccessiveDivision, deathRate, numberRsuccDiv = maxRRankAll, numberMotherPerLineage, proportionRmother)
# numberSuccessiveDivision is not correct, use the number from 'MaxRRank' which will be generated below
maxPerCol <- t(matrix(rep(apply(dataLineages[,c(2,5,7,8)], MARGIN = 2, FUN = max), length= nrow(dataLineages)*4), nrow = 4, ncol = nrow(dataLineages)))
dataLineagesNormalized <- dataLineages[,c(2,5,7,8)]/ maxPerCol
myDistance <- dist(dataLineagesNormalized, method = "euclidean")
clusters <- hclust(myDistance, method="ward.D") #ward.D

plot(clusters, main = "Hierachical clustering of lineages")

clusterCut <- cutree(clusters, 4)

print("Number of lineages per cluster")
print(table(clusterCut))
dataLineages$clusters <- as.character(clusterCut)

res.pca  <- prcomp(dataLineagesNormalized, scale. = T, center = TRUE)

autoplot(res.pca, data = dataLineages, colour = 'clusters', label = TRUE, loadings = TRUE, loadings.colour = 'yellow',
         loadings.label = TRUE, loadings.label.size = 5)

# Kmeans
# myKmeans <- kmeans(myDistance, centers = 5)
# kmeanCluster <- myKmeans$cluster
# mycol <- c("pink","green","blue","purple","cyan")
# autoplot(res.pca, data = dataLineages, colour = mycol[kmeanCluster], label = TRUE, loadings = TRUE, loadings.colour = 'yellow',
#          loadings.label = TRUE, loadings.label.size = 5, main = "Kmean clustering result")

# plot(dataLineages$maxCellNumberLineages ~ dataLineages$deathRate, pch = 16, cex = 0.8)
boxplot(dataLineages$maxCellNumberLineages ~ clusterCut, beside = TRUE, main = "Max cell number", xlab = "cluster number")
beeswarm(dataLineages$maxCellNumberLineages ~ clusterCut, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

boxplot(dataLineages$finalCellNumberLineages ~ clusterCut,  main = "Final cell number", xlab = "cluster number")
beeswarm(dataLineages$finalCellNumberLineages ~ clusterCut, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

boxplot(dataLineages$deathRate ~ clusterCut, main = "Death rate", xlab = "cluster number")
beeswarm(dataLineages$deathRate ~ clusterCut, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

boxplot(dataLineages$numberSuccessiveDivision ~ clusterCut, beside = TRUE,  main = "Number successive divisions", xlab = "cluster number")
beeswarm(dataLineages$numberSuccessiveDivision ~ clusterCut, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

boxplot(dataLineages$numberRsuccDiv ~ clusterCut, beside = TRUE,  main = "Number R successive divisions", xlab = "cluster number")
beeswarm(dataLineages$numberRsuccDiv ~ clusterCut,  add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

boxplot(dataLineages$proportionRmother ~ clusterCut, beside = TRUE,  main = "Proportion R divisions", xlab = "cluster number")
beeswarm(dataLineages$proportionRmother ~ clusterCut, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

dataLineages$mouseID <- as.character(dataLineages$treeID)
splitLineageID <- strsplit(as.character(dataLineages$treeID), split ="\\.") 
for(i in seq_along(splitLineageID)){
  dataLineages$mouseID[i] <-   splitLineageID[[i]][1]
}

barplot(with(dataLineages, table(mouseID, clusters)), col= rainbow(length(unique(dataLineages$mouseID))), xlab = "Cluster number", 
        ylab = "Number of lineages", beside =T, main = "Distribution of Mouse ID \nin the different clusters")
legend("top", legend = unique(dataLineages$mouseID), fill = rainbow(length(unique(dataLineages$mouseID))), bty = "n", ncol = 2, cex = 0.8)


################################################################################################################################
# Final clone size
################################################################################################################################

#############################################
## Final and max number of cells per clone (mean)
#############################################
cloneFinalComposition <- cloneOutcome(dataTp) # Composition of the clone at the last observed time point 

# Proportion active vs quiescent clones: get only lineages with R cells
Rlineages <- dataMax$Clone_uniqueID[RlineageID]

numberRClones <- length(Rlineages)
numberQuiescentClones <- length(quiescentClonesID)
numberActiveClones <- numberRClones - numberQuiescentClones 
myplot <- barplot(c(numberQuiescentClones, numberActiveClones)/numberRClones*100, names.arg = c("quiescent","active"), ylab = "Proportion of R lineages")
text(myplot, 5, labels =  c(numberQuiescentClones, numberActiveClones))

#############################################
## Number of successive divisions in each clone
#############################################

# Quiescent cells have rank 1 (first generation), but division 0, hence the - 1
numberSuccessiveDivision <- with(dataMax, tapply(X = Rank, INDEX = Clone_uniqueID, FUN = max)) -1 
# Double check whether this value is correct

reorderLineages <- match(cloneFinalComposition$Clone_uniqueID, names(numberSuccessiveDivision)) 
numberSuccessiveDivision <- numberSuccessiveDivision[reorderLineages]
cloneFinalComposition$numberSuccessiveDivision <- numberSuccessiveDivision

barplot(table(cloneFinalComposition$numberSuccessiveDivision), ylab = "# lineages", xlab = "# successive divisions in each clone")
myCor <- cor(cloneFinalComposition$numberSuccessiveDivision, cloneFinalComposition$Imaging_time_DPI)

plot(cloneFinalComposition$numberSuccessiveDivision ~ cloneFinalComposition$Imaging_time_DPI, xlab = "Imaging time [DPI]", ylab ="# successive divisions in each clone", 
     pch = 16, cex =0.8, main = paste("No relation between #successive divisions \n and imaging length cor = ",round(myCor,digits = 2)), cex.main = 0.8)

#############################################
## Clone final size
#############################################

cloneFinalComposition  <- cloneFinalComposition[!(cloneFinalComposition$Clone_uniqueID %in% quiescentClonesID),] # Remove the quiescent clones
myCor <- cor(cloneFinalComposition$TotalCellNumber , cloneFinalComposition$Imaging_time_DPI)

hist(cloneFinalComposition$TotalCellNumber, col = "white", ylab = "number of clones", xlab = "Number of cells at the final time point", cex.axis = 0.8, 
     main = paste("Distribution of final clone size n= ", dim(cloneFinalComposition)[1],"\n mean = ",round(mean(cloneFinalComposition$TotalCellNumber),2),"sd = ", 
     round(sd(cloneFinalComposition$TotalCellNumber),2), " quiescent clones excluded"), cex.main = 0.9, breaks = seq(0,max(cloneFinalComposition$TotalCellNumber,1)))
abline(v = mean(cloneFinalComposition$TotalCellNumber), col = "red", lwd = 2)

# par(mfrow = c(1,2))
beeswarm(cloneFinalComposition$TotalCellNumber, ylab = "Final number of cells per clone")
bxplot(cloneFinalComposition$TotalCellNumber, add = T, axis.lty = 1)

barPlotBeeswarm(cloneFinalComposition$TotalCellNumber, ylab = "Final number of cells per clone", cex.main = 0.8)

# par(mfrow = c(1,1))
meanCloneSize <- round(mean(cloneFinalComposition$TotalCellNumber),2)
semCloneSize <- round(sd(cloneFinalComposition$TotalCellNumber)/sqrt(length(cloneFinalComposition$TotalCellNumber)),2)

plot(cloneFinalComposition$TotalCellNumber ~ cloneFinalComposition$Imaging_time_DPI, xlab = "Imaging time [DPI]", ylab ="Final cell number in clones", pch = 16, cex =0.8, 
     main = paste("No relation between final clone size \n and imaging length cor = ",round(myCor,digits = 2)), cex.main = 0.8)


################################################################################################################################
# Maximun clone size
################################################################################################################################

maxCellNumberLineages <- tapply(X = dataTp$TotalCellNumber, as.factor(dataTp$Clone_uniqueID),  max)
maxCellNumberLineages

# hist(maxCellNumberLineages, xlab = "maximum number of cells in a lineage", 
# ylab = "Number of lineages", main = "Distribution of maximun number of cells per lineage", cex.main = 0.9)
# abline(v = mean(maxCellNumberLineages), col = "red", lwd = 2)

beeswarm(maxCellNumberLineages, ylab = "maximum number of cells in a lineage", xlab = "Number of lineages", main = "Distribution of maximun number of cells per lineage", cex.main = 0.9)
bxplot(maxCellNumberLineages, add = T, col = "red")

dataMaxR <- dataMax[dataMax$CellType2 == 1 & dataMax$UncertaintyType2 == 1 & !dataMax$Clone_uniqueID %in% quiescentClonesID ,]
# Try:
# dataMaxR <- dataMax[dataMax$CellType2 == 1 & !dataMax$Clone_uniqueID %in% quiescentClonesID ,]

#############################################
## Number of R successive divsision, vs maximum number of cells in the lineage
#############################################

# Maximum rank of a R cell in every lineage
# Need to distinguish persisting R cells (the rank would be one more)
# For R-R divsion, need to know whether it randomly choose the number from one sublineage
maxRRank <- with(dataMaxR, tapply(Rank, as.factor(Clone_uniqueID), max)) 

maxCellNumberLineagesR <- maxCellNumberLineages[match(names(maxRRank), names(maxCellNumberLineages))]

# hist(maxCellNumberLineages, xlab = "maximum number of cells in a lineage", ylab = "Number of lineages", 
#      main = "Distribution of maximun number of cells per lineage", cex.main = 0.9)

dataRankR <- data.frame(maxRRank, maxCellNumberLineagesR)

kt <- kruskal.test(maxCellNumberLineagesR ~ maxRRank)
boxplot(maxCellNumberLineagesR ~ maxRRank, ylab = "Max number of cells in lineage", xlab = "Number of successive R divisions in the lineage", pch = 16, cex = 0.6, 
        main = paste("Max number of cells vs \nnumber of R sucessive division in the lineage \np = ", round(kt$p.value,4)), cex.main = 0.9)
beeswarm(maxCellNumberLineagesR ~ maxRRank, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")
tab <- table(maxRRank)
for (i in seq_along(tab)){
  text(x = i, y = 2, labels = tab[i], cex = 0.8)  
}

#############################################
# Without the quiescent cells:
#############################################

dataTpNoQuiescent <- dataTp[!dataTp$Clone_uniqueID %in% quiescentClonesID, ] 
maxCellNumberLineagesNoQuiescent <- tapply(X = dataTpNoQuiescent$TotalCellNumber, as.factor(dataTpNoQuiescent$Clone_uniqueID),  max)

# pdf(file = "S3_MaxNumberOfCellPerLineage.pdf", width = 3, height = 4)
barPlotBeeswarm(maxCellNumberLineagesNoQuiescent, ylab = "Max number of cells per clone")
# dev.off()

barPlotBeeswarm(maxRRank, ylab = "Max number of cells per clone")

#############################################
# Comparison final and maximal number of cells
#############################################

# Final clone composition of non quiescent clones
cloneFinalCompositionNonQuiescent <- cloneFinalComposition[!cloneFinalComposition$Clone_uniqueID %in% quiescentClonesID , ]
finalCellNumberNonQuiescent <- cloneFinalCompositionNonQuiescent$TotalCellNumber

maxCellNonQuiescent <- maxCellNumberLineages[!names(maxCellNumberLineages) %in% quiescentClonesID ]
finalCellNumberNonQuiescent <- finalCellNumberNonQuiescent[match(cloneFinalCompositionNonQuiescent$Clone_uniqueID, names(maxCellNonQuiescent))]

cloneSize <- rbind(data.frame(cellNumber = maxCellNonQuiescent, type = rep("max", length(maxCellNonQuiescent))),
                   data.frame(cellNumber = finalCellNumberNonQuiescent , type = rep("final", length(finalCellNumberNonQuiescent))))

cloneSize <- data.frame(maxCellNonQuiescent, finalCellNumberNonQuiescent)

plot(0, 0, pch = 1, xlim = c(0.5,2.5), ylim = c(0,25), ylab = "# cells", 
     xlab = "Max number = 1, final number = 2", las = 1, main = "Relationship final and maximum \nnumber of cells in clones")

# Threshold for how many cells lost:  30% of total number of cells = average cell death Proportion in the dataset
deathRate <- 0.5
count <- 0
mycolor <- c("black")
cellLoss <- c(1:length(maxCellNonQuiescent))
for (i in seq_along(maxCellNonQuiescent)){
  cellLoss[i] <- (maxCellNonQuiescent[i] - finalCellNumberNonQuiescent[i]) / maxCellNonQuiescent[i]
  if(cellLoss[i] > deathRate){ 
    mycolor <- "red"
    count <- count + 1
  } else mycolor <- "black"  
  points(x= c(1,2), y = c(maxCellNonQuiescent[i],finalCellNumberNonQuiescent[i]), type = "o", pch = 18, col = mycolor)
}
legend("topright",legend = c("less than 50% loss", "more than 50% loss"), bty = "n", col = c("black","red"), lty = 1, lwd =2, cex = 0.8)


## Correlation cell loss in lineage and max number of cells
mycor <- cor.test(cellLoss, maxCellNonQuiescent)
myLm <- lm(maxCellNonQuiescent ~ cellLoss)
plot(cellLoss, jitter(maxCellNonQuiescent), xlab = "Loss of cells [% max]", ylab = "Max number of cells in the lineage", 
     main = paste("Relation cell loss and max size of lineage \n Pearson R = ", round(mycor$estimate,3), ", p = ", round(mycor$p.value,3)))
abline(myLm, col = "red", lty = 2)

## Correlation cell loss in lineage and final number of cells
mycor <- cor.test(cellLoss, finalCellNumberNonQuiescent, method = "pearson", exact = F)
cellLossPercent <- cellLoss*100
myLm <- lm(finalCellNumberNonQuiescent ~ cellLossPercent)
plot(cellLossPercent, jitter(finalCellNumberNonQuiescent), xlab = "Loss of cells [% max]", ylab = "Final number of cells in the lineage", 
     main = paste("Relation cell loss and final size of lineage \n Pearson correlation R = ", round(mycor$estimate,3), ", p = ", round(mycor$p.value,12)))
abline(myLm, col = "red", lty = 2)

################################################################################################################################
# Final clone composition
################################################################################################################################

myColors <- c("red","orange","cyan", "white", "grey")
myLegend <-  c("R","NR","N","Cell_Death","Uncertain")
rownames(cloneFinalComposition) <- cloneFinalComposition$Clone_uniqueID
barchart(as.matrix(cloneFinalComposition[cloneFinalComposition$numberSuccessiveDivision <= 3, 6:10]) ,  col = myColors, 
         xlab = "number of cells at the last recording, 3 or less successive divisions", key= list(text = list(myLegend), 
         space= "top", columns = 3, rectangles=list(col=myColors)))
barchart(as.matrix(cloneFinalComposition[cloneFinalComposition$numberSuccessiveDivision > 3, 6:10]) ,  col = myColors, 
         xlab = "number of cells at the last recording 4 or more divisions", key= list(text = list(myLegend), 
         space= "top", columns = 3, rectangles=list(col=myColors)))

################################################################################################################################
# Evolution of the number of cells over time
# Clone displayed according to their cluster determined by the hierarchical clustering earlier
################################################################################################################################

#############################################
## For simplified cell types
#############################################
## Not used anymore
## Sort the lineages by the max number of cells at any time point, and then define like 3 categories so that the plots are at the right scale
# maxCellNumberLineages <- tapply(X = dataTp$TotalCellNumber, as.factor(dataTp$Clone_uniqueID),  max)
# largeLineages <- levels(as.factor(dataTp$Clone_uniqueID))[which(maxCellNumberLineages>=10)]
# dataTpLarge <-  dataTp[dataTp$Clone_uniqueID %in% largeLineages,]
# 
# mediumLineages <- levels(as.factor(dataTp$Clone_uniqueID))[which(maxCellNumberLineages>=5 & maxCellNumberLineages<10 )]
# dataTpMedium <-  dataTp[dataTp$Clone_uniqueID %in% mediumLineages,]
# 
# smallLineages <- levels(as.factor(dataTp$Clone_uniqueID))[which(maxCellNumberLineages>1 & maxCellNumberLineages<5 )]
# dataTpSmall <-  dataTp[dataTp$Clone_uniqueID %in% smallLineages,]
# 
# quiencentLineages <- levels(as.factor(dataTp$Clone_uniqueID))[which(maxCellNumberLineages==1)]
# dataTpQuiescent <-  dataTp[dataTp$Clone_uniqueID %in% quiencentLineages,]


#############################################
## Alternatively use the clusters found by hierarchical clustering
#############################################

cluster1TreeID <- as.character(dataLineages$treeID[dataLineages$clusters==1])
cluster2TreeID <- as.character(dataLineages$treeID[dataLineages$clusters==2])
cluster3TreeID <- as.character(dataLineages$treeID[dataLineages$clusters==3])
cluster4TreeID <- as.character(dataLineages$treeID[dataLineages$clusters==4])

dataCluster1 <- dataTp[dataTp$Clone_uniqueID %in% cluster1TreeID,]
dataCluster2 <- dataTp[dataTp$Clone_uniqueID %in% cluster2TreeID,]
dataCluster3 <- dataTp[dataTp$Clone_uniqueID %in% cluster3TreeID,]
dataCluster4 <- dataTp[dataTp$Clone_uniqueID %in% cluster4TreeID,]

myColors <- c("red","orange","cyan","green","black", "purple")

my.settings <- list(superpose.symbol=list(col=myColors),
                    superpose.line=list(col=myColors, lwd = 2),
                    strip.background=list(col="black"))

xyplot(R + NR ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster1, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster1", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2,cex= 0.8),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(R + NR ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster2, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster2", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(R + NR ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster3, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster3", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(R + NR + N + Cell_Death  ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster1, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster1", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2, cex= 0.8),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       }) #+ TotalCellNumber


xyplot(R + NR + N + Cell_Death  ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster2, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster2", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       }) #+ TotalCellNumber

xyplot(R + NR + N + Cell_Death  ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster3, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster3", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       }) #+ TotalCellNumber

xyplot(R + NR + N + Cell_Death  ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster4, 
       ylab = "# of cells", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster4", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2, cex = 0.8),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

#############################################
# Total cell number
#############################################

myColours2 <- rainbow(15)
my.settings <- list(superpose.symbol=list(col=myColours2),
                    superpose.line=list(col=myColours2, lwd = 2)
)
xyplot(TotalCellNumber ~ Imaging_time_DPI, groups =  Clone_uniqueID, data = dataCluster1, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Total cell number cluster1",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(TotalCellNumber ~ Imaging_time_DPI ,groups =  Clone_uniqueID, data = dataCluster2, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Total cell number cluster2",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(TotalCellNumber ~ Imaging_time_DPI ,groups =  Clone_uniqueID, data = dataCluster3, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Total cell number cluster3",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(TotalCellNumber ~ Imaging_time_DPI ,groups =  Clone_uniqueID, data = dataCluster4, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Total cell number cluster4",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })



