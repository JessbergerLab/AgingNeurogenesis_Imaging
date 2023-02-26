################################################################################################################################
# Cell death
#############################################
# % of cell death is the number of cells undergoing cell death on the total number of cells (including dividing progenitors)
# Cell Death Proportion is the number of cells undergoing cell death on the number of leave cells (only the last cells in the tree not the intermediates)
# When not stated the cell death refers to all cell death (early and late), 
# cell death is considered early if it happens within `r lateDeathThreshold` days afer the cell birth.
################################################################################################################################

#############################################
## When does cell death occurs, counted as the number of days after induction (DPI) 
#############################################

# Because need to take the rate of death of leave cells
cellDeathFrequency <- getDeathRate(dataMax) 

dataDeath <- dataTp[dataTp$Cell_Death != 0,]
tabDpi <- table(dataDeath$Cell_Death, dataDeath$Imaging_time_DPI)
cellDeathFactor <- matrix(rep(as.numeric(rownames(tabDpi)), dim(tabDpi)[2]), ncol = dim(tabDpi)[2], nrow = length(rownames(tabDpi)))  
tabDeath <- tabDpi * cellDeathFactor
cellDeath <- colSums(tabDeath)

par(mfrow = c(1,1))
myBy <- 2
myTo <- round(max(dataDeath$Imaging_time_DPI)/myBy)*myBy
myBreaks <- seq(from= 0, to = myTo, by = myBy)
bins <- .bincode(as.numeric(names(cellDeath)), breaks = myBreaks)

binedDeath <- data.frame(as.numeric(names(cellDeath)), bins, NbCells = as.vector(cellDeath))

barplot(tapply(X = binedDeath$NbCells,INDEX = bins, FUN = sum), names.arg = myBreaks[unique(binedDeath$bins)], 
        ylab= "Number of cell death", xlab = "DPI intervals", main = paste("Timing of cell death (all types) n= ", 
        sum(cellDeath), "\n % of cell death (all cells)",round(cellDeathFrequency,1) , sep = " "), cex.names = 0.7, cex.main = 0.9)

#############################################
## All Cell death: types of cells
#############################################

# Table with info about cells that die 
infoCellDeath <- dataMax[dataMax$cellDeath == 0 , ] 

dataMaxSureType <- dataMax[dataMax$UncertaintyType3 == 1,]
infoCellDeathSureTypeBeforeDeath <- infoCellDeath[infoCellDeath$UncertaintyType3 == 1, ]

# Proportion of death for every cell type
tabDeath <- table(infoCellDeathSureTypeBeforeDeath$CellType3)

#tabDeath[4] <- 0
tabAll <- table(dataMaxSureType$CellType3)

par(mfrow = c(1,2))
my.hist <- barplot(as.vector(tabDeath/sum(tabDeath))*100, names.arg = c("R", "NR", "N"), col = c("red", "orange","cyan"), 
                   main = paste("Type of dying cells \nonly sure type last time point before cell death \n n = ", sum(tabDeath)), 
                   ylab = "% of dying cells", cex.main =  0.8)

# When there is no R-death, try: 
# my.hist <- barplot(as.vector(tabDeath/sum(tabDeath))*100, names.arg = c("NR", "N"), col = c("orange","cyan"), 
#                    main = paste("Type of dying cells \nonly sure type last time point before cell death \n n = ", sum(tabDeath)), 
#                    ylab = "% of dying cells", cex.main =  0.8)

for(i in seq_along(tabDeath)){
  text(x = my.hist[i], y = 5, labels = tabDeath[i], cex = 0.7)
}

deathPerType <- tabDeath/tabAll
# Different number of column, need to conform into the same structure


barplot(deathPerType*100, names.arg = c("R", "NR", "N"), col = c("red", "orange","cyan"), main = "Proportion of cell death per type", ylab = "% of dying cells in each type", cex.main = 0.9)

names(deathPerType) <- c("R", "NR", "N")
print("Proportion of each cell type undergoing cell death:")
print(round(deathPerType*100,2))

#############################################
## Early Cell death: types of cells
#############################################

infoCellDeathEarly <- dataMax[dataMax$cellDeathEarly == 0 , ] # table with info about cells that die 

infoCellDeathEarlySureTypeBeforeDeath <- infoCellDeathEarly[infoCellDeathEarly$UncertaintyType3 == 1, ]

# Proportion of early death on all cells for every cell type
tabDeathEarly <- table(infoCellDeathEarlySureTypeBeforeDeath$CellType3)

par(mfrow = c(1,3))
my.hist <- barplot(as.vector(tabDeathEarly/sum(tabDeath))*100, names.arg = c("R", "NR", "N"), col = c("red", "orange","cyan"), 
                   main = paste("Type of early dying cells \nonly sure type last time point before cell death \n n = ", sum(tabDeathEarly)), 
                   ylab = "% of dying cells", cex.main =  0.8)
for(i in seq_along(tabDeath)){
  text(x = my.hist[i], y = 5, labels = tabDeathEarly[i], cex = 0.7)
}


tabAll <- table(dataMaxSureType$CellType3)

deathPerTypeEarly <- tabDeathEarly/tabAll
barplot(deathPerTypeEarly * 100, names.arg = c("R", "NR", "N"), col = c("red", "orange","cyan"), main = "Proportion of early cell death per type", ylab = "% of early dying cells in each type", cex.main = 0.9)

names(deathPerTypeEarly) <- c("R", "NR", "N")
print("Proportion of each cell type undergoing early cell death:")
print(round(deathPerTypeEarly*100,2))

# Proportion of early death 
propEarlyDeath <- tabDeathEarly/tabDeath
barplot(as.vector(propEarlyDeath * 100), names.arg = c("R", "NR", "N"), col = c("red", "orange","cyan"), main = "Proportion of early cell death over all cell death", ylab = "% of early dying cells over all cell death", cex.main = 0.9)

#############################################
# Survival time after last division
#############################################

# All cells
# cellDeathID <- dataMinMax$Cell_uniqueID[dataMinMax$cellDeath == 0]
# dataDeath <- dataMinMax[dataMinMax$Cell_uniqueID %in% cellDeathID,]
# 
# deathTimingAfterDiv <- c(1:length(cellDeathID))
# for (i in seq_along(cellDeathID)){
#   deathDPI <- dataDeath$DPI[dataDeath$Cell_uniqueID == cellDeathID[i]]
#   deathTimingAfterDiv[i] <- abs(deathDPI[2] - deathDPI[1])
# }
# 
# hist(deathTimingAfterDiv, breaks = seq(from =0, to= 60, by = 1), main = "Cell death timing after last division", xlab = "survival time [Days]", ylab = "# cells")

# For certain neurons only
cellDeathID <- dataMinMax$Cell_uniqueID[dataMinMax$cellDeath == 0 & dataMinMax$CellType3 == 4 & dataMinMax$UncertaintyType3 == 1]
cellDeathID <- unique(cellDeathID)
dataDeath <- dataMinMax[dataMinMax$Cell_uniqueID %in% cellDeathID,]

deathTimingAfterDiv <- c(1:length(cellDeathID))
for (i in seq_along(cellDeathID)){
  deathDPI <- dataDeath$DPI[dataDeath$Cell_uniqueID == cellDeathID[i]]
  deathTimingAfterDiv[i] <- abs(deathDPI[2] - deathDPI[1])
}

par(mfrow = c(1,1))
hist(deathTimingAfterDiv, breaks = seq(from =0, to= 60, by = 1), main = paste("Cell death timing after last division, only certain Neurons,\n n = ",
                          length(deathTimingAfterDiv)), xlab = "survival time [Days]", ylab = "# Neurons",freq = T)
abline(v = lateDeathThreshold, lty = 2, col = "red")
# Adjust x- $ y-axis scale
# hist(deathTimingAfterDiv, breaks = seq(from =0, to= 60, by = 1), main = paste("Cell death timing after last division, only certain Neurons,\n n = ",
#                           length(deathTimingAfterDiv)), xlab = "survival time [Days]", ylab = "# Neurons",freq = T, xlim = c(0, 40), ylim = c(0, 80))

# pdf(file = "4B_beeswarm.pdf", width = 4.5, height = 4)
# beeswarm(as.matrix(deathTimingAfterDiv), method = "center", ylab = "survival time [Days]", xlab = "# Neurons", cex = 0.55, corral = "wrap")
# abline(h = lateDeathThreshold, lty = 2, col = "red")
# dev.off()

NeuronDeathTable <- dataMax[dataMax$Cell_uniqueID %in% cellDeathID,]
NeuronDeathTable$deathTimingAfterDiv <- deathTimingAfterDiv

DeathlineageIDWave1 <- NeuronDeathTable[NeuronDeathTable$deathTimingAfterDiv < lateDeathThreshold,]

DeathlineageIDWave2 <- NeuronDeathTable[NeuronDeathTable$deathTimingAfterDiv >= lateDeathThreshold,]


LineageIDFirstNotSecondWave <- unique(DeathlineageIDWave1$Clone_uniqueID)[! unique(DeathlineageIDWave1$Clone_uniqueID) %in% unique(DeathlineageIDWave2$Clone_uniqueID)]

# Table with only infos of the lineages with cell death only in the first and not second cell death wave
decodedDataLineageFirstNotSecondWave <- decodedData[decodedData$Clone_uniqueID %in% LineageIDFirstNotSecondWave,] #unique(DeathlineageIDWave1$Clone_uniqueID),]
# plotLineage2(decodedDataLineageFirstNotSecondWave,myTime = "DPI", saveGraph = TRUE, myFileName = "LineagesCellDeathFirstWaveNotSecondDeath.pdf")

# Plot of the neuron survival time within each lineage
plot(NeuronDeathTable$deathTimingAfterDiv ~ as.factor(NeuronDeathTable$Clone_uniqueID), col = "white", las =3, cex.axis = 0.5, xlab = "", ylab = "Cell death timing after last division")
beeswarm(deathTimingAfterDiv ~ as.factor(Clone_uniqueID), data = NeuronDeathTable, 
         add = T, pch = 20, col = add.alpha('red',alpha = 0.6), cex.axis = 0.7, cex = 0.8, corral = "wrap")
abline(h = lateDeathThreshold, lty = 2)

###################3
## Fitting of Cell death timing after last division: bimodal distribution?
####  Quite a bit of work


# Fit an exponential distrib compute the Kolmogorov smirnov statistic compared to data
# Sample from this distribution (like a lot) and do ks.test compared to fitted distribution and see whether the simulated distributions are different from original observed one
# 
# deathTimingAfterDiv <- deathTimingAfterDiv[deathTimingAfterDiv > 0]
# 
# library(fitdistrplus)
# 
# FitExp <- fitdist(deathTimingAfterDiv,distr = "exp")
# plot(FitExp)
# summary(FitExp)
# 
# ks.test(deathTimingAfterDiv, "pexp", 0.1414862, alternative = c("two.sided"))
# mu <- FitExp$estimate
# x <- seq(0, to = 10, by = 0.01)
# y <- mu * exp(-mu*x)
# plot(x,log(y))
# 
# FitLogNormal <- fitdist(log(deathTimingAfterDiv), distr = "norm")
# plot(FitLogNormal)
# 
# install.packages("mixtools")
# library(mixtools)
# mixture <- normalmixEM(log(deathTimingAfterDiv), k=1, fast=TRUE)
# plot(mixture)
# plot.mixEM(mixture



################################################################################################################################
# Occurence of cell death over time
################################################################################################################################

myColors <- "black"

my.settings <- list(superpose.symbol=list(col=myColors),
                    superpose.line=list(col=myColors, lwd = 2),
                    strip.background=list(col="black"),
                    strip.border=list(col="black"))

xyplot(cumCellDeath  ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster1, 
       ylab = "# of cell death", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes cluster1", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2, cex= 0.8), main = "Cumulative Cell Death cluster1",
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(cumCellDeath ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster2, 
       ylab = "# of cell death", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster2", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2), main = "Cumulative Cell Death cluster2",
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(cumCellDeath ~ Imaging_time_DPI | Clone_uniqueID, data = dataCluster3, 
       ylab = "# of cell death", type='b', cex= 0.5, 
       auto.key= list(space= "top", columns = 2, rectangle = FALSE, title = "CellTypes Cluster3", lines = TRUE, points = FALSE),
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2), main = "Cumulative Cell Death cluster3",
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

# Change the color
myColors <- rainbow(15)

my.settings <- list(superpose.symbol=list(col=myColors),
                    superpose.line=list(col=myColors, lwd = 2),
                    strip.background=list(col="black"),
                    strip.border=list(col="black"))

xyplot(cumCellDeath ~ Imaging_time_DPI , groups =  Clone_uniqueID, data = dataCluster1, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Cumulative Cell Death Cluster1",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(cumCellDeath ~ Imaging_time_DPI , groups =  Clone_uniqueID, data = dataCluster2, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Cumulative Cell Death cluster2",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })

xyplot(cumCellDeath ~ Imaging_time_DPI , groups =  Clone_uniqueID, data = dataCluster3, 
       ylab = "# of cells", type='b', cex= 0.5, main = "Cumulative Cell Death cluster3",
       auto.key= FALSE,
       par.settings = my.settings, 
       par.strip.text=list(col="white", font=2),
       panel=function(x,y,...){
         panel.grid(h=-1, v=0); 
         panel.xyplot(x,y,...)
       })


################################################################################################################################
# Neuronal Cell death and rank of division
#############################################
# I thought it might be also relevant to use the NR rank (For the NR cells rank 1 corresponds to the first apparition of an RN cell generated by a R cell in the branch). 
# NR rank 2 means that the neuron has been generated after one division of NR cells, rank 3 : 2 successive divisions of NR cells before neuron generation ... 
# Only neurons generated after 6 or 7 NR divisions tend to have a higher death rate.
# One has to keep in mind that the rank can be sligthly inaccurate due to lineage uncertainty (some divisions might have been missed).
################################################################################################################################

#############################################
## Frequency of Neuron cell death according to the rank of appearance in the tree
#############################################

neuronDataMax <- dataMax[dataMax$CellType3 == 4 & dataMax$UncertaintyType3 == 1, ]
tab <- with(neuronDataMax, table(cellDeath, Rank))
propTab  <- prop.table(tab,2)
marginTab <- margin.table(tab,2)
Ftest <- fisher.test(tab[,2:3])
myPlot <- barplot(propTab[1,]*100, xlab = "Rank of neuron production", ylab= "Neuron death Proportion [%]", 
                  main = paste("Neuron cell death according to their rank of appearance \n red line = overall neuron death Proportion, numbers in bars = # of neurons produced at each rank \n rank 3 vs 4 p = ", 
                  round(Ftest$p.value,4)), cex.main = 0.85)
abline(h = deathPerType[3]*100, col = "red", lwd = 2)
text(x = myPlot, y = 3, labels = marginTab, cex = 0.7)

#############################################
## Unsing the NR rank, including only the NR cells
#############################################

neuronMotherID <- neuronDataMax$Mother_uniqueID[neuronDataMax$UncertaintyLineage == 1]
neuronDataMaxNRMother <- neuronDataMax[dataMax$CellType2[match(neuronMotherID,dataMax$Cell_uniqueID)] == 2,]
tab <- with(neuronDataMaxNRMother, table(cellDeath, RankNR))
propTab  <- prop.table(tab,2)
marginTab <- margin.table(tab,2)
#print("Statistical test for influence of NR rank on cell death Proportion")
Ftest <- chisq.test(tab[,1:2])
myPlot <- barplot(propTab[1,]*100, xlab = "Rank of neuron production \n(NR rank, only NR mother cells certain transitions)", ylab= "Neuron death Proportion [%]", main = paste("Neuron cell death according to their rank of appearance (NR rank) \n Rank 2 vs 3  p = ", round(Ftest$p.value,4)), cex.main = 0.8, cex.lab = 0.8)
abline(h = deathPerType[3]*100, col = "red", lwd = 2)
text(x = myPlot, y = 3, labels = marginTab, cex = 0.7)

# print("Statistical test for difference rank 2 vs 3")
# chisq.test(tab[,1:2])  # Prints test summary
# print("Statistical test for difference rank 3 vs 4")
# chisq.test(tab[,2:3])  # Prints test summary

#############################################
## With only early cell death:
#############################################

tab <- with(neuronDataMax, table(cellDeathEarly, Rank))
propTab  <- prop.table(tab,2)
marginTab <- margin.table(tab,2)
Ftest <- fisher.test(tab[,2:3])
myPlot <- barplot(propTab[1,]*100, xlab = "Rank of neuron production", ylab= "Neuron early death Proportion [%]", 
                  main = paste("Neuron early cell death according to their rank of appearance \n rank 3 vs 4 p = ", round(Ftest$p.value,4)), cex.main = 0.8)
abline(h = deathPerTypeEarly[3]*100, col = "red", lwd = 2)
text(x = myPlot, y = 3, labels = marginTab, cex = 0.7)

#############################################
# NR rank NR cells
#############################################

tab <- with(neuronDataMaxNRMother, table(cellDeathEarly, RankNR))
propTab  <- prop.table(tab,2)
marginTab <- margin.table(tab,2)
#print("Statistical test for influence of NR rank on cell death Proportion")
Ftest <- chisq.test(tab[,1:2])
myPlot <- barplot(propTab[1,]*100, xlab = "Rank of neuron production \n(NR rank, only NR mother cells certain transitions)", ylab = "Neuron early death Proportion [%]", 
                  main = paste("Neuron cell death according to their rank of appearance (NR rank) \n red line = overall neuron early death Proportion \n Rank 2 vs 3  p = ", 
                  round(Ftest$p.value,4)), cex.main = 0.85, cex.lab = 0.8)
abline(h = deathPerTypeEarly[3]*100, col = "red", lwd = 2)
text(x = myPlot, y = 3, labels = marginTab, cex = 0.7)


################################################################################################################################
# Symmetry of survival for sister neurons
################################################################################################################################

# Might need to filter out the trees with particular short recording time, 
# because we are not sure about the ultimate survival of neurons: this bias is removed when we look only at early cell death.

# Generates the division table using CellType3, simplified cell types but does not indicate dead cells
divisionTable <- transformTreeIntoDivision(decodedData, cellType = "simplified no death")   

# Filters for only sure types and transitions
filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyLineage = c("1")) 

# Filters for 2 daughters = neurons
filteredDivisionTableNeurons <- filteredDivisionTable[filteredDivisionTable$D1CellType == 4 & filteredDivisionTable$D2CellType == 4 & filteredDivisionTable$D1UncertaintyType == 1 & filteredDivisionTable$D2UncertaintyType == 1,] 

# Rank of the division initially, need to add one to the rank of the mother cell to get the rank of production of the neuron
filteredDivisionTableNeurons$NeuronRank <- as.factor(as.numeric(as.character(filteredDivisionTableNeurons$DivisionRank))+1) 

# Rank of the division initially, need to add one to the rank of the mother cell the get the rank of production of the neuron
filteredDivisionTableNeurons$NeuronRankNR <- as.factor(as.numeric(as.character(filteredDivisionTableNeurons$DivisionRankNR))+1) 

## Adds a column survival for the neurons
# Creates a new column that gives the joined fate of both daughter neurons
filteredDivisionTableNeurons$Survival <- filteredDivisionTableNeurons$D1Death  

bothDeathID <- NULL
bothSurvivalID <- NULL
MixedID <- NULL 

for (i in seq_along(filteredDivisionTableNeurons$D1ID)){
  if (filteredDivisionTableNeurons$D1Death[i] == 0 & filteredDivisionTableNeurons$D2Death[i] == 0){
    filteredDivisionTableNeurons$Survival[i] <- "2_Deaths"
    bothDeathID <- c(bothDeathID, filteredDivisionTableNeurons$D1ID[i], filteredDivisionTableNeurons$D2ID[i])
  } else if(filteredDivisionTableNeurons$D1Death[i] == 1 & filteredDivisionTableNeurons$D2Death[i] == 0 | filteredDivisionTableNeurons$D1Death[i] == 0 & filteredDivisionTableNeurons$D2Death[i] == 1){
    filteredDivisionTableNeurons$Survival[i] <- "1Death_1Survival"
    MixedID <- c(MixedID, filteredDivisionTableNeurons$D1ID[i], filteredDivisionTableNeurons$D2ID[i])
  } else {
    filteredDivisionTableNeurons$Survival[i] <- "2_Survivals"
    bothSurvivalID <- c(bothSurvivalID, filteredDivisionTableNeurons$D1ID[i], filteredDivisionTableNeurons$D2ID[i])
  }
}

## Adds a column early survival for the neurons, survives longer than the defined threshold for early cell death
# Creates a new column that gives the joined fate of both daughter neurons
filteredDivisionTableNeurons$SurvivalEarly <- filteredDivisionTableNeurons$D1Death  


bothDeathIDEarly <- NULL
bothSurvivalIDEarly <- NULL
MixedIDEarly <- NULL 

for (i in seq_along(filteredDivisionTableNeurons$D1ID)){
  if (filteredDivisionTableNeurons$D1DeathEarly[i] == 0 & filteredDivisionTableNeurons$D2DeathEarly[i] == 0){
    filteredDivisionTableNeurons$SurvivalEarly[i] <- "2_Deaths"
    bothDeathIDEarly <- c(bothDeathIDEarly, filteredDivisionTableNeurons$D1ID[i], filteredDivisionTableNeurons$D2ID[i])
  } else if(filteredDivisionTableNeurons$D1DeathEarly[i] == 1 & filteredDivisionTableNeurons$D2DeathEarly[i] == 0 | filteredDivisionTableNeurons$D1DeathEarly[i] == 0 & filteredDivisionTableNeurons$D2DeathEarly[i] == 1){
    filteredDivisionTableNeurons$SurvivalEarly[i] <- "1Death_1Survival"
    MixedIDEarly <- c(MixedIDEarly, filteredDivisionTableNeurons$D1ID[i], filteredDivisionTableNeurons$D2ID[i])
  } else {
    filteredDivisionTableNeurons$SurvivalEarly[i] <- "2_Survivals"
    bothSurvivalIDEarly <- c(bothSurvivalIDEarly, filteredDivisionTableNeurons$D1ID[i], filteredDivisionTableNeurons$D2ID[i])
  }
}


#############################################
# Proportion of survival of neurons pairs according the mother type
#############################################

tab <- table(filteredDivisionTableNeurons$Survival, filteredDivisionTableNeurons$MotherCellType)
barplot(t(prop.table(tab)*100), cex.names = 0.7, main = paste("Survival pattern of sister neurons,\n n = ", sum(tab)," pairs"), ylab = "% of neuron pairs", col = c("red", "orange"), las = 0)
legend("topright",legend = c("Mother R", "Mother NR"), fill = c("red", "orange"), bty = "n")

tab <- table(filteredDivisionTableNeurons$SurvivalEarly, filteredDivisionTableNeurons$MotherCellType)
barplot(t(prop.table(tab)*100), cex.names = 0.7, main = paste("Early Survival pattern of sister neurons,\n n = ", sum(tab)," pairs"), ylab = "% of neuron pairs", col = c("red", "orange"), las = 0)
legend("topright",legend = c("Mother R", "Mother NR"), fill = c("red", "orange"), bty = "n")

#############################################
### @@@!!@!@!@ Here check if the assomption of independence between early cell death and survival is satisfied. 
### It seems in the first place that it is not
#############################################

#############################################
# Proportion of survival of neurons pairs according the rank of neuron generation
#############################################

tab <- table(filteredDivisionTableNeurons$Survival, filteredDivisionTableNeurons$NeuronRank)
myPlot <- barplot(prop.table(tab,2)*100, cex.names = 0.8, main = paste("Survival pattern of sister neurons,\n n = ", sum(tab)," pairs, numbers on the graph show \nthe numbers of pairs at each rank"), ylab = "% of neuron pairs at each rank", col = c("gray", "black", "white"), beside = TRUE, cex.main = 0.8, xlab = "Rank of neuron generation")
legend("top", legend = rownames(tab), fill =  c("gray", "black", "white"), bty = "n", cex = 0.8)
text(myPlot[2,]+ 1, 3, margin.table(tab,2), cex = 0.8)

#############################################
# Proportion of survival of neurons pairs according the NR rank of neuron generation
#############################################

# Mother cell NR to use the NR rank
filteredDivisionTableNeuronsNRMother <- filteredDivisionTableNeurons[filteredDivisionTableNeurons$MotherCellType ==2 ,]
tab <- table(filteredDivisionTableNeuronsNRMother$Survival, filteredDivisionTableNeuronsNRMother$NeuronRankNR)

myPlot <- barplot(prop.table(tab,2)*100, cex.names = 0.8, main = paste("Survival pattern of sister neurons,\n n = ", sum(tab)," pairs, numbers on the graph show \nthe numbers of pairs at each rank"), ylab = "% of neuron pairs at each rank", col = c("gray", "black", "white"), beside = TRUE, cex.main = 0.8, xlab = "NR rank neuron generation")
legend("top", legend = rownames(tab), fill =  c("gray", "black", "white"), bty = "n", cex = 0.8)
text(myPlot[2,] +1, 3, margin.table(tab,2), cex = 0.8)


#############################################
## Early cell death only
#############################################
#############################################
# Proportion of survival of neurons pairs according the rank of neuron generation
#############################################

tab <- table(filteredDivisionTableNeurons$SurvivalEarly, filteredDivisionTableNeurons$NeuronRank)
myPlot <- barplot(prop.table(tab,2)*100, cex.names = 0.8, main = paste("Early Survival pattern of sister neurons,\n n = ", sum(tab)," pairs, numbers on the graph show \nthe numbers of pairs at each rank"), ylab = "% of neuron pairs at each rank", col = c("gray", "black", "white"), beside = TRUE, cex.main = 0.8, xlab = "Rank of neuron generation")
legend("top", legend = rownames(tab), fill =  c("gray", "black", "white"), bty = "n", cex = 0.8)
text(myPlot[2,]+ 1, 3, margin.table(tab,2), cex = 0.8)

#############################################
# Proportion of survival of neurons pairs according the NR rank of neuron generation
#############################################

# Mother cell NR to use the NR rank
tab <- table(filteredDivisionTableNeuronsNRMother$SurvivalEarly, filteredDivisionTableNeuronsNRMother$NeuronRankNR)

myPlot <- barplot(prop.table(tab,2)*100, cex.names = 0.8, main = paste("Early Survival pattern of sister neurons,\n n = ", sum(tab)," pairs, numbers on the graph show \nthe numbers of pairs at each rank"), ylab = "% of neuron pairs at each rank", col = c("gray", "black", "white"), beside = TRUE, cex.main = 0.8, xlab = "NR rank neuron generation")
legend("top", legend = rownames(tab), fill =  c("gray", "black", "white"), bty = "n", cex = 0.8)
text(myPlot[2,] +1, 3, margin.table(tab,2), cex = 0.8)



################################################################################################################################
## Are the cell death probabilities of sister neurons independent from one another?
#############################################

# No evidence to exclude that the observed outcome could have been generated by a mechanism where the cell death of sister neurons 
# is independent from one another both for all cell death and early cell death (cf analysis with David's help). 
# There is no evidence to reject the hypothesis that the observed distribution can be generated from a binomial distribution 
# (assumption of independence, proba of early cell death p = 0.38, proba of overall cell death p = 0.49, maximum likelihood estimation followed by Chi square goodness of fit run in Mathematica).

```{r neuron_DyingSurvivingpair_Simulation, include = TRUE, echo=FALSE,  fig.width=4, fig.height=4, fig.show='hold'} 
# Simulation of death vs survival type for pair of neurons using the indepence hypothesis.
# simulationNeuronPair <- function(probDying = p1, iteration = 100, npair = 96, dataTable = tab, myMain = "Surviving Dying of sister neurons"){
#   
#   # early = 1 ; late = 0
#   # dying = 0; surviving = 1
#   
#   tabSim <- matrix(nrow= iteration, ncol = 3)
#   
#   for(iter in c(1:iteration)){
#     pair <- data.frame(D1 = vector("numeric", length = npair), D2 = vector("numeric", length = npair), combined = vector("numeric", length = npair))
#     for (i in c(1:npair)){
#       randomDeath <- runif(2)
#       if(randomDeath[1]< probDying){pair[i,1] <- 0} else {pair[i,1] <- 1}
#       if(randomDeath[2]< probDying){pair[i,2] <- 0} else {pair[i,2] <- 1}
#     }
#     
#     pair[,3] <- rowSums(pair,na.rm = T) 
#     pair$combined <- factor(x = pair$combined, levels = c(0,1,2))
#     
#     tabSim[iter,] <- table(pair[,3])
#   }
#   
#   meanSim <- apply(tabSim, MARGIN = 2, FUN = mean)
#   sdSim <- apply(tabSim/rowSums(tabSim,1), MARGIN = 2, FUN = sd)
#   globalDiff <- fisher.test(matrix(c(round(meanSim,0), dataTable) ,nrow = 3, ncol = 2))
#   
#   # In the simulation how many iteration are significantly different from observed data
#   fisherTestFun <- function(x, tab = dataTable){
#     test <- fisher.test(matrix(c(x,tab),nrow = 3, ncol = 2))
#     test$p.value
#   }
#   pValueSimVsObserved <- apply(tabSim, MARGIN = 1, FUN = fisherTestFun, tab = dataTable)
#   differenceWithObserved <- sum(p.adjust(pValueSimVsObserved, method = "bonferroni") < 0.5)
#   
# 
#   propTab <- prop.table(dataTable)
#   propTabRecap <- rbind(propTab*100, prop.table(meanSim)*100)
#   myPlot <- barplot(propTabRecap, beside = T,cex.names = 0.9, main = paste(myMain,"\nproba cell death = ", round(probDying,2)), ylab = "% of neuron pairs", col = c("white","gray"), las = 0, cex.main = 0.8)
#   legend("topleft", legend = c("Observed", "Simulation"), fill = c("white","gray"), bty = "n")
#   error.bar(x = myPlot[2,], y = propTabRecap[2,], upper = sdSim)
#   
#   cost <- sum((propTab - meanSim/100)^2)
#   list(cost,globalDiff$p.value, differenceWithObserved)
# }
```


```{r independenceAssumption, echo=FALSE, fig.width=4, fig.height=4}
# #Assuming the 2 neurons of the pair behave totally independently. Here is what the outcome would be if we would simulate neuron pairs.
# tab <- table(filteredDivisionTableNeurons$Survival)
# propTab <- prop.table(tab)
# p1 <- sqrt(propTab[2])
# p2 <- 1-sqrt(propTab[3])
# delta <- 4 - 4 * 2 * propTab[1]
# p3.1 <- (2 - sqrt(delta))/4
# p3.2 <- (2 + sqrt(delta))/4
# 
# # tab for simulation must be 1st both dying, 2nd 1death 1 surviving, 3rd both surviving
# tab <- tab[c(2,1,3)]
# simulationNeuronPair(probDying = p1,iteration = 100, npair = 96, dataTable = tab, myMain = "Sister neurons survival vs death \n simulation independence")
# simulationNeuronPair(probDying = p2,iteration = 100, npair = 96, dataTable = tab, myMain = "Sister neurons survival vs death \nsimulation independence")
# simulationNeuronPair(probDying = p3.2,iteration = 100, npair = 96, dataTable = tab, myMain = "Sister neurons survival vs death \nsimulation independence")

# 
# This hypothesis corresponds to assuming no correlation between the neurons in the pair.
# Let pose $p$ the probability of a neuron dying, the probability of a neuron surviving would then be $1-p$
# 
#   * The probability of having 2 sister neurons dying is  $p^2$
#   
#   * The probability of having 2 sister neurons surviving is  $(1-p)^2$
# 
#   * The probability of having 1 neuron dying and 1 surviving is $2p(1-p)$
# 
# From the observations :
# $$p^2 = 0.33  (1)$$
# $$p = 0.58$$ 
# 
# $$(1-p)^2 = 0.21  (2)$$
# $$p = 0.54$$ 
# 
# $$2p(1-p) = 0.46  (3)$$
# $$p = 0.64$$  or  $$p = 0.35$$
# 
# We can find a $p$ that satisfies the 3 equations, so we do not have evidence against the independence of the death or survival of sister neurons. There is no correlation between the fate of sister neurons.
```


```{r independenceAssumptionEarlyDeath, echo=FALSE, fig.width=4, fig.height=4}
## Inpedendence of sister neurons, only early cell death
# tab <- table(filteredDivisionTableNeurons$SurvivalEarly)
# propTab <- prop.table(tab)
# p1 <- sqrt(propTab[2])
# p2 <- 1-sqrt(propTab[3])
# delta <- 4 - 4 * 2 * propTab[1]
# p3.1 <- (2 - sqrt(delta))/4
# p3.2 <- (2 + sqrt(delta))/4
# 
# This hypothesis corresponds to assuming no correlation between the neurons in the pair.
# Let pose $p$ the probability of a neuron dying, the probability of a neuron surviving would then be $1-p$
# 
#   * The probability of having 2 sister neurons dying is  $p^2$
#   
#   * The probability of having 2 sister neurons surviving is  $(1-p)^2$
# 
#   * The probability of having 1 neuron dying and 1 surviving is $2p(1-p)$
# 
# From the observations :
# $$p^2 = 0.18  (1)$$
# $$p = 0.42$$ 
# 
# $$(1-p)^2 = 0.42  (2)$$
# $$p = 0.35$$ 
# 
# $$2p(1-p) = 0.41  (3)$$
# $$p = 0.28$$  or  $$p = 0.78$$

```


```{r independenceAssumptionEarlyDeath_Simulation, echo=FALSE, fig.width=4, fig.height=4}
# tab for simulation must be 1st both dying, 2nd 1death 1 surviving, 3rd both surviving
# tab <- tab[c(2,1,3)]
# simulationNeuronPair(probDying = p1,iteration = 100, npair = 96, dataTable = tab, myMain = "Sister neurons early survival vs death \n simulation independence")
# simulationNeuronPair(probDying = p2,iteration = 100, npair = 96, dataTable = tab, myMain = "Sister neurons early survival vs death \n simulation independence")
# simulationNeuronPair(probDying = p3.1,iteration = 100, npair = 96, dataTable = tab, myMain = "Sister neurons early survival vs death \n simulation independence")
# #simulationNeuronPair(probDying = p3.2,iteration = 100, npair = 100, dataTable = tab, myMain = "Sister neurons early survival vs death \n simulation independence")



################################################################################################################################
# Early vs late neuronal death
################################################################################################################################

#############################################
# From the distribution of timing of neuronal death, one can detect 2 distributions corresponding to early cell death 
# (occurs within 'r lateDeathThreshold' days after birth), late cell death (occurs more than 10 days after birth).
#############################################

deathType <- c("early", "late")
NeuronDeathTable$DeathType <- deathType[(NeuronDeathTable$deathTimingAfterDiv >= lateDeathThreshold)*1 +1]

# Take only the number of early death: 0 means both daugther late death, 1 means 1 early 1 late, 2 means both early
tab <- table(NeuronDeathTable$DeathType) 

barplot(t(prop.table(tab)*100), cex.names = 0.9, main = paste("Death types of all neurons,\n n = ", sum(tab)," neurons"), 
        ylab = "% of dying neurons", col = c("white"), las = 0)

#############################################
# Cell death type and survival duration and time before death for mixed pairs
#############################################

deathIDofMixed <- match(MixedID, NeuronDeathTable$Cell_uniqueID)
deathIDofMixed <- deathIDofMixed[complete.cases(deathIDofMixed)]
dataMixedDeath  <- NeuronDeathTable[deathIDofMixed,]

dataMixedDeath$DeathType <- deathType[(dataMixedDeath$deathTimingAfterDiv >= lateDeathThreshold)*1 + 1]

# Take only the number of early death: 0 means both daugther late death, 1 means 1 early 1 late, 2 means both early
tabMixed <- table(dataMixedDeath$DeathType) 

tabComparison <- rbind(tab, tabMixed)
Xsq <- chisq.test(tabComparison)


barplot(t(prop.table(tabMixed)*100), cex.names = 0.9, 
        main = paste("Death type of sister neurons mixed survival (1 death 1 survival),\n n = ", 
        sum(tabMixed)," pairs \n different proportion than full dataset p = ", round(Xsq$p.value,4)), 
        ylab = "% of dying neurons", col = c("white"), las = 0, cex.main = 0.8)

#############################################
## Survival time
#############################################

dataMixedPair <- dataMinMax[dataMinMax$Cell_uniqueID %in% MixedID,] 

endTimingAfterDiv <- c(1:length(MixedID))
for (i in seq_along(MixedID)){
  endDPI <- dataMixedPair$DPI[dataMixedPair$Cell_uniqueID == MixedID[i]]
  endTimingAfterDiv[i] <- abs(endDPI[2] - endDPI[1])
}

MixedID <- t(matrix(MixedID, ncol = length(MixedID)/2))
endTimingAfterDiv <- t(matrix(endTimingAfterDiv, ncol = length(endTimingAfterDiv)/2))

# Need to make sure that all the cells that die (column 1 shorter survival duration) or survive (column 2) are in the same column. 
for (i in c(1:nrow(endTimingAfterDiv))){
  endTimingAfterDiv[i,1:2] <- endTimingAfterDiv[i,order(endTimingAfterDiv[i,1:2], decreasing = F)]
}

plot(endTimingAfterDiv[,1], endTimingAfterDiv[,2], xlab = "Survival time dying cell [Days]", 
     ylab = "Survival time surviving cell [Days]", pch = 16, cex = 0.7, xlim = c(0,35), ylim = c(0,35), 
     main = "Survival time of mixed neurons pairs", cex.main = 0.9)

#############################################
# Cell death type and survival duration of pairs with 2 deaths
#############################################

dataBothDeath  <- NeuronDeathTable[match(bothDeathID, NeuronDeathTable$Cell_uniqueID),]

dataBothDeath$DeathType <- deathType[(dataBothDeath$deathTimingAfterDiv >= lateDeathThreshold)*1 +1]

tab <- table(table(dataBothDeath$DeathType, dataBothDeath$Mother_uniqueID)[1,]) # take only the number of early death: 0 means both daugther late death, 1 means 1 early 1 late, 2 means both early
names(tab) <- c("both late", "1 early 1 late", "both early")
barplot(t(prop.table(tab)*100), cex.names = 0.9, main = paste("Death type of sister neurons,\n n = ", sum(tab)," pairs"), 
        ylab = "% of neuron pairs", col = c("white"), las = 0)

bothDeathID <- t(matrix(dataBothDeath$Cell_uniqueID, ncol = length(dataBothDeath$Cell_uniqueID)/2))
deathTimingAfterDiv <- t(matrix(dataBothDeath$deathTimingAfterDiv, ncol = length(dataBothDeath$deathTimingAfterDiv)/2))

# For a better visualization, cell 1 has the shortest survival time and Cell 2 the longest
for (i in c(1:nrow(deathTimingAfterDiv))){
  deathTimingAfterDiv[i,1:2] <- deathTimingAfterDiv[i,order(deathTimingAfterDiv[i,1:2], decreasing = F)]
}

plot(deathTimingAfterDiv[,1], deathTimingAfterDiv[,2], xlab = "Shorter Survival time cell 1 [Days]", 
     ylab = "Longer Survival time cell 2 [Days]", pch = 16, cex = 0.7, xlim = c(0,35), ylim = c(0,35), 
     main = "Survival time of both dying neurons pairs", cex.main = 0.9)

#############################################
# Assuming the 2 neurons of the both dying pair behave totally independently
#############################################

```{r independenceAssumptionDeathType, echo=FALSE, fig.width=8, fig.height=4}
# tab <-  table(table(dataBothDeath$DeathType, dataBothDeath$Mother_uniqueID)[1,]) # take only the number of early death: 0 means both daugther late death,
# propTab <- prop.table(tab)
# names(propTab) <- c("both late", "1 early 1 late", "both early")
# p1 <- sqrt(propTab[3])
# p2 <- 1-sqrt(propTab[1])
# delta <- 4 - 4 * 2 * propTab[2]
# p3.1 <- (2 - sqrt(delta))/4
# p3.2 <- (2 + sqrt(delta))/4


# This hypothesis corresponds to assuming no correlation between the neurons in the pair.
# Let pose $p$ the probability of a neuron dying early, the probability of a neuron undergoing a late cell death would then be $1-p$
# 
#   * The probability of having 2 sister neurons dying early  $p^2$
#   
#   * The probability of having 2 sister neurons dying late  $(1-p)^2$
# 
#   * The probability of having 1 neuron dying early and 1 dying late is $2p(1-p)$
# 
# From the observations :
# $$p^2 = 0.56  (1)$$
# $$p = 0.71$$ 
# 
# $$(1-p)^2 = 0.093  (2)$$
# $$p = 0.60$$ 
# 
# $$2p(1-p) = 0.34  (3)$$
# $$p = 0.78$$ or $$p = 0.22$$
# 
# We can find a $p$ that satisfies the 3 equations, so we do not have evidence against the independence of the death types of sister neurons. This result is supported by the simulation of pairs of neuron death type using the independence assumption and $p = `r p1`$.

```


```{r neuron_Dyingpair_Simulation, include = TRUE, echo=FALSE,  fig.width=4, fig.height=4, fig.show='hold'} 
# # Simulation of death type for pair of dying neurons using the indepence hypothesis.
# 
# # early = 1 ; late = 0
# iteration <- 100
# tabSim <- matrix(nrow=100, ncol = 3)
# npair <- 32
# 
# 
# for(iter in c(1:iteration)){
#   pair <- data.frame(D1 = vector("numeric", length = npair), D2 = vector("numeric", length = npair), combined = vector("numeric", length = npair))
#   for (i in c(1:npair)){
#     randomDeath <- runif(2)
#     if(randomDeath[1]< p1){pair[i,1] <- 1} else {pair[i,1] <- 0}
#     if(randomDeath[2]< p1){pair[i,2] <- 1} else {pair[i,2] <- 0}
#   }
#   
#   pair[,3] <- rowSums(pair,na.rm = T) 
#   pair$combined <- factor(x = pair$combined, levels = c(0,1,2))
#   
#   tabSim[iter,] <- table(pair[,3])
# }
# 
# meanSim <- apply(tabSim, MARGIN = 2, FUN = mean)
# sdSim <- apply(tabSim, MARGIN = 2, FUN = sd)
# 
# propTabRecap <- rbind(propTab*100, prop.table(meanSim)*100)
# myPlot <- barplot(propTabRecap, beside = T,cex.names = 0.9, main = paste("Death type of sister neurons,\n n = ", sum(tab)," pairs"), ylab = "% of neuron pairs", col = c("white","gray"), las = 0)
# legend("topleft", legend = c("Observed", "Simulation"), fill = c("white","gray"), bty = "n")
# error.bar(x = myPlot[2,], y = propTabRecap[2,], upper = sdSim)


################################################################################################################################
# Neural pairs where both neurons survive
################################################################################################################################

#############################################
# Survival duration of pairs with 2 survivors
#############################################

dataBothSurvival <- dataMinMax[dataMinMax$Cell_uniqueID %in% bothSurvivalID,]

survivalTimingAfterDiv <- c(1:length(bothSurvivalID))
for (i in seq_along(bothSurvivalID)){
  survivalDPI <- dataBothSurvival$DPI[dataBothSurvival$Cell_uniqueID %in% bothSurvivalID[i]]
  survivalTimingAfterDiv[i] <- abs(survivalDPI[2] - survivalDPI[1])
}

bothSurvivalID <- t(matrix(bothSurvivalID, ncol = length(bothSurvivalID)/2))
survivalTimingAfterDiv <- t(matrix(survivalTimingAfterDiv, ncol = length(survivalTimingAfterDiv)/2))

cloneFinalComposition <- cloneOutcome(dataTp) # composition of the clone at the last observed time point 


bothSurvivalTreeID <- dataMax$Clone_uniqueID[match(bothSurvivalID[,1],dataMax$Cell_uniqueID)]
bothSurvivalTreeMaxDPI <- vector(length = length(bothSurvivalTreeID))
for (i in seq_along(bothSurvivalID[,1])){
  bothSurvivalTreeMaxDPI[i] <- cloneFinalComposition$Imaging_time_DPI[cloneFinalComposition$Clone_uniqueID == bothSurvivalTreeID[i]]
}

plot(survivalTimingAfterDiv[,1], bothSurvivalTreeMaxDPI, xlab = "Survival time pairs [Days]", ylab = "Imaging time [Days]", pch = 16, cex = 0.7, main = "Survival time of both surviving neurons pairs \nvs imaging time", cex.main = 0.9, xlim = c(0,60))
beeswarm(survivalTimingAfterDiv[,1], main = "Distribution of the survival time of survival pairs", xlab = "# of pairs", ylab = "Survival time [Days]", cex.main = 0.9)
bxplot(survivalTimingAfterDiv[,1], col = "red", add = T)

################################################################################################################################
# Cell death Proportion per lineage and sublineages (all cell death included early + late)
################################################################################################################################

# There is a high heterogeneity of cell death Proportions 
# (% of cell death in the leaf cells) in the different lineages, with a distribution skewed towards the high death Proportions.
# For the analysis in the different subtrees. Every tree is first divided in secondary subtrees (trees rank 2, might be more than 2 in case of uncertain transitions), 
# every secondary tree is divided into tertiary subtrees (trees rank 3). The death Proportion of the different subtree is then computed.

nonQuiescentID <- as.character(dataLineages$treeID[!dataLineages$treeID %in% quiescentClonesID ])
dataLineagesNoQuiescent <- dataLineages[!dataLineages$treeID %in% quiescentClonesID,]

par(mfrow =c(1,2))
beeswarm(dataLineagesNoQuiescent$deathRate, xlab = "Number of lineages", ylab = "Cell death Proportion [%]", main = paste("Cell death Proportion per lineage, quiescent cells removed \n overall cell death Proportion of leave cells =", round(cellDeathFrequency,2),"%"), cex.main = 0.8)
bxplot(dataLineagesNoQuiescent$deathRate, col = "red", add = T)
barPlotBeeswarm(dataLineagesNoQuiescent$deathRate, xlab = "Number of lineages", ylab = "Cell death Proportion [%]", cex.main = 0.8)

meanDeathRate <- mean(dataLineagesNoQuiescent$deathRate)
semDeathRate <- sd(dataLineagesNoQuiescent$deathRate)/sqrt(length(dataLineagesNoQuiescent$deathRate))

par(mfrow =c(1,1))
#pdf("Fig4A_hist.pdf", width = 4, height = 4)
hist(dataLineagesNoQuiescent$deathRate, ylab = "Number of lineages", xlab = "Cell death Proportion [%]", main = paste("Cell death Proportion per lineage, quiescent cells removed \n mean =", round(meanDeathRate,2),"%, sem = ", round(semDeathRate,2)), breaks = seq(0,100,10))
abline(v = mean(dataLineagesNoQuiescent$deathRate), col = "red", lwd = 2)
arrows(x0 = meanDeathRate - semDeathRate,y0 = 6, x1 = meanDeathRate + semDeathRate,y1 = 6,code = 3, angle = 90, length = 0.05, lwd = 2)
#dev.off()

boxplot(dataLineagesNoQuiescent$deathRate ~ dataLineagesNoQuiescent$mouseID, xlab = "Mouse ID", ylab = "Cell death Proportion [%]")
beeswarm(dataLineagesNoQuiescent$deathRate ~ dataLineagesNoQuiescent$mouseID, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

# library(dunn.test)
# library(rcompanion)
# dt <- dunn.test(dataLineagesNoQuiescent$deathRate, as.factor(dataLineagesNoQuiescent$mouseID), method = 'bonferroni')
# try(cldList(comparison = dt$comparisons, p.value    = dt$P.adjusted, threshold  = 0.05))


# Table with cell death Proportion in whole tree and in different subtrees, non quiescent lineages
listCellDeath <- data.frame() 
for(i in seq_along(nonQuiescentID)){
  myTree <- dataMax[dataMax$Clone_uniqueID == nonQuiescentID[i],]
  listCellDeathTemp <- data.frame(treeID = nonQuiescentID[i], getDeathRateInSubtrees(myTree))
  listCellDeath <- rbind(listCellDeath, listCellDeathTemp)
}
listCellDeath$rank <- as.numeric(listCellDeath$rank)

myColours2 <- c(add.alpha("blue",0.5), add.alpha("red",0.5), add.alpha("green",0.5), add.alpha("orange",0.5))
my.settings <- list(superpose.symbol=list(col=myColours2)) 

bgColors <- c("black")
txtColors <- c("white")

# Create a function to be passes to "strip=" argument of xyplot
myStripStyle <- function(which.panel, factor.levels, par.strip.text,
                         custBgCol=par.strip.text$custBgCol,
                         custTxtCol=par.strip.text$custTxtCol,...)     {
  panel.rect(0, 0, 1, 1,
             col = custBgCol,
             border = 1)
  panel.text(x = 0.5, y = 0.5,
             font=0.3,
             lab = factor.levels[which.panel],
             col = custTxtCol)
}

with(listCellDeath, xyplot(deathRate ~ rank | treeID, groups = motherRoot, pch = 16, col = myColours2, cex = 0.7, 
                           xlab = "subtree rank (colors indicate common original root for subtrees)", ylab = "death Proportion (% leaf cells)",
                           par.settings = my.settings,  jitter.x = T, jitter.y = T,
                           par.strip.text=list(custBgCol=bgColors,custTxtCol=txtColors, font=2,cex= 0.8),
                           strip=myStripStyle
)
) 
# par.settings = my.settings, par.strip.text=list(col="white", font=2,cex= 0.8),


# To be able to compare the variability of cell death Proportions in subtrees compared to the overall tree, 
# I computed the weighted standard deviation (wSD) of the subtrees at rank 2 and 3 weighted by the number of leaves of the subtree 
# (because trees are mostly asymmetric, a subtree with only 2 leaf cells has a lower influence compared to a subtree with 10 leaf cells for instance). 
# If a wSD of rank 2 subtrees is 10% with an overall death Proportion of 60%, it means that the death Proportions of the subtrees at rank 2 are ~70% and ~50% respectively. 
# The wSD is 0 if the death Proportions of rank 2 subtrees is exactly equal to the overal death Proportion of the lineage. 
# The wSD is a measure of the variability of death Proportions compared to the global cell death of the whole lineage. 


## Weighted Standard deviation
g <- function(y){sqrt(wtd.var(y[,1], y[,2]))}


#############################################
## Using all rank2 and rank3
#############################################

# varDeathRate <- with(listCellDeath, summarize(cbind(deathRate, LeafNumber), llist(treeID, rank), g, stat.name='w.SD'))
# NumberMaxLeaf <- with(listCellDeath, tapply(LeafNumber, INDEX = list(treeID, rank), FUN = sum ))[,1]
# 
# listCellDeath$maxLeaves <- NumberMaxLeaf[match(as.character(listCellDeath$treeID), names(NumberMaxLeaf))]
# listCellDeath$propLeaves <- listCellDeath$LeafNumber/listCellDeath$maxLeaves
# 
# varDeathRate <- varDeathRate[varDeathRate$rank > 1, ]
# varDeathRateRank2 <- varDeathRate[varDeathRate$rank == 2,]
# propLeavesRank2 <- listCellDeath[listCellDeath$rank == 2,]
# 
# #with(varDeathRate,densityplot(~w.SD, groups = rank))
# hist(varDeathRate$w.SD[varDeathRate$rank ==3], col = add.alpha("red",0.5), breaks = seq(0, max(varDeathRate$w.SD+5), by = 5), freq = T, ylim = c(0,15), main = "Histogram of weighted Standard deviation \nof cell death rate in the subtrees at rank 2 and 3", cex.main = 0.9, ylab = "Number of lineages", xlab = "weighted Standard Deviation of subtree death rates (%)")
# hist(varDeathRate$w.SD[varDeathRate$rank ==2], col = add.alpha("blue",0.5), add = T, freq = T)
# legend("topright", legend = c("Subtrees rank 2", "Subtrees rank 3"), fill = c(add.alpha("blue",0.5), add.alpha("red",0.5)))
# 
# wilcoxNonpaired <- wilcox.test(varDeathRate$w.SD[varDeathRate$rank ==3], varDeathRate$w.SD[varDeathRate$rank ==2])
# 
# beeswarm(varDeathRate$w.SD ~ varDeathRate$rank, ylab = "weighted Standard deviation of sublineages", xlab = "Rank", main = paste("Variation of death rate between sublineages at different ranks\n Wilcoxon test pvalue = ", round(wilcoxNonpaired$p.value,5)), cex.main = 0.7)
# bxplot(varDeathRate$w.SD ~ varDeathRate$rank, add = T)

#############################################
## If takes only trees with exactly 2 subtrees at rank 2 and 4 or less rank 3 subtrees (remove trees with too much uncertainties)
#############################################

tab <- table(listCellDeath$rank, listCellDeath$treeID)
treeId2Rank2 <- dimnames(tab)[[2]][tab[2,] == 2]
listCellDeathCertainRank2 <- listCellDeath[listCellDeath$rank == 2 & listCellDeath$treeID %in% treeId2Rank2,]

treeId2Rank3 <- dimnames(tab)[[2]][tab[3,] <= 4]
listCellDeathCertainRank3 <- listCellDeath[listCellDeath$rank == 3 & listCellDeath$treeID %in% treeId2Rank2 ,] # taking all the rank3 subtrees as long as the rank2 are a pair

listCellDeathCertain <- rbind(listCellDeath[listCellDeath$rank == 1 ,], listCellDeathCertainRank2, listCellDeathCertainRank3)


varDeathRateCertain <- with(listCellDeathCertain, summarize(cbind(deathRate, LeafNumber), llist(treeID, rank), g, stat.name='w.SD'))
NumberMaxLeafCertain <- with(listCellDeathCertain, tapply(LeafNumber, INDEX = list(treeID, rank), FUN = sum ))[,1]

listCellDeathCertain$maxLeaves <- NumberMaxLeafCertain[match(as.character(listCellDeathCertain$treeID), names(NumberMaxLeafCertain))]
listCellDeathCertain$propLeaves <- listCellDeathCertain$LeafNumber/listCellDeathCertain$maxLeaves

varDeathRateCertain <- varDeathRateCertain[varDeathRateCertain$rank > 1, ]
varDeathRateRank2Certain <- varDeathRateCertain[varDeathRateCertain$rank == 2,]
propLeavesRank2Certain <- listCellDeathCertain[listCellDeathCertain$rank == 2,]

# With(varDeathRate,densityplot(~w.SD, groups = rank))
hist(varDeathRateCertain$w.SD[varDeathRateCertain$rank ==3], col = add.alpha("red",0.5), breaks = seq(0, max(varDeathRateCertain$w.SD+5), by = 5), freq = T, ylim = c(0,15), main = "Histogram of weighted Standard deviation \nof cell death Proportion in the subtrees at rank 2 and 3", cex.main = 0.9, ylab = "Number of lineages", xlab = "weighted Standard Deviation of subtree death Proportions (%)")
hist(varDeathRateCertain$w.SD[varDeathRateCertain$rank ==2], col = add.alpha("blue",0.5), add = T, freq = T)
legend("topright", legend = c("Subtrees rank 2", "Subtrees rank 3"), fill = c(add.alpha("blue",0.5), add.alpha("red",0.5)))

wilcoxNonpairedCertain <- wilcox.test(varDeathRateCertain$w.SD[varDeathRateCertain$rank ==3], varDeathRateCertain$w.SD[varDeathRateCertain$rank ==2], exact = F)

beeswarm(varDeathRateCertain$w.SD ~ varDeathRateCertain$rank, ylab = "weighted Standard deviation of sublineages", xlab = "Rank", main = paste("Variation of death Proportion between sublineages at different ranks\n Rank2 n= ",length(varDeathRateCertain$w.SD[varDeathRateCertain$rank ==2]),"Rank 3 n = ",length(varDeathRateCertain$w.SD[varDeathRateCertain$rank ==3]),"Wilcoxon test pvalue = ", round(wilcoxNonpairedCertain$p.value,5)), cex.main = 0.7)
bxplot(varDeathRateCertain$w.SD ~ varDeathRateCertain$rank, add = T)


###########################################
# Weighted standard deviation branch wise for rank 3
##########################################

varDeathRatePair <- with(listCellDeathCertain, summarize(cbind(deathRate, LeafNumber), llist(treeID, rank, motherRoot), g, stat.name='w.SD'))
varDeathRatePair <- varDeathRatePair[varDeathRatePair$rank >1, ]

hist(varDeathRatePair$w.SD[varDeathRatePair$rank ==3], col = add.alpha("red",0.5), breaks = seq(0, max(varDeathRatePair$w.SD+5), by = 5), freq = T, ylim = c(0,25), main = "Histogram of weighted Standard deviation \nof cell death Proportion in the subtrees at rank 2 and 3 (pair wise)", cex.main = 0.9, ylab = "Number of lineages", xlab = "weighted Standard Deviation of subtree death Proportions (%)")
hist(varDeathRatePair$w.SD[varDeathRatePair$rank ==2], col = add.alpha("blue",0.5), add = T, freq = T)
legend("topright", legend = c("Subtrees rank 2", "Subtrees rank 3"), fill = c(add.alpha("blue",0.5), add.alpha("red",0.5)))

wilcoxPaired <- wilcox.test(varDeathRatePair$w.SD ~ varDeathRatePair$rank, exact =F)

beeswarm(varDeathRatePair$w.SD ~ varDeathRatePair$rank, ylab = "weighted Standard deviation of paired sublineages", xlab = "Rank", main = paste("Variation of death Proportion between paired sublineages at different ranks\n Wilcoxon test pvalue = ", round(wilcoxPaired$p.value,3)), cex.main = 0.7)
bxplot(varDeathRatePair$w.SD ~ varDeathRatePair$rank, add = T)

barPlotBeeswarm(varDeathRatePair[,c(4,2)], ylab = "weighted Standard deviation of paired sublineages", xlab = "Rank", cex.main = 0.7)


# Relationship between overall death Proportion and death Proportions of rank 2 subtrees. 
# For the analysis only trees that have exactly 2 rank2 subtrees are considered and subtree 1 is arbitrarily set as the one with the highest death Proportion. 


### Here are considered only the lineages which have 2 rank2 subtrees

listCellDeath$treeID <- as.character(listCellDeath$treeID)

# To plot the death Proportion of the whole lineage vs the death Proportion of secondary subtrees
onlyPrimaryLineageID <- names(table(listCellDeathCertain$treeID))[as.vector(table(listCellDeathCertain$treeID) == 1)] # ID of lineages with only death Proportion of whole lineage

listDeathRate <- listCellDeath[listCellDeath$rank %in% c(1,2) & !listCellDeath$treeID %in% onlyPrimaryLineageID,] # take only death Proportion of primary lineage and secondary subtrees

# Remove trees with more than 2 secondary trees
tab <- table(listDeathRate$treeID, listDeathRate$rank)
twoSecondaryTreeID <- rownames(tab)[tab[,2]==2]
listDeathRate <- listDeathRate$deathRate[listDeathRate$treeID %in% twoSecondaryTreeID]

# Transform the vector into a matrix with 1st column, whole lineage death rate, 2nd subtree1 death rate, 3rd column subtree2 death rate
listDeathRate <- t(matrix(listDeathRate, ncol = length(listDeathRate)/3, nrow = 3))

# Order the subtrees so that the second column contains the highest death rate of the two
for (i in c(1:nrow(listDeathRate))){
  listDeathRate[i,2:3] <- listDeathRate[i,order(listDeathRate[i,2:3],decreasing = TRUE)+1]
}
listDeathRate <- as.data.frame(listDeathRate + 0.001) # will not work if there are lines of 0!
names(listDeathRate) <- c("WholeLineage", "Subtree1", "Subtree2")

#############################################
# dotplot
#############################################
plot(0, 0, pch = 1, xlim = c(0.5,3.5), ylim = c(0,100), ylab = "Cell death Proportion", 
     xlab = "Whole lineage = 1, First subtree = 2, Second subtree = 3", las = 1, main = "Relationship cell death between \nwhole tree and rank 2 subtrees", cex.main = 0.9, cex.lab = 0.9)

thresholdDeathRateVariation <- quantile((listDeathRate[,2]-listDeathRate[,3]))[3]

#  mean(varDeathRate$w.SD[varDeathRate$rank ==2])/100 + 2 * sd(varDeathRate$w.SD[varDeathRate$rank ==2])/100

for (i in seq_along(twoSecondaryTreeID)){
  mycolor <- add.alpha("black",0.5)
  if((listDeathRate[i,2]-listDeathRate[i,3]) > thresholdDeathRateVariation){ mycolor = add.alpha("red",0.8)} #/listDeathRate[i,1]
  points(x= c(1,2,3), y = listDeathRate[i,1:3], type = "o", pch = 18, col = mycolor)
}
legend("bottomleft", legend= c(paste("> ",round(thresholdDeathRateVariation,0),"% difference (median)"), paste("< ",round(thresholdDeathRateVariation,0),"% difference")), fill = c(add.alpha("red",0.5),add.alpha("black",0.5)), cex = 0.7, bty = "n")


beeswarm((listDeathRate[,2]-listDeathRate[,3]), ylab = "Difference of Death Proportion between 2 rank2 subtress", main = "Distribution of cell death variation \nbetween the 2 rank2 subtrees", cex.main = 0.8)
bxplot((listDeathRate[,2]-listDeathRate[,3]), add = T, col= "red")

#############################################
## scatter plot 
#############################################

with(listDeathRate, plot(jitter(Subtree1) ~ jitter(Subtree2), pch = 3, xlim = c(0,100), ylim = c(0,100), cex = 0.8, xlab = "Death Proportion rank 2 subtree2 [%]", ylab = "Death Proportion rank 2 subtree1 [%]", main = paste("Comparison death Proportions of rank 2 subtrees \n trees between the red lines have less \nthan ", round(thresholdDeathRateVariation,0) ,"% difference (mean wSD at rank 2)"), cex.main = 0.8))
abline(a = 0, b = 1, col = "black", lwd = 2, lty = 2)
abline(a =  thresholdDeathRateVariation, b = 1, col = "red", lwd = 2)

####
# Triangle plot
###
# The triangle plot enables to visualize 3 dimensional data, in this case the death rate of the whole lineage, and the 2 rank 2 subtrees. The death rates are transformed into percentages of the sum of the 3 death rates, showing which one accounts the most for the death rate. If death rate is symetrically inherited the contribution of each tree is the same : 0.33 (represented by the red circle). 
# If overall death rate is 55%, rank 2 subtree1 80%, rank 2 subtree 2 30%, the coordinates for the triangle plot will be 55/(55+80+30)=0.33, 80/(55+80+30)=0.48, 30/(55+80+30)/0.18.
# triangle.plot(listDeathRate, labeltriangle = TRUE, addmean = T, scale = FALSE, show.position = FALSE, draw.line = F) #, min3 = c(0.1,0.1,0.1),max3 = c(0.9,0.9,0.9), pch = 10)
# symbols(x = 0, y = 0, circles = c(0.07), col= add.alpha("red", alpha = 0.5), add = TRUE, inches = F, bg = add.alpha("red",alpha = 0.3))
# #draw.circle(x = 0, y = 0, radius =  0.05, col= add.alpha("red",alpha = 0.5), border = add.alpha("red",alpha = 0.5))
# abline(a = 0, b = 1.7, lty =2)
# abline(a = 0, b = -1.7, lty =2)
# abline(h =0, lty =2)


################################################################################################################################
# Cell death Proportion focusing on early cell death
################################################################################################################################

#############################################
## Cell death Proportion only early death
#############################################

# Late cell death are considered as surviving neurons, the threshold early vs late cell death is `r lateDeathThreshold` days.

dataMax$UncertaintyType2 <- as.character(dataMax$UncertaintyType2)
listCellDeathEarly <- data.frame() # table with cell death rate in whole tree and in different subtrees, non quiescent lineages
for(i in seq_along(nonQuiescentID)){
  myTree <- dataMax[dataMax$Clone_uniqueID == nonQuiescentID[i],]
  listCellDeathTemp <- data.frame(treeID = nonQuiescentID[i], getDeathRateInSubtrees(myTree, myDeathType = "early"))
  listCellDeathEarly <- rbind(listCellDeathEarly, listCellDeathTemp)
}
listCellDeathEarly$rank <- as.numeric(listCellDeathEarly$rank)

listCellDeathEarly$motherRoot <- paste(listCellDeathEarly$motherRoot, listCellDeathEarly$treeID, sep = ".")
listCellDeathEarly$root <- paste(listCellDeathEarly$root, listCellDeathEarly$treeID, sep = ".")

lineageDeathRate <- data.frame(deathType = "all", deathRate = listCellDeath$deathRate[listCellDeath$rank == 1])
lineageEarlyDeathRate <- listCellDeathEarly$deathRate[listCellDeathEarly$rank == 1]
lineageDeathRate <- rbind(lineageDeathRate, data.frame(deathType = "early", deathRate = lineageEarlyDeathRate))

meanDeathRateall <- mean(lineageDeathRate$deathRate[lineageDeathRate$deathType =="all"])
semDeathRateall <- sd(lineageDeathRate$deathRate[lineageDeathRate$deathType =="all"])/sqrt(sum(lineageDeathRate$deathType =="all"))

meanDeathRateEarly <- round(mean(lineageDeathRate$deathRate[lineageDeathRate$deathType =="early"]),2)
semDeathRateEarly <- round(sd(lineageDeathRate$deathRate[lineageDeathRate$deathType =="early"])/sqrt(sum(lineageDeathRate$deathType =="early")),2)

# boxplot(lineageDeathRate$deathRate ~ lineageDeathRate$deathType, ylab = "Death Rate [%]")
beeswarm(lineageDeathRate$deathRate ~ lineageDeathRate$deathType, col = "red", corral = "wrap", ylab = "Death rate [%]", xlab = "type of cell death")
bxplot(lineageDeathRate$deathRate ~ lineageDeathRate$deathType, add = T,ylab = "Death Rate [%]")

wilcox.test(lineageDeathRate$deathRate ~ lineageDeathRate$deathType)

hist(lineageEarlyDeathRate, col = "white", ylab = "Number of lineages", main = paste("Early death Proportion \n mean = ", round(mean(lineageEarlyDeathRate),2), "sem =", round(semDeathRateEarly,2)), cex.main = 0.8)
abline(v = mean(lineageEarlyDeathRate), col = "red", lwd = 2)


#shapiro.test(lineageEarlyDeathRate)
#shapiro.test(listCellDeath$deathRate[listCellDeath$rank == 1])


myColours2 <- c(add.alpha("blue",0.5), add.alpha("red",0.5), add.alpha("green",0.5), add.alpha("orange",0.5))
my.settings <- list(superpose.symbol=list(col=myColours2)) 

bgColors <- c("black")
txtColors <- c("white")

# Create a function to be passes to "strip=" argument of xyplot
myStripStyle <- function(which.panel, factor.levels, par.strip.text,
                         custBgCol=par.strip.text$custBgCol,
                         custTxtCol=par.strip.text$custTxtCol,...)     {
  panel.rect(0, 0, 1, 1,
             col = custBgCol,
             border = 1)
  panel.text(x = 0.5, y = 0.5,
             font=0.3,
             lab = factor.levels[which.panel],
             col = custTxtCol)
}

with(listCellDeathEarly, xyplot(deathRate ~ rank | treeID, groups = motherRoot, pch = 16, col = myColours2, cex = 0.7, 
                                xlab = "subtree rank (colors indicate common original root for subtrees)", ylab = "early death Proportion (% leaf cells)",
                                par.settings = my.settings,  jitter.x = T, jitter.y = T,
                                par.strip.text=list(custBgCol=bgColors,custTxtCol=txtColors, font=2,cex= 0.8),
                                strip=myStripStyle
)
) # par.settings = my.settings, par.strip.text=list(col="white", font=2,cex= 0.8),

####################################################
## Weighted Standard deviation
####################################################

g <- function(y){sqrt(wtd.var(y[,1], y[,2]))}

# varDeathRateEarly <- with(listCellDeathEarly, summarize(cbind(deathRate, LeafNumber), llist(treeID, rank), g, stat.name='w.SD'))
# NumberMaxLeafEarly <- with(listCellDeathEarly, tapply(LeafNumber, INDEX = list(treeID, rank), FUN = sum ))[,1]
# 
# listCellDeathEarly$maxLeaves <- NumberMaxLeafEarly[match(as.character(listCellDeathEarly$treeID), names(NumberMaxLeafEarly))]
# listCellDeathEarly$propLeaves <- listCellDeathEarly$LeafNumber/listCellDeathEarly$maxLeaves
# 
# varDeathRateEarly <- varDeathRateEarly[varDeathRateEarly$rank >1, ]
# #varDeathRateEarly <- varDeathRateEarly[varDeathRateEarly$rank == 2,]
# propLeavesRank2Early <- listCellDeathEarly[listCellDeathEarly$rank == 2,]
# 
# #with(varDeathRate,densityplot(~w.SD, groups = rank))
# hist(varDeathRateEarly$w.SD[varDeathRateEarly$rank ==3], col = add.alpha("red",0.5), breaks = seq(0, max(varDeathRateEarly$w.SD)+5, by = 5), freq = T, ylim = c(0,15), main = "Histogram of weighted Standard deviation \nof early cell death rate in the subtrees at rank 2 and 3", cex.main = 0.9, ylab = "Number of lineages", xlab = "weighted Standard Deviation of subtree early death rates (%)")
# hist(varDeathRateEarly$w.SD[varDeathRateEarly$rank ==2], col = add.alpha("blue",0.5), add = T, freq = T)
# legend("topright", legend = c("Subtrees rank 2", "Subtrees rank 3"), fill = c(add.alpha("blue",0.5), add.alpha("red",0.5)))
# 
# wilcoxNonpaired <- wilcox.test(varDeathRateEarly$w.SD[varDeathRateEarly$rank ==3], varDeathRateEarly$w.SD[varDeathRateEarly$rank ==2], exact = F)
# 
# beeswarm(varDeathRateEarly$w.SD ~ varDeathRateEarly$rank, ylab = "weighted Standard deviation of sublineages", xlab = "Rank", main = paste("Variation of death rate between sublineages at different ranks\n Wilcoxon test pvalue = ", round(wilcoxNonpaired$p.value,5)), cex.main = 0.7)
# bxplot(varDeathRateEarly$w.SD ~ varDeathRateEarly$rank, add = T)

####################################################
## Test if we taken only the subtrees where mother have only exactly 2 daughters
####################################################

#numberSubtreePerMother <- table(listCellDeathEarly$motherRoot)
#only2SubtreesID <- names(numberSubtreePerMother)[numberSubtreePerMother == 2]

#listCellDeathCertainTest <- listCellDeathEarly[listCellDeathEarly$motherRoot %in% only2SubtreesID,]

#listCellDeathCertainEarly <- listCellDeathCertainTest


####################################################
## If takes only trees with exactly 2 subtrees at rank 2 a listCellDeathCertainEarly and 4 or less rank 3 subtrees (remove trees with too much uncertainties)
####################################################

tab <- table(listCellDeathEarly$rank, listCellDeathEarly$treeID)
treeId2Rank2 <- dimnames(tab)[[2]][tab[2,] == 2]
listCellDeathCertainRank2Early <- listCellDeathEarly[listCellDeathEarly$rank == 2 & listCellDeathEarly$treeID %in% treeId2Rank2,]

# treeId2Rank3 <- dimnames(tab)[[2]][tab[3,] <= 4]
listCellDeathCertainRank3Early <- listCellDeathEarly[listCellDeath$rank == 3 & listCellDeathEarly$treeID %in% treeId2Rank2 ,] # & listCellDeath$treeID %in% treeId2Rank2

listCellDeathCertainEarly <- rbind(listCellDeathEarly[listCellDeathEarly$rank == 1 ,], listCellDeathCertainRank2Early, listCellDeathCertainRank3Early)


varDeathRateEarlyCertain <- with(listCellDeathCertainEarly, summarize(cbind(deathRate, LeafNumber), llist(treeID, rank), g, stat.name='w.SD'))
NumberMaxLeafEarlyCertain <- with(listCellDeathCertainEarly, tapply(LeafNumber, INDEX = list(treeID, rank), FUN = sum ))[,1]

listCellDeathCertainEarly$maxLeaves <- NumberMaxLeafEarlyCertain[match(as.character(listCellDeathCertainEarly$treeID), names(NumberMaxLeafEarlyCertain))]
listCellDeathCertainEarly$propLeaves <- listCellDeathCertainEarly$LeafNumber/listCellDeathCertainEarly$maxLeaves

varDeathRateEarlyCertain <- varDeathRateEarlyCertain[varDeathRateEarlyCertain$rank >1, ]
#varDeathRateEarly <- varDeathRateEarly[varDeathRateEarly$rank == 2,]
propLeavesRank2EarlyCertain <- listCellDeathCertainEarly[listCellDeathCertainEarly$rank == 2,]

# with(varDeathRate,densityplot(~w.SD, groups = rank))
hist(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==3], col = add.alpha("red",0.5), breaks = seq(0, max(varDeathRateEarlyCertain$w.SD)+5, by = 5), freq = T, ylim = c(0,15), main = "Histogram of weighted Standard deviation \nof early cell death Proportion in the subtrees at rank 2 and 3", cex.main = 0.9, ylab = "Number of lineages", xlab = "weighted Standard Deviation of subtree early death Proportions (%)")
hist(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==2], col = add.alpha("blue",0.5), add = T, freq = T)
legend("topright", legend = c("Subtrees rank 2", "Subtrees rank 3"), fill = c(add.alpha("blue",0.5), add.alpha("red",0.5)))

wilcoxNonpairedCertain <- wilcox.test(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==3], varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==2], exact = F)

meanVariationRank2 <- round(mean(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==2]),2)
semVariationRank2 <- round(sd(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==2])/ length(sqrt(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==2])),2)

meanVariationRank3 <- round(mean(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==3]),2)
semVariationRank3 <- round(sd(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==3])/ length(sqrt(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==3])),2)


beeswarm(varDeathRateEarlyCertain$w.SD ~ varDeathRateEarlyCertain$rank, ylab = "weighted Standard deviation of sublineages", xlab = "Rank", main = paste("Variation of early death Proportion between sublineages at different ranks\n Rank2 n= ",length(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==2]),"Rank 3 n = ",length(varDeathRateEarlyCertain$w.SD[varDeathRateEarlyCertain$rank ==3]),"Wilcoxon test pvalue = ", round(wilcoxNonpairedCertain$p.value,5)), cex.main = 0.7)
bxplot(varDeathRateEarlyCertain$w.SD ~ varDeathRateEarlyCertain$rank, add = T)


barPlotBeeswarm(varDeathRateEarlyCertain[,c(3,2)], ylab = "weighted Standard deviation of cell death sublineages", xlab = "rank of sublineages")
text(x = 1, y = 50, paste("Wilcoxon test pvalue = ", round(wilcoxNonpairedCertain$p.value,5)))

####################################################
## Weighted standard deviation branch wise for rank 3
###################################################

varDeathRatePairEarly <- with(listCellDeathCertainEarly , summarize(cbind(deathRate, LeafNumber), llist(treeID, rank, motherRoot), g, stat.name='w.SD'))

varDeathRatePairEarly <- varDeathRatePairEarly[varDeathRatePairEarly$rank >1, ] # interested in differences for rank > 1
hist(varDeathRatePairEarly$w.SD[varDeathRatePairEarly$rank ==3], col = add.alpha("red",0.5), breaks = seq(0, max(varDeathRatePairEarly$w.SD+5), by = 5), freq = T, ylim = c(0,25), main = "Histogram of weighted Standard deviation \nof early cell death Proportion in the subtrees at rank 2 and 3 (pair wise)", cex.main = 0.9, ylab = "Number of lineages", xlab = "weighted Standard Deviation of subtree early death Proportions (%)")
hist(varDeathRatePairEarly$w.SD[varDeathRatePairEarly$rank ==2], col = add.alpha("blue",0.5), add = T, freq = T)
legend("topright", legend = c("Subtrees rank 2", "Subtrees rank 3"), fill = c(add.alpha("blue",0.5), add.alpha("red",0.5)))

wilcoxNonpaired <- wilcox.test(varDeathRatePairEarly$w.SD[varDeathRatePairEarly$rank ==3], varDeathRatePairEarly$w.SD[varDeathRatePairEarly$rank ==2])
beeswarm(varDeathRatePairEarly$w.SD ~ varDeathRatePairEarly$rank, ylab = "weighted Standard deviation of sublineages", xlab = "Rank", main = paste("Variation of early death Proportion between paired sublineages at different ranks\n Wilcoxon test pvalue = ", round(wilcoxNonpaired$p.value,5)), cex.main = 0.7)
bxplot(varDeathRatePairEarly$w.SD ~ varDeathRatePairEarly$rank, add = T)



### Here are considered only the lineages which have 2 rank2 subtrees

listCellDeathEarly$treeID <- as.character(listCellDeathEarly$treeID)

# To plot the death rate of the whole lineage vs the death rate of secondary subtrees
onlyPrimaryLineageID <- names(table(listCellDeathEarly$treeID))[as.vector(table(listCellDeathEarly$treeID) == 1)] # ID of lineages with only death rate of whole lineage

listCellDeathEarly2 <- listCellDeathEarly[listCellDeathEarly$rank %in% c(1,2) & !listCellDeathEarly$treeID %in% onlyPrimaryLineageID,] # take only death rate of primary lineage and secondary subtrees


# Remove trees with more than 2 secondary trees
tab <- table(listCellDeathEarly2$treeID, listCellDeathEarly2$rank)
twoSecondaryTreeID <- rownames(tab)[tab[,2]==2]
listDeathRateEarly <- listCellDeathEarly2$deathRate[listCellDeathEarly2$treeID %in% twoSecondaryTreeID]

# Transform the vector into a matrix with 1st column, whole lineage death rate, 2nd subtree1 death rate, 3rd column subtree2 death rate
listDeathRateEarly <- t(matrix(listDeathRateEarly, ncol = length(listDeathRateEarly)/3, nrow = 3))

# Order the subtrees so that the second column contains the highest death rate of the two
for (i in c(1:nrow(listDeathRateEarly))){
  listDeathRateEarly[i,2:3] <- listDeathRateEarly[i,order(listDeathRateEarly[i,2:3],decreasing = TRUE)+1]
}
listDeathRateEarly <- as.data.frame(listDeathRateEarly + 0.001) # will not work if there are lines of 0!
names(listDeathRateEarly) <- c("WholeLineage", "Subtree1", "Subtree2")

####################################################
# dotplot
####################################################

plot(0, 0, pch = 1, xlim = c(0.5,3.5), ylim = c(0,100), ylab = "Cell death Proportion", 
     xlab = "Whole lineage = 1, First subtree = 2, Second subtree = 3", las = 1, main = "Relationship early death Proportion between \nwhole tree and rank 2 subtrees", cex.main = 0.9, cex.lab = 0.9)

thresholdDeathRateVariation <- quantile((listDeathRateEarly[,2]-listDeathRateEarly[,3]))[3]
#  mean(varDeathRate$w.SD[varDeathRate$rank ==2])/100 + 2 * sd(varDeathRate$w.SD[varDeathRate$rank ==2])/100

for (i in seq_along(twoSecondaryTreeID)){
  mycolor <- add.alpha("black",0.5)
  if((listDeathRateEarly[i,2]-listDeathRateEarly[i,3]) > thresholdDeathRateVariation){ mycolor = add.alpha("red",0.8)} #/listDeathRate[i,1]
  points(x= c(1,2,3), y = listDeathRateEarly[i,1:3], type = "o", pch = 18, col = mycolor)
}
legend("bottomleft", legend= c(paste("> ",round(thresholdDeathRateVariation,0),"% difference (median)"), paste("< ",round(thresholdDeathRateVariation,0),"% difference")), fill = c(add.alpha("red",0.5),add.alpha("black",0.5)), cex = 0.7, bty = "n")


beeswarm((listDeathRateEarly[,2]-listDeathRateEarly[,3]), ylab = "Difference of Early Death Proportion between 2 rank2 subtress", main = "Distribution of cell death variation \nbetween the 2 rank2 subtrees", cex.main = 0.8)
bxplot((listDeathRateEarly[,2]-listDeathRateEarly[,3]), add = T, col= "red")


####################################################
## scatter plot
####################################################

with(listDeathRateEarly, plot(jitter(Subtree1) ~ jitter(Subtree2), pch = 3, xlim = c(0,100), ylim = c(0,100), cex = 0.8, xlab = "Early Death Proportion rank 2 subtree2 [%]", ylab = "Early death Proportion rank 2 subtree1 [%]", main = paste("Comparison early death Proportions of rank 2 subtrees \n trees between the red lines have less \nthan ", round(thresholdDeathRateVariation,0) ,"% difference (mean wSD at rank 2)"), cex.main = 0.8))
abline(a = 0, b = 1, col = "black", lwd = 2, lty = 2)
abline(a =  thresholdDeathRateVariation, b = 1, col = "red", lwd = 2)


treeVisualizationDeathRate(listCellDeath, saveGraph = TRUE, myFileName = "DeathRateTreeVisualization.pdf", sizeAmplif =10)
treeVisualizationDeathRate(listCellDeathEarly, saveGraph = T, myFileName = "EarlyDeathRateTreeVisualization_CellTypes.pdf", sizeAmplif =10)


#############################################
## Only late cell death
#############################################


dataMax$cellDeathLate <- rep(1,n = length(dataMax$cellDeathLate))

for (i in seq_along(dataMax$cellDeath)){
  if(dataMax$cellDeath[i] == 0 & dataMax$observationTime[i] > 7){
    dataMax$cellDeathLate[i] <- 0
  }
}
numberLineages4LateDeath <- sum(table(dataMax$cellDeathLate, dataMax$Clone_uniqueID)[1,] > 3)

hist(table(dataMax$cellDeathLate, dataMax$Clone_uniqueID)[1,], ylab = "Number of lineages", xlab = "Number of late cell death per lineage", main = paste("Distribution of number of late cell death per lineage,\n",numberLineages4LateDeath,"lineages with 4 or more late cell death"), cex.main = 0.8)


