################################################################################################################################
# Lineage asymmetry
#############################################

# The asymmetry index of the subtrees rank2 indicates the proportion difference of generated leaves between subtrees. 
# The asymmetry index is between 2 extrems 0 = no asymmetry (50% of leaves in each subtree), 100 = completely asymmetric (not feasible : 100% on one side, 0 on the other).

#############################################
## Subtree asymmetry
#############################################
#listCellDeathRank2 <- listCellDeathCertain[listCellDeathCertain$rank == 2 & listCellDeathCertain$treeID %in% twoSecondaryTreeID,]
#varDeathRateRank2Only2Subtrees <- varDeathRateCertain[varDeathRateCertain$rank == 2 & varDeathRateCertain$treeID %in% twoSecondaryTreeID,]
#numberLeavesRank2 <- with(listCellDeathRank2, tapply(LeafNumber, INDEX = list(treeID), sum))
# 
# propSubtreeRank2 <- with(listCellDeathRank2, tapply(propLeaves, INDEX = list(treeID), max))
# asymIndex <- (propSubtreeRank2-0.5)*2 # transformation of the proportion of leaves of bigger subtree to get an index between 2 extrems 0 = no asymmetry (50% of leaves in each subtree), 100 = completely asymmetric (not feasible : 100% on one side, 0 on the other).
# 
# hist(asymIndex, xlab = "% Asymmetry Index between rank2 subtrees", ylab = "Number of lineages", main = "Asymmetry of sublineages at Rank 2")
# 
# plot(varDeathRateRank2Only2Subtrees$w.SD ~ propSubtreeRank2, main = "Sublineages rank 2 :\nRelationship cell death rate and lineage asymmetry", cex.main = 0.8, xlab = "Proportion of the bigger subtree rank2", ylab = "weighted Standard Deviation of subtree death rates (%)")
# cor.test(varDeathRateRank2Only2Subtrees$w.SD, propSubtreeRank2, method = "spearman")
# 
# #plot(varDeathRateRank2Only2Subtrees$w.SD ~ numberLeavesRank2, main = "Sublineages rank 2 :\nRelationship cell death rate and lineage size", cex.main = 0.8, xlab = "Number of leaves in the lineages", ylab = "weighted Standard Deviation of subtree death rates (%)")
# #cor.test(varDeathRateRank2Only2Subtrees$w.SD, numberLeavesRank2, method = "spearman")
# 
# #plot(propSubtreeRank2 ~ numberLeavesRank2, main = "Sublineages rank 2 :\nRelationship sublineage asymmetry and lineage size", cex.main = 0.8, xlab = "Number of leaves in the lineages", ylab = "Sublineage Asymmetry")
# #cor.test(propSubtreeRank2, numberLeavesRank2, method = "spearman")
# 
#############################################
# ## Subtree asymmetry vs early cell death
#############################################
# listEarlyCellDeathRank2 <- listCellDeathEarly[listCellDeathEarly$rank == 2 & listCellDeathEarly$treeID %in% twoSecondaryTreeID,]
# varEarlyDeathRateRank2Only2Subtrees <- varDeathRateEarlyCertain[varDeathRateEarlyCertain$rank == 2 & varDeathRateEarlyCertain$treeID %in% twoSecondaryTreeID,]
# numberLeavesRank2 <- with(listEarlyCellDeathRank2, tapply(LeafNumber, INDEX = list(treeID), sum))

#myCor <- cor.test(varEarlyDeathRateRank2Only2Subtrees$w.SD, numberLeavesRank2, method = "spearman",exact = F)

#plot(varEarlyDeathRateRank2Only2Subtrees$w.SD ~ numberLeavesRank2, main = paste("Sublineages rank 2 :\nRelationship early cell death rate and lineage size "), cex.main = 0.8, xlab = "Number of leaves in the lineages", ylab = "weighted SD of subtree early death rates [%]") #cor pvalue = ,round(myCor$p.value,2)


################################################################################################################################
# Progeny of Radial and Non Radial progenitors
################################################################################################################################

#############################################
# Cell type 2:
## - R = Type 1 : radial progenitor
## - NR = Type2/Neuroblast that divides (Non radial progenitor)
## - N = Neuron / Neuroblast that does not divide
## - A = Astrocyte
## - Cell death

# For this analysis I selected only divisions for which cell type is certain for mother and daughters and both lineage transitions towards sisters are also certain. 
# There are 195 divisions for RN, 68 for R.
#############################################

divisionTable <- transformTreeIntoDivision(decodedData, cellType = "simplified no death")  

filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyType = c("1"), uncertaintyLineage = c("1"))

resultsRank <- getTransitionDiagram(filteredDivisionTable, cellType = "simplified no death", rank = "successive", detailedSymmetry = "TRUE")

# Try: 
# filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyType = c("1", "0-5"), uncertaintyLineage = c("1", "0-5"))
# resultsRank <- getTransitionDiagram(divisionTable, cellType = "simplified no death", rank = "successive", detailedSymmetry = "TRUE")

symRank <- resultsRank[[1]]
nRank <-  resultsRank[[2]]
namesRank <-  resultsRank[[3]]
namesDiv <- resultsRank[[4]]
namesCellType <- resultsRank[[5]]

my.hist <- barplot(symRank[1:4,]*100, names.arg = namesRank, main = paste("R: Proportion of division type according to division rank"),
                   ylab = "% division type", cex.main = 0.8, xlab = "rank of division", col = rainbow(4), beside = FALSE)
text(x = as.matrix(my.hist), y = 5,  labels = colSums(nRank[1:4,]), cex = .8)

my.hist <- barplot(symRank[5:8,]*100, names.arg = namesRank, main = paste("NR: Proportion of division type according to division rank"),
                   ylab = "% division type", cex.main = 0.8, xlab = "rank of division", col = rainbow(4), beside = FALSE)
text(x = as.matrix(my.hist), y = 5,  labels = colSums(nRank[5:8,]), cex = .8)
legend("right",legend = namesDiv, fill = rainbow(4), cex = 0.8)

# Statistical analysis

#### R

RDivTypeProportions <- nRank[1:4,]

cat("Regular rank : R rank 1 vs 2")
print(fisher.test(cbind(RDivTypeProportions[,2],RDivTypeProportions[,3])))

cat("Regular rank : R rank 2 vs 3")
print(fisher.test(cbind(RDivTypeProportions[,3],RDivTypeProportions[,4])))

## Sym vs Asym div vs rank
RDivTypeAsymSym <- rbind(colSums(RDivTypeProportions[1:2,]), colSums(RDivTypeProportions[3:4,]))
fisher.test(RDivTypeAsymSym[,2:4])


### NR
NRDivTypeProportions <- nRank[5:8,]

cat("Regular rank : NR rank 2 vs 3")
print(fisher.test(cbind(NRDivTypeProportions[,3],NRDivTypeProportions[,4])))

cat("Regular rank : NR rank 3 vs 4")
print(fisher.test(cbind(NRDivTypeProportions[,4],NRDivTypeProportions[,5])))

cat("Regular rank : NR rank 4 vs 5")
print(fisher.test(cbind(NRDivTypeProportions[,5],NRDivTypeProportions[,6])))

## Sym vs Asym div vs rank
NRDivTypeAsymSym <- rbind(colSums(NRDivTypeProportions[1:2,]), colSums(NRDivTypeProportions[3:4,]))
fisher.test(NRDivTypeAsymSym[,3:8])

#############################################
## NR rank
#############################################
# Here the NR rank is taken : 
# for R cell nothing is changes, for NR cells rank 1 corresponds to the first apparition of an RN cell generated by a R cell in the branch.

resultsRank <- getTransitionDiagram(filteredDivisionTable, cellType = "simplified no death", rank = "NR", detailedSymmetry = "TRUE")
# Try:
# resultsRank <- getTransitionDiagram(divisionTable, cellType = "simplified no death", rank = "NR", detailedSymmetry = "TRUE")


symRank <- resultsRank[[1]]
nRank <-  resultsRank[[2]]
namesRank <-  resultsRank[[3]]
namesDiv <- resultsRank[[4]]
namesCellType <- resultsRank[[5]]

my.hist <- barplot(symRank[5:8,]*100, names.arg = namesRank, main = paste("NR: Proportion of division type according to NR division rank"),
                   ylab = "% division type", cex.main = 0.8, xlab = "rank of division", col = rainbow(4), beside = FALSE)
text(x = as.matrix(my.hist), y = 5,  labels = colSums(nRank[5:8,]), cex = .8)

# Statistical analysis

### NR
NRDivTypeProportions <- nRank[5:8,]

cat("NR rank : NR rank 1 vs 2")
print(fisher.test(cbind(NRDivTypeProportions[,2],NRDivTypeProportions[,3])))

cat("NR rank : NR rank 2 vs 3")
print(fisher.test(cbind(NRDivTypeProportions[,3],NRDivTypeProportions[,4])))

cat("NR rank : NR rank 3 vs 4")
print(fisher.test(cbind(NRDivTypeProportions[,4],NRDivTypeProportions[,5])))

cat("NR rank : NR rank 4 vs 5")
print(fisher.test(cbind(NRDivTypeProportions[,5],NRDivTypeProportions[,6])))
# 

## Sym vs Asym div vs rank
NRDivTypeAsymSym <- rbind(colSums(NRDivTypeProportions[1:2,]), colSums(NRDivTypeProportions[3:4,]))
#fisher.test(NRDivTypeAsymSym[,2:7])


################################################################################################################################
# Source of Neurons
################################################################################################################################

neuronCertainMotherID <- dataMax$Mother_uniqueID[dataMax$CellType3 == 4 & dataMax$UncertaintyType3 == 1]

neuronMotherType <- vector(mode = "numeric", length= length(neuronCertainMotherID))

for(i in seq_along(neuronCertainMotherID)){
  if(dataMax$UncertaintyType2[dataMax$Cell_uniqueID == neuronCertainMotherID[i]] == 1){
    neuronMotherType[i] <- dataMax$CellType2[dataMax$Cell_uniqueID == neuronCertainMotherID[i]]
  } else {neuronMotherType[i] <-  NA}
}


tab <- table(neuronMotherType)
numberNeuron <- length(neuronCertainMotherID)
propTab <- round(prop.table(table(neuronMotherType)) * 100)
propNRneuronMother <- propTab[2]


################################################################################################################################
# R last action
################################################################################################################################


RdeathID <- dataMax$Cell_uniqueID[dataMax$CellType3 == 1 & dataMax$cellDeath == 0]
lastDivisionR <- filteredDivisionTable[filteredDivisionTable$MotherCellType == 1 & filteredDivisionTable$D1CellType != 1 &  filteredDivisionTable$D2CellType != 1,]

percentDeath <- round(length(RdeathID)/ (length(RdeathID) + length(lastDivisionR$MotherID)),2)*100

getTransitionDiagram(lastDivisionR, cellType = "simplified no death", rank = NULL, detailedSymmetry = "FALSE")

divisionTable <- transformTreeIntoDivision(decodedData, cellType = "simplified")  
filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyType = c("1"), uncertaintyLineage = c("1"))
lastDivisionR <- filteredDivisionTable[filteredDivisionTable$MotherCellType == 1 & filteredDivisionTable$D1CellType != 1 &  filteredDivisionTable$D2CellType != 1,]

#pdf(file = "S3_RLastActionWithDeath.pdf",width = 3, height = 3)
getTransitionDiagram(lastDivisionR, cellType = "simplified", rank = NULL, detailedSymmetry = "FALSE")



################################################################################################################################
# Clone activity duration
################################################################################################################################

dataMaxProg <- dataMax[dataMax$Cell_uniqueID %in% dataMax$Mother_uniqueID,]

getActivityDuration <- function(time){
  activityDuration <- max(time)-min(time)
  activityDuration
}

activityDurationClone <- with(dataMaxProg, tapply(DPI, INDEX = Clone_uniqueID, FUN = getActivityDuration))
activityDurationCloneNonZero <- activityDurationClone[activityDurationClone > 0]

hist(activityDurationCloneNonZero, xlab = "Activity duration [Days]", ylab = "Number of lineages", 
     main = paste("Duration of clone activity \n quiescent clones, one division clones excluded  \n mean = ", 
     round(mean(activityDurationCloneNonZero),2), " sem = ", round(sd(activityDurationCloneNonZero)/sqrt(length(activityDurationCloneNonZero)),1), "n =", length(activityDurationCloneNonZero)), cex.main = 0.8 )
abline(v = mean(activityDurationCloneNonZero), col = "red", lwd = 2)


beeswarm(activityDurationCloneNonZero, ylab = "Activity duration [Days]", xlab = "Number of lineages", 
         main = paste("Duration of clone activity \n quiescent clones, one division clones excluded  \n mean = ", 
         round(mean(activityDurationCloneNonZero),2), " sem = ", round(sd(activityDurationCloneNonZero)/sqrt(length(activityDurationCloneNonZero)),1), "n =", length(activityDurationCloneNonZero)), cex.main = 0.8 )
bxplot(activityDurationCloneNonZero, add = T, col = "red")

# pdf(file = "S3_CloneActivity.pdf", width = 3,height = 4)
barPlotBeeswarm(activityDurationCloneNonZero, ylab = "Activity duration [Days]", xlab = "Number of lineages", cex.main = 0.8)
text(1, 40, labels = "NB: 1 lineage only 1 division\n = activity duration is 0 --> removed", cex =0.5)
# dev.off()


# Yicheng
#############################################
# Try with 'Zero': 
# hist(activityDurationClone, xlab = "Activity duration [Days]", ylab = "Number of lineages", 
#      main = paste("Duration of clone activity  \n mean = ", 
#      round(mean(activityDurationClone),2), " sem = ", round(sd(activityDurationClone)/sqrt(length(activityDurationClone)),1), "n =", length(activityDurationClone)), cex.main = 0.8 )
# abline(v = mean(activityDurationClone), col = "red", lwd = 2)

# beeswarm(activityDurationClone, ylab = "Activity duration [Days]", xlab = "Number of lineages", 
#         main = paste("Duration of clone activity  \n mean = ", 
#         round(mean(activityDurationClone),2), " sem = ", round(sd(activityDurationClone)/sqrt(length(activityDurationClone)),1), "n =", length(activityDurationClone)), cex.main = 0.8 )
# bxplot(activityDurationClone, add = T, col = "red")

# barPlotBeeswarm(activityDurationClone, ylab = "Activity duration [Days]", xlab = "Number of lineages", cex.main = 0.8)
#############################################


################################################################################################################################
# R Self-renewal duration
################################################################################################################################

#############################################
# Time between first division of the R cell and the exhaustion of R cells, either depleting division or cell death. 
# But in a few cases the R cell is still present at the end of the recording meaning that the self-renewing duration can be bigger than the activity time.
#############################################

barplot(table(maxRRank), ylab = "Number of lineages", xlab = "Number of successive R divisions", main = "Lineages with R cells self-renewal", col = "white", cex.lab = 0.7, cex.main = 0.7)

# Maximum rank of a R cell in every lineage
maxRRank <- with(dataMaxR, tapply(Rank, as.factor(Clone_uniqueID), max)) 

SelfRenawlDuration <- c(1:length(maxRRank))

# Need to substract the time before the first division
for (l in seq_along(maxRRank)){
  dataMaxRLineage <- dataMaxR[dataMaxR$Clone_uniqueID == names(maxRRank)[l],]
  SelfRenawlDuration[l] <- max(dataMaxRLineage$DPI[dataMaxRLineage$Rank == maxRRank[l]] - dataMax$DPI[dataMax$Clone_uniqueID == names(maxRRank)[l] & dataMax$Rank == 1]) 
}
SelfRenawlDurationNonZero <- SelfRenawlDuration[SelfRenawlDuration>0]

hist(SelfRenawlDurationNonZero, ylab = "Number of lineages", xlab = "Self-renewal duration from 1st division [Days]", 
     main = paste("Duration of R cell self-renewal \nmean = ", round(mean(SelfRenawlDurationNonZero),2)," SD = ", 
     round(sd(SelfRenawlDurationNonZero),1)), cex.lab = 0.7, cex.main = 0.7)
abline(v = mean(SelfRenawlDurationNonZero), col = "red", lwd = 2)

beeswarm(SelfRenawlDurationNonZero, xlab = "Number of lineages", ylab = "Self-renewal duration from 1st division [Days]", 
         main = paste("Duration of R cell self-renewal \nmean = ",round(mean(SelfRenawlDurationNonZero),2)," SD = ", 
         round(sd(SelfRenawlDurationNonZero),1)), cex.lab = 0.7, cex.main = 0.7)
bxplot(SelfRenawlDurationNonZero, add = T, col = "red")

# pdf(file = "F2C_RSelfrenewal.pdf", width = 3,height = 4)
barPlotBeeswarm(SelfRenawlDurationNonZero,xlab = "Number of lineages", ylab = "Self-renewal duration from 1st division [Days]", cex.lab = 0.7, cex.main = 0.7)
text(1, 30, labels = "NB: 3 lineage only 1 R division\n = R selfrenewal duration is 0 --> removed", cex =0.5)
# dev.off()

# Yicheng
#############################################
# Try with 'Zero': 
# hist(SelfRenawlDuration, ylab = "Number of lineages", xlab = "Self-renewal duration [Days]", 
#      main = paste("Duration of R cell self-renewal \nmean = ", round(mean(SelfRenawlDuration),2)," SD = ", 
#      round(sd(SelfRenawlDuration),1)), cex.lab = 0.7, cex.main = 0.7)
# abline(v = mean(SelfRenawlDuration), col = "red", lwd = 2)

# beeswarm(SelfRenawlDuration, xlab = "Number of lineages", ylab = "Self-renewal duration [Days]", 
#          main = paste("Duration of R cell self-renewal \nmean = ",round(mean(SelfRenawlDuration),2)," SD = ", 
#          round(sd(SelfRenawlDuration),1)), cex.lab = 0.7, cex.main = 0.7)
# bxplot(SelfRenawlDuration, add = T, col = "red")
#############################################



boxplot(SelfRenawlDuration ~ maxRRank,  pch = "+", ylab = "Self-renewal duration from 1st division [Days]", xlab = "Number of successive R cell division in the lineage", cex.lab = 0.7)
beeswarm(SelfRenawlDuration ~ maxRRank, data = NeuronDeathTable, 
         add = T, pch = 20, col = add.alpha('red',alpha = 0.6), cex =1, corral = "wrap")

typeDuration <- c(rep("RSelfrenewal",length(SelfRenawlDuration)), rep("activityClone", length(activityDurationCloneNonZero)))

boxplot(c(SelfRenawlDuration, activityDurationCloneNonZero) ~ typeDuration, ylab = "Duration [Days]")
beeswarm(c(SelfRenawlDuration, activityDurationCloneNonZero) ~ typeDuration, add = T,pch = 20, col = add.alpha('red',alpha = 0.6), corral = "wrap")

RenewalActivity <- data.frame(SelfRenawlDuration, activityDuration = activityDurationClone[match(names(maxRRank), names(activityDurationClone))])
# Can be saved for GraphPad analysis

plot(RenewalActivity, xlab = "Self-renewing duration of R cell in each clone [days]", ylab = "activity duration of the clone [days]")
abline(a = 0, b = 1, col = "red", lwd = 2)



# dataMaxR <- dataMax[dataMax$CellType2 == 1 & dataMax$UncertaintyType2 == 1,]

# lineageMin2Rdiv <- with(dataMaxR, tapply(DPI, Clone_uniqueID, length))>1 # min 2 successive divisions of R cell
# lineageMin2RdivID <- names(lineageMin2Rdiv)[lineageMin2Rdiv]
# dataMaxRMin2RDiv <- dataMaxR[dataMaxR$Clone_uniqueID %in% lineageMin2RdivID,]
# dotplot(as.factor(dataMaxRMin2RDiv$Clone_uniqueID) ~ dataMaxRMin2RDiv$DPI, xlab = "Time (Days post induction)", type = "b")



# At least one division of R cell but remove quiescent clones
# Get the ID of the lineage with a sure rank 1 R cell that divides
IdRsureDivides <- dataMaxR$Clone_uniqueID[dataMaxR$CellType2 == 1 & dataMaxR$UncertaintyType2 == 1 & dataMaxR$Rank == 1]
dataMaxRMin2RDiv <- dataMaxR[dataMaxR$Clone_uniqueID %in% IdRsureDivides,]

dataMaxRMin2RDiv$Rank <- as.factor(dataMaxRMin2RDiv$Rank)

## Gets the right order of the lineages for the plot
dataMaxRMin2RDivRrank1 <- dataMaxRMin2RDiv[dataMaxRMin2RDiv$Rank == 1, ]
dataMaxRMin2RDivRrank1$SelfrenewalDuration <- SelfRenawlDuration[match(as.character(dataMaxRMin2RDivRrank1$Clone_uniqueID),names(maxRRank))]
dataMaxRMin2RDivRrank1$Clone_uniqueID <- factor(dataMaxRMin2RDivRrank1$Clone_uniqueID, levels(as.factor(dataMaxRMin2RDivRrank1$Clone_uniqueID))[ match( dataMaxRMin2RDivRrank1$Clone_uniqueID,levels(as.factor(dataMaxRMin2RDivRrank1$Clone_uniqueID)))])
newLevels <- match(dataMaxRMin2RDivRrank1$DPI, sort(dataMaxRMin2RDivRrank1$DPI, decreasing = F)) 

# Deals with ex aequo, here use the selfrenewal time
for(i in c(1:length(newLevels))){
  if(sum(newLevels == i) > 1){
    exequoID <- which(newLevels == i)
    
    #for (j in c(1:length(exequoID))){
    # newLevels[exequoID[j]] <- newLevels[exequoID[j]]+ (j-1)
    #}
    newOrder <- match(dataMaxRMin2RDivRrank1$SelfrenewalDuration[exequoID],sort(dataMaxRMin2RDivRrank1$SelfrenewalDuration[exequoID], decreasing = F))
    newLevels[exequoID] <- newLevels[exequoID] + (newOrder-1)
    if(sum(duplicated(newLevels[exequoID])) > 0){
      temp <- newLevels[exequoID]
      SelfRenawlExequoID <- which(duplicated(newLevels[exequoID]) == T)
      
      for (j in c(1:length(SelfRenawlExequoID))){
        newLevels[exequoID[SelfRenawlExequoID[j]]] <- temp[SelfRenawlExequoID[j]] + (1)
      }
    } #else {newLevels[exequoID] <- newLevels[exequoID] + (newOrder-1)}
    
  }
}

dataMaxRMin2RDiv$Clone_uniqueID <- factor(x = dataMaxRMin2RDiv$Clone_uniqueID, levels = levels(as.factor(dataMaxRMin2RDivRrank1$Clone_uniqueID))[match(c(1:length(IdRsureDivides)), newLevels)])

ggplot(dataMaxRMin2RDiv, aes(Clone_uniqueID, DPI)) +
  geom_point() +
  coord_flip() +
  geom_line(aes(group = Clone_uniqueID)) +
  geom_point(aes(color = Rank)) + 
  ggtitle("Self-renewal duration of R cells") +
  labs(y = "Duration [Days]", x = "Lineage ID")


################################################################################################################################
# Types of cells, rank of division
################################################################################################################################

#############################################
# Proportions of R and NR progenitors in the whole dataset. At which generation the R and NR progenitors appear in the recorded lineage trees. 
# The rank of division is the generation of the lineage tree the division occured. It also gives an idea how many successive divisions (how long) are in the lineage trees.

# Only 3 astocytes, and 1 astrocyte division, leaving it out for the Tc

# For the NR cells rank 1 corresponds to the first apparition of an RN cell generated by a R cell in the branch.
#############################################

# Sara
#############################################
# divisionTableMotherSure <- divisionTable[divisionTable$MotherUncertaintyType == 1,]
#############################################

# Yicheng: add MotherUncertaintyType == "NA"
#############################################
# filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyType = c("1","0.5","0", "NA"), uncertaintyLineage = c("1"))
#############################################

filteredDivisionTable <- filteredDivisionTable[filteredDivisionTable$MotherCellType != "5",] 
tab <- table(filteredDivisionTable$MotherCellType)
barplot(prop.table(tab), col = c("red", "orange"), ylab = "frequency of mother cell type", names.arg = c("R", "NR"), 
        main = paste("Proportions of progenitor types \nin the whole dataset, n= ", sum(tab), " divisions"), cex.main = 0.8, cex.axis = 0.7)

tab <- table(filteredDivisionTable$DivisionRank)
barplot(prop.table(tab), col = c("white"), ylab = "frequency of division", main = paste("Proportions of division ranks \nin the whole dataset, n= ", 
        sum(tab), " divisions"), cex.main = 0.8, cex.axis = 0.7, xlab = "rank of division")

tab <- table(filteredDivisionTable$MotherCellType, filteredDivisionTable$DivisionRank)
barplot(prop.table(tab,margin = 2), col = c("red", "orange"), ylab = "proportions of progenitors", xlab = "rank of division")

maxRank <- max(as.numeric(as.character(filteredDivisionTable$DivisionRank)))
tabRank <- table(filteredDivisionTable$DivisionRankNR, filteredDivisionTable$MotherCellType)
barplot(tabRank, beside = TRUE, names.arg = c("R","NR special rank"), col = rainbow(maxRank-1), ylab = "Number of divisions")
legend(x= 4, y= 60, legend = c(1:maxRank), fill = rainbow(maxRank-1), ncol = 2, title = "rank of division", cex = 0.7, bty = 'n')

#############################################
# Number of NR successive divisions
#############################################

dataMaxNR <- dataMax[dataMax$CellType2 == 2 & dataMax$UncertaintyType2 ==1,]
maxNRRank <- with(dataMaxNR, tapply(RankNR, as.factor(Clone_uniqueID), max))  
tab <- table(maxNRRank)
barplot(tab, xlab = "Number of NR successive divisions", ylab = "Number of lineages", col = "white", main = paste("mean = ", 
        round(mean(maxNRRank),2), "sem = ", round(sd(maxNRRank)/sqrt(length(maxNRRank)),2), " n = ", length(maxNRRank)))


beeswarm(maxNRRank, ylab = "Number of NR successive divisions", xlab = "Number of lineages", main = paste("mean = ", round(mean(maxNRRank),2), "sem = ", 
         round(sd(maxNRRank)/sqrt(length(maxNRRank)),2), " n = ", length(maxNRRank)))
bxplot(maxNRRank, col = "red", add = T)

# Don't understand
# barPlotBeeswarm(maxNRRank,  ylab = "Number of NR successive divisions", xlab = "Number of lineages", cex.main  = 0.8)
# text(0.75, 6, "one lineage no NR cells")



################################################################################################################################
# Cell cycle length
################################################################################################################################

#############################################
## So far the cell cycle lenght (Tc) is measured only for cells for which transitions from mother and to daughters are certain. 
## Which cells to take into account is still to fully determine. For the details of Tc for R or NR progenitors, only cells for which the type is certain are taken.

## Each division happens during a time window, the last time point the mother is seen and the first time point the daughters are seen, this timing uncertainty can be in some cases long. 
## For every potential Tc, the timing uncertainties of the 2 successives divisions are added and the Tc is discarted if this period is greater than 6 DPI.

## So far a minimum Tc is calculated which is the time (in days) between cell birth and cell lastly seen before division, and a maximum Tc time between mother last seen and daughters first seen. 
## The mean between these 2 values is taken as the Tc.

## Data is displayed as mean +/- sem.
#############################################


divisionTable <- transformTreeIntoDivision(decodedData, cellType = "simplified")   

filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyType = c("1","0.5","0"), uncertaintyLineage = c("1"))
filteredDivisionTable <- filteredDivisionTable[filteredDivisionTable$MotherCellType != "5",]  # Remove the astrocyte mother cell

tabTc <- getTc(filteredDivisionTable) # Every row corresponds to one daughter that is also a mother = whose Tc can be measured

# Yicheng: add MotherUncertaintyType == "NA"
#############################################
# filteredDivisionTable  <- filterDivisionTable(divisionTable, uncertaintyType = c("1","0.5","0", "NA"), uncertaintyLineage = c("1"))
# filteredDivisionTable <- filteredDivisionTable[filteredDivisionTable$MotherCellType != "5",]

# tabTc <- getTc(filteredDivisionTable)
#############################################

hist(tabTc$totalWindowTc, main = paste("Distributions of time uncertainty, n = ",length(tabTc$meanTc)), ylab = "# of Tc", xlab = "Uncertainty [Days]", cex.main = 0.9, breaks = seq(from=0, to= 20, by= 2))
thresholdDivisionWindow <- 6 # Maximum time uncertainty in the 2 divisions taken together
abline(v = thresholdDivisionWindow, col = "red", lwd = 2)
tabTc <- tabTc[tabTc$totalWindowTc < thresholdDivisionWindow, ]
cellIDTc <- tabTc$MotherID

hist(tabTc$meanTc, main = paste("Distributions of Tc, n = ",length(tabTc$meanTc),"\n mean = ", round(mean(tabTc$meanTc),1),"days"), ylab = "# of Tc", xlab = "Tc [Days]", cex.main = 0.9)
abline(v = mean(tabTc$meanTc), col = "red", lwd = 2)
# plot line which corresponds to th mean


kt <- kruskal.test(tabTc$meanTc ~ tabTc$DivisionRank)
tc <- plotMeanTc(tabTc = tabTc, cond = list(tabTc$DivisionRank), color = rainbow(7), my.main = paste("Tc for every rank of division, \nall progenitors, pvalue =", round(kt$p.value,4),"\nKruskal-Wallis test" ))


################################################################################################################################
# Sister Tc
################################################################################################################################

#############################################
## Tc of sister cells
#############################################

MotherSisterTc <- tabTc$MotherID[duplicated(tabTc$MotherID)] # Mother of 2 sisters with Tc --> 2 occurence in motherID

tabTcSisters <- tabTc[tabTc$MotherID %in% MotherSisterTc,]    # Select the table of sister Tc

TcSisters <- t(matrix(tabTcSisters$meanTc, ncol = length(tabTcSisters$meanTc)/2))


# Puts the shorter Tc in the first column
for (i in c(1:nrow(TcSisters ))){
  TcSisters [i,1:2] <- TcSisters [i,order(TcSisters [i,1:2], decreasing = F)]
}

myCor <- cor.test(TcSisters[,1], TcSisters[,2], method = "spearman", exact = F)
myLm <- lm(TcSisters[,2] ~ TcSisters[,1])
plot(jitter(TcSisters[,1],factor = 1), jitter(TcSisters[,2], factor = 1), pch = 3, main = paste("Correlation Tc of sister cells \n n = ", 
            nrow(TcSisters), " pairs, Spearman correlation rho = ", round(myCor$estimate,4),"\n p = ", round(myCor$p.value,9)), ylab = "longer sister Tc [days]", xlab = "shorter sister Tc [days]", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)

# log log plots and correlations:
myCor <- cor.test(log(TcSisters[,1]), log(TcSisters[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSisters[,2]) ~ log(TcSisters[,1]))
plot(jitter(log(TcSisters[,1]),factor = 1), jitter(log(TcSisters[,2]), factor = 1), pch = 3, main = paste("Correlation log(Tc) of sister cells \n n = ", 
            nrow(TcSisters), " pairs, Spearman correlation rho = ", round(myCor$estimate,3),"\n p = ", round(myCor$p.value,9)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)


#############################################
## Correlation test with simulation
## Take randomly pairs of Tc within the measured Tc vector and check whether we get a correlation
#############################################

# Method1 takes randomly 2 observed Tc (not sisters) and check for a potential correlation

#####
# Kind of weird, with ceiling the random numbers should never be 0 but it does happen often... cannot find out why


npair <- dim(TcSisters)[1]  # Number of pairs to be taken at every simulation (same number as the observed sister Tc)
iteration <- 1000            # Repeats the simulation 100 times  

# Prealocates variables
randomRatio <- vector("numeric", length= iteration * npair)
randomRho <- vector("numeric", length = iteration)
randomPvalue <- vector("numeric", length = iteration)

for(iter in c(1:iteration)){        # Loops in the simulations 
  randomPair <- matrix(nrow = npair, ncol = 2) 
  for(i in c(1:npair)){             # Loops in the pairs
    pair <- ceiling(runif(2)*length(tabTc$meanTc))  # Randomly draws 2 ID of observed Tc
    
    while(pair[1] == pair[2]){ # to make sure we are not selecting sister Tc (#abs(pair[2]-pair[1]) == 1)
      pair <- round(runif(2)*length(tabTc$meanTc))
    }
    
    if(sum(pair == 0) > 0){ # If the ID is 0
      #print(paste(pair[1],pair[2], "in ID \n"))
      pair[pair == 0] <- 1
      #print(paste(pair))
    }
    
    randomPair[i, ] <- tabTc$meanTc[pair]
  }
  # Orders the pairs first column = shorter Tc
  for (i in c(1:nrow(randomPair))){
    randomPair[i,1:2] <- randomPair[i,order(randomPair[i,1:2], decreasing = F)]
  }
  
  randomRatio[((i-1)*npair + 1) : (i*npair)] <- randomPair[,2]/randomPair[,1]
  
  randomCor <- cor.test(randomPair[,1], randomPair[,2], method = "spearman", exact = F)
  randomRho[iter] <- randomCor$estimate
  randomPvalue[iter] <- randomCor$p.value
  
  # Vizualization of the simulation results
  #randomLm <- lm(randomPair[,2] ~ randomPair[,1])
  #plot(randomPair[,1], randomPair[,2], pch = 3, main = paste("Correlation Tc of random pair cells \n n = ", nrow(randomPair), " pairs, Spearman Correlation rho = ", round(randomCor$estimate,4),", p = ", round(randomCor$p.value,5)), ylab = "longer Tc", xlab = "shorter Tc", cex.main = 0.9)
  #abline(randomLm, col = "red", lty = 2)
}

hist(randomRho, main = "Histogram of random Tc pairs (100 times 43 pairs) \n Observed Sister Tc correlation is outside the 95% condidence interval of \n correlation between random Tc pairs \n the sister Tc are statistically more correlated than random Tc pairs", 
     breaks = seq(-0.2,0.9,0.05), xlab = "Spearman Correlation Rho for unrelated pair Tc", cex.main = 0.8)
abline(v = myCor$estimate, col = "red", lwd = 2)
abline(v = quantile(randomRho, 0.95), col = "green", lwd = 2)
# abline(v = c(mean(randomRho)+sd(randomRho)/sqrt(iteration), mean(randomRho)-sd(randomRho)/sqrt(iteration)), lwd = 2, col = "green")
text(x = 0.8, y= 0.1*iteration, labels = "95% data before this point", cex = 0.8, col = "green")
text(x = 0.8, y= 0.15*iteration, labels = "Observed sister Tc \ncorrelation", cex = 0.8, col = "red")


#############################################
## Method with fit
#############################################

# # fit the distribution of the Tc
# library(fitdistrplus)
# 
# fitTcExp <- fitdist(tabTc$meanTc,distr = "exp", method = "mme")
# lamda <- fitTcExp$estimate
# 
# # Check the validity of the exponential fit
# # plot(fitTcExp)
# # simData <- rexp(40000, lamda)
# # simData <- matrix(simData, nrow = 1000, ncol = 40)
# # meanSimData <- apply(simData,1,FUN = mean); meanSimData
# # sd(meanSimData)
# # hist(meanSimData)
# # se <- sd(meanSimData)/40
# # lower <- mean(meanSimData)- se
# # upper <- mean(meanSimData)+ se
# 
# 
# 
# #totalSisterTc <- c(TcSisters[,1], TcSisters[,2])
# randomRho <- vector("numeric", length = 100)
# randomPvalue <- vector("numeric", length = 100)
# 
# for (iteration in c(1:100)){
#   simTc <- rexp(92, lamda)
#   randomPairs <- NULL
#   
#   for (i in seq(1,length(simTc)/2,1)){
#     pair <- runif(2)*length(simTc)
#     #print(round(pair))
#     randomPairs <- rbind(randomPairs, simTc[round(pair)])
#   }
#   
#   # orders the pairs first column = shorter Tc
#   for (i in c(1:nrow(randomPairs))){
#     randomPairs[i,1:2] <- randomPairs[i,order(randomPairs[i,1:2], decreasing = F)]
#   }
#   randomCor <- cor.test(randomPairs[,1], randomPairs[,2], method = "spearman", exact = F)
#   randomRho[iteration] <- randomCor$estimate
#   randomPvalue[iteration] <- randomCor$p.value
# }
# 
# #randomLm <- lm(randomPairs[,2] ~ randomPairs[,1])
# #plot(randomPairs[,1], randomPairs[,2], pch = 3, main = paste("Correlation Tc of random pair cells \n n = ", nrow(randomPairs), " pairs, Spearman Correlation rho = ", round(randomCor$estimate,4),", p = ", round(randomCor$p.value,10)), ylab = "longer Tc", xlab = "shorter Tc")
# #abline(randomLm, col = "red", lty = 2)
# 
# # log log plots and correlations:
# #randomCor <- cor.test(log(randomPairs[,1]), log(randomPairs[,2]), method = "spearman", exact = F)
# #randomLm <- lm(log(randomPairs[,2]) ~ log(randomPairs[,1]))
# # plot(jitter(log(randomPairs[,1]), factor = 1), jitter(log(randomPairs[,2]), factor=1), pch = 3, 
# #      main = paste("Correlation log(Tc) of random pair cells \n n = ", nrow(randomPairs), " pairs, Spearman Correlation rho = ", 
# #      round(randomCor$estimate,4),", p = ", round(randomCor$p.value,10)), ylab = "longer Tc", xlab = "shorter Tc")
# #abline(randomLm, col = "red", lty = 2)
# 
# hist(randomRho, main = "Correlation between random pairs")
# hist(randomPvalue, main = "Pvalue Correlation between random pairs")

#############################################
## Sister Tc vs mother type
#############################################

tabTcSistersMotherTypeSure <- tabTcSisters[tabTcSisters$MotherUncertaintyType == 1,]
motherSisterType <- tabTcSistersMotherTypeSure$MotherCellType[seq(1, dim(tabTcSistersMotherTypeSure)[1],2)]

TcSisters <- t(matrix(tabTcSistersMotherTypeSure$meanTc, ncol = length(tabTcSistersMotherTypeSure$meanTc)/2))

# Need to make sure that shorter Tc in the first column
for (i in c(1:nrow(TcSisters ))){
  TcSisters [i,1:2] <- TcSisters [i,order(TcSisters [i,1:2], decreasing = F)]
}
myCol <- c("red", "orange")

myCor <- cor.test(TcSisters[,1], TcSisters[,2], method = "spearman", exact = F)
myLm <- lm(log(TcSisters[,2]) ~ log(TcSisters[,1]))

plot(jitter(log(TcSisters[,1]),factor = 1), jitter(log(TcSisters[,2]), factor = 1), pch = 3, col = myCol[motherSisterType], 
     main = paste("Correlation Tc of sister cells / mother type \n n = ", nrow(TcSisters), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = c("R mother","NR mother"), fill = myCol, cex = 0.8)


#############################################
## Only R mothers
#############################################

TcSistersRmother <- TcSisters[motherSisterType == 1,]

# Need to make sure that all the cells that die (column 1 shorter survival duration) or survive (column 2) are in the same column. 
for (i in c(1:nrow(TcSistersRmother))){
  TcSistersRmother[i,1:2] <- TcSistersRmother[i,order(TcSistersRmother[i,1:2], decreasing = F)]
}

myCorRmother <- cor.test(log(TcSistersRmother[,1]), log(TcSistersRmother[,2]), method = "spearman", exact = F)
myLmRmother <- lm(log(TcSistersRmother[,2]) ~ log(TcSistersRmother[,1]))


#############################################
## Only NR mothers
#############################################

TcSistersNRmother <- TcSisters[motherSisterType == 2,]

# Need to make sure that all the cells that die (column 1 shorter survival duration) or survive (column 2) are in the same column. 
for (i in c(1:nrow(TcSistersNRmother))){
  TcSistersNRmother[i,1:2] <- TcSistersNRmother[i,order(TcSistersNRmother[i,1:2], decreasing = F)]
}

myCorNRmother <- cor.test(log(TcSistersNRmother[,1]), log(TcSistersNRmother[,2]), method = "spearman", exact = F)
myLmNRmother <- lm(log(TcSistersNRmother[,2]) ~ log(TcSistersNRmother[,1]))


plot(jitter(log(TcSisters[,1]),factor = 1), jitter(log(TcSisters[,2]), factor = 1), pch = 3, col = myCol[motherSisterType], 
     main = paste("Correlation Tc of sister cells / mother type \n n = ", nrow(TcSisters), " pairs, Spearman correlation for \nR daughters rho = ", 
     round(myCorRmother$estimate,3),"p = ", round(myCorRmother$p.value,5), "\nfor NR daughters rho = ", round(myCorNRmother$estimate,3),"p = ", round(myCorNRmother$p.value,5)), 
     ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.7)
abline(myLmRmother, col = "red", lty = 2)
abline(myLmNRmother, col = "orange", lty = 2)

legend("bottomright", legend = c("R mother","NR mother"), fill = myCol, cex = 0.8)

# Figure 3 CD
# pdf("Figure3CD_sameAxes.pdf", width = 8, height = 4)
## par(mfrow = c(1,2))
plot(jitter(log(TcSistersRmother[,1]),factor = 1), jitter(log(TcSistersRmother[,2]), factor = 1), pch = 3, 
     main = paste("Correlation log(Tc) of sister cells R mothers \n n = ", nrow(TcSistersRmother), " pairs, Spearman correlation rho = ", 
     round(myCorRmother$estimate,3),"\n p = ", round(myCorRmother$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8, xlim = c(-0.5, 2.8), ylim = c(-0.5, 2.8))
abline(myLmRmother, col = "red", lty = 2)

plot(jitter(log(TcSistersNRmother[,1]),factor = 1), jitter(log(TcSistersNRmother[,2]), factor = 1), pch = 3, 
     main = paste("Correlation log(Tc) of sister cells NR mothers \n n = ", nrow(TcSistersNRmother), " pairs, Spearman correlation rho = ", 
     round(myCorNRmother$estimate,3),"\n p = ", round(myCorNRmother$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8, xlim = c(-0.5, 2.8), ylim = c(-0.5, 2.8))
abline(myLmNRmother, col = "red", lty = 2)

# dev.off()

###########################
# Sister Tc vs sister types
##########################

tabTcSistersSisterTypeSure <- tabTcSisters[tabTcSisters$DaughterUncertaintyType == 1,]

# Take only pairs where both sisters are sure, be careful if 2 pairs have only 1 cell sure, then no error!!
tabTcSistersSisterTypeSure <- tabTcSistersSisterTypeSure[tabTcSistersSisterTypeSure$MotherID %in% tabTcSistersSisterTypeSure$MotherID[duplicated(tabTcSistersSisterTypeSure$MotherID)],]
sisterTypes <- colSums(matrix(tabTcSistersSisterTypeSure$DaughterCellType,nrow = 2, ncol = length(tabTcSistersSisterTypeSure$DaughterID)/2))-1


TcSisters <- t(matrix(tabTcSistersSisterTypeSure$meanTc, ncol = length(tabTcSistersSisterTypeSure$meanTc)/2))
for (i in c(1:nrow(TcSisters ))){
  TcSisters [i,1:2] <- TcSisters [i,order(TcSisters [i,1:2], decreasing = F)]
}
myCol <- c("red", "blue","orange")
myLegend <- c("R-R sisters","R-NR sisters", "NR-NR sisters")
myPch <- c(0,1,2)

# TcSisters <- TcSisters[sisterTypes ==3,]

myCor <- cor.test(log(TcSisters[,1]), log(TcSisters[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSisters[,2]) ~ log(TcSisters[,1]))


plot(jitter(log(TcSisters[,1]),factor = 2), jitter(log(TcSisters[,2]), factor = 2), col = myCol[sisterTypes], pch = myPch[sisterTypes], 
     main = paste("Correlation Tc of sister cells /sister types \n n = ", nrow(TcSisters), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = myLegend, fill = myCol)

##########################
# NR-NR pairs
##########################

TcNRsisters <- TcSisters[sisterTypes == 3,]

myCor <- cor.test(log(TcNRsisters[,1]), log(TcNRsisters[,2]), method = "spearman", exact = F)
plot(jitter(log(TcNRsisters[,1]),factor = 2), jitter(log(TcNRsisters[,2]), factor = 2), pch = 3, 
     main = paste("Correlation Tc of sister cells / NR sisters \n n = ", nrow(TcNRsisters), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)

##########################
# R-NR pairs
##########################

TcR.NRsisters <- TcSisters[sisterTypes != 3,]
myCor <- cor.test(log(TcR.NRsisters[,1]), log(TcR.NRsisters[,2]), method = "spearman", exact = F)
plot(jitter(log(TcR.NRsisters[,1]),factor = 2), jitter(log(TcR.NRsisters[,2]), factor = 2), pch = 3, 
     main = paste("Correlation Tc of sister cells / R.R or R.NR sisters  \n n = ", nrow(TcR.NRsisters), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)

##########################
# Type of division with mother info
##########################

tabTcSistersMotherSisterTypeSure <- tabTcSisters[tabTcSisters$DaughterUncertaintyType == 1 & tabTcSisters$MotherUncertaintyType == 1,]

# Take only pairs where both sisters are sure, be careful if 2 pairs have only 1 cell sure, then no error!!
tabTcSistersMotherSisterTypeSure <- tabTcSistersMotherSisterTypeSure[tabTcSistersMotherSisterTypeSure$MotherID %in% tabTcSistersMotherSisterTypeSure$MotherID[duplicated(tabTcSistersMotherSisterTypeSure$MotherID)],]
sisterTypes <- colSums(matrix(tabTcSistersMotherSisterTypeSure$DaughterCellType,nrow = 2, ncol = length(tabTcSistersMotherSisterTypeSure$DaughterID)/2))-1
motherTypes <- tabTcSistersMotherSisterTypeSure$MotherCellType[seq(1, dim(tabTcSistersMotherSisterTypeSure)[1],2)]

TcSisters <- t(matrix(tabTcSistersMotherSisterTypeSure$meanTc, ncol = length(tabTcSistersMotherSisterTypeSure$meanTc)/2))
for (i in c(1:nrow(TcSisters ))){
  TcSisters [i,1:2] <- TcSisters [i,order(TcSisters [i,1:2], decreasing = F)]
}
myCol <- c("red", "orange")
myLegend <- c("R mother","NR mother")

##########################
##  NR-NR pairs
##########################

TcSistersNR <- TcSisters[sisterTypes ==3,]
motherTypesNR <- motherTypes[sisterTypes ==3]

myCor <- cor.test(log(TcSistersNR[,1]), log(TcSistersNR[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSistersNR[,2]) ~ log(TcSistersNR[,1]))

plot(jitter(log(TcSistersNR[,1]),factor = 2), jitter(log(TcSistersNR[,2]), factor = 2), col = myCol[motherTypesNR], pch = 3, 
     main = paste("Correlation Tc of sister cells /NR sister types \n n = ", nrow(TcSistersNR), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = c("R mother","NR mother"), fill = myCol, cex = 0.8)

##########################
##  R-R NR-R pairs
##########################

TcSistersR <- TcSisters[sisterTypes != 3,]
motherTypesR <- motherTypes[sisterTypes != 3]

myCor <- cor.test(log(TcSistersR[,1]), log(TcSistersR[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSistersR[,2]) ~ log(TcSistersR[,1]))

plot(jitter(log(TcSistersR[,1]),factor = 2), jitter(log(TcSistersR[,2]), factor = 2), col = myCol[motherTypesR], pch = 3, 
     main = paste("Correlation Tc of sister cells / R-R R-NR sister types \n n = ", nrow(TcSistersR), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = c("R mother","NR mother"), fill = myCol, cex = 0.8)

##########################
## R mothers, colors sister cell types
##########################

TcSistersMotherR <- TcSisters[motherTypes == 1,]
sisterTypesMotherR <- sisterTypes[motherTypes == 1]

myCor <- cor.test(log(TcSistersMotherR[,1]), log(TcSistersMotherR[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSistersMotherR[,2]) ~ log(TcSistersMotherR[,1]))

myCol <- c("red", "blue","orange")
myLegend <- c("R-R sisters","R-NR sisters", "NR-NR sisters")
myPch <- c(0,1,2)

plot(jitter(log(TcSistersMotherR[,1]),factor = 2), jitter(log(TcSistersMotherR[,2]), factor = 2), col = myCol[sisterTypesMotherR], pch = myPch[sisterTypesMotherR], 
     main = paste("Correlation Tc of sister cells /sister types for R mothers \n n = ", nrow(TcSistersMotherR), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = myLegend, fill = myCol, cex = 0.8)


##########################
## Tc according to cell types
##########################

tabTc <- tabTc[tabTc$DaughterUncertaintyType ==1, ] # filters for Tc whose cell type is certain
kt <- kruskal.test(tabTc$meanTc ~ tabTc$DaughterCellType)
tc <- plotMeanTc(tabTc = tabTc, cond = list(tabTc$DaughterCellType), myLegend = c("R","NR"), my.main = paste("Tc for progenitor types, \n pvalue =", round(kt$p.value,4),"Kruskal-Wallis test"))

# plot the histogramms of Tc for R and NR progenitors
listTc <- with(tabTc, tapply(X = meanTc, INDEX = as.factor(DaughterCellType), FUN = identity))

TcR <- unlist(listTc[[1]])
hist(TcR, main = paste("Histogram of R Tc \n mean =", round( mean(TcR),1), "days, n =", length(TcR)), xlab = "Tc [Days]", ylab = "# of Tc", cex.main = 0.9)
abline(v = mean(TcR), col = "red", lwd = 2)
barPlotBeeswarm(TcR, ylab = "Tc R cells")

TcNR <- unlist(listTc[[2]])
hist(TcNR, main = paste("Histogram of NR Tc \n mean =", round( mean(TcNR),1), "days, n =", length(TcNR)), xlab = "Tc [Days]", ylab = "# of Tc", cex.main = 0.9)
abline(v = mean(TcNR), col = "red", lwd = 2)
barPlotBeeswarm(TcNR, ylab = "Tc NR cells")


# my.settings <- list(plot.symbols = list(col=c("red", "orange")),
#                     box.rectangle=list(col=c("red", "orange")),
#                     box.dot = list(col=c("black")),
#                     box.umbrella = list(col=c("black")))
# 
# with(tabTc, bwplot(DaughterCellType ~ meanTc, notch = TRUE, par.settings = my.settings, fill = c("red", "orange")))

with(tabTc, boxplot(meanTc ~ DaughterCellType,  col = c("red", "orange"), varwidth = F, horizontal = T, notch = F, arg.names = c("R","NR"), xlab = "mean Tc [Days]", ylab = "Cell type"))

# beeswarn(meanTc ~ DaughterCellType,  data = tabTc, 
#     add = T, pch = 20, col = add.alpha('blue',alpha = 0.6), cex =1, corral = "wrap")

tabTcR <- tabTc[tabTc$DaughterCellType ==1, ] # filters for  R
tabTcR$DivisionRank <- droplevels(tabTcR$DivisionRank)
kt <- with(tabTcR, kruskal.test(meanTc ~ DivisionRank))
tcRRank <- plotMeanTc(tabTc = tabTcR, cond = list(tabTcR$DivisionRank), color = rainbow(7), my.main = paste("Tc for every rank of division \nR progenitors, p =", round(kt$p.value,4)))

tabTcNR <- tabTc[tabTc$DaughterCellType ==2, ] # filters for  NR
kt <- with(tabTcNR, kruskal.test(meanTc ~ DivisionRankNR))
tc <- plotMeanTc(tabTc = tabTcNR, cond = list(tabTcNR$DivisionRankNR), color = rainbow(7), my.main = paste("Tc for every rank of division \nNR progenitors, p =", round(kt$p.value,4)), ylim = c(0,10))


cloneMeanTc <- with(tabTc, tapply(meanTc, INDEX =  factor(Clone_uniqueID), FUN = mean))
cloneSDTc <- with(tabTc, tapply(meanTc, INDEX =  factor(Clone_uniqueID), FUN = sd))
numberSuccessiveDivisionCloneTc <-  numberSuccessiveDivision[match(names(cloneMeanTc), names(numberSuccessiveDivision))]

boxplot(cloneMeanTc ~ numberSuccessiveDivisionCloneTc, pch = 16, cex = 0.8, ylab = "mean Tc per clone [Days]", xlab = "number of successive divisions per clone", main = "mean Tc per clone")

hist(cloneSDTc, xlab = "SD Tc per clone [Days]", main = "Tc variations within a clone")

boxplot(cloneSDTc ~ numberSuccessiveDivisionCloneTc, pch = 16, cex = 0.8, ylab = "SD Tc per clone [Days]", xlab = "number of successive divisions per clone", main = "Variation Tc per clone")
beeswarm(cloneSDTc ~ numberSuccessiveDivisionCloneTc, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), cex =1, corral = "wrap")


##########################
## Tc per clone
##########################

my.settings <- list(strip.background=list(col="black"),
                    strip.border=list(col="black"))

dotplot(factor(tabTc$DaughterCellType)  ~ jitter(tabTc$meanTc,factor = 0.5) | factor(tabTc$Clone_uniqueID), cex =0.8, col = add.alpha("blue", 0.6),
        par.settings = my.settings,  par.strip.text=list(col="white", font=1), box.width = 1/3, xlab = "Mean Tc [Days]")


##########################
# Here division types are defined as follows:
# Si : symmetric proliferative (for NR cells it is really inheriting, for R cells it is only in the rare cases when R -> 2R but mainly there are R -> 2NR)
# Sd : symmetric differentiative (producing 2 neurons)
# Ap : producing 2 progenitors of different types (R -> R + NR)
# Ad : generating one neuron and one progenitor
##########################

# Death not taken into account, only fate of last time point before death or div
divisionTable2 <- transformTreeIntoDivision(decodedData, cellType = "simplified no death")   

filteredDivisionTable2  <- filterDivisionTable(divisionTable2, uncertaintyType = c("1"), uncertaintyLineage = c("1"))
filteredDivisionTable2 <- getDivisionType(filteredDivisionTable2) # adds the division type as last column



tabTc <- getTc(filteredDivisionTable2) # every row corresponds to one daughter that is also a mother = whose Tc can be measured 
tabTc$DivisionType <- filteredDivisionTable2$DivType[match(tabTc$DaughterID, filteredDivisionTable2$MotherID)]

kt <- kruskal.test(tabTc$meanTc ~ tabTc$DivisionType)
detailsTc <- plotMeanTc(tabTc, cond = list(tabTc$DivisionType), myLegend = FALSE, my.main = paste("Tc for division types all progenitor types \n p = ", round(kt$p.value,4)), color = "white")
boxplot(tabTc$meanTc ~ tabTc$DivisionType, main = "Tc for division types all progenitor types", ylab = "Tc [Days]")
beeswarm(tabTc$meanTc ~ tabTc$DivisionType,  add = T, pch = 20, col = add.alpha('red',alpha = 0.6), cex =1, corral = "wrap")


tabTcR2 <- tabTc[tabTc$DaughterCellType == 1,]
kt <- kruskal.test(tabTcR2$meanTc ~ tabTcR2$DivisionType)
detailsTc <- plotMeanTc(tabTcR2, cond = list(tabTcR2$DivisionType), myLegend = FALSE, my.main = paste("Tc for division types R cells \np = ", round(kt$p.value,4)), color = "white")
boxplot(tabTcR2$meanTc ~ tabTcR2$DivisionType, main = "Tc for division types R cells", ylab = "Tc [Days]")
beeswarm(tabTcR2$meanTc ~ tabTcR2$DivisionType,  add = T, pch = 20, col = add.alpha('red',alpha = 0.6), cex =1, corral = "wrap")


tabTcNR <- tabTc[tabTc$DaughterCellType == 2,]
kt <- kruskal.test(tabTcNR$meanTc ~ tabTcNR$DivisionType)
detailsTc <- plotMeanTc(tabTcNR, cond = list(tabTcNR$DivisionType), myLegend = FALSE, my.main = paste("Tc for division types NR cells \np = ", round(kt$p.value,4)), color = "white")
boxplot(tabTcNR$meanTc ~ tabTcNR$DivisionType, main = "Tc for division types NR cells", ylab = "Tc [Days]")
beeswarm(tabTcNR$meanTc ~ tabTcNR$DivisionType, add = T, pch = 20, col = add.alpha('red',alpha = 0.6), cex =1, corral = "wrap")

##########################
## Sister Tc division type
##########################

MotherSisterTc <-  tabTc$MotherID[duplicated(tabTc$MotherID)]
tabTcSisters <- tabTc[tabTc$MotherID %in% MotherSisterTc,]
motherTypes <- tabTcSisters$MotherCellType[seq(1, dim(tabTcSisters)[1],2)]

TcSisters <- t(matrix(tabTcSisters$meanTc, ncol = length(tabTcSisters$meanTc)/2))

sisterDivTypeMatrix <- matrix(tabTcSisters$DivisionType, nrow = 2, ncol = length(tabTcSisters$DaughterID)/2)
f <- function(x){paste(x, collapse = "-")}
sisterDivType <- apply(X = sisterDivTypeMatrix, MARGIN = 2, FUN = f)
sisterDivTypeSimple <- rep("Diff", length(sisterDivType))

for (i in seq_along(sisterDivType)) {
  if (sisterDivType[i] == "Ap-Ad") {
    sisterDivType[i] <- "Ad-Ap"
    sisterDivTypeSimple[i] <- "Diff"
  }
  if (sisterDivType[i] == "Sd-Ad") {
    sisterDivType[i] <- "Ad-Sd"
    sisterDivTypeSimple[i] <- "Diff"
  }
  if (sisterDivType[i] == "Si-Ad") {
    sisterDivType[i] <- "Ad-Si"
    sisterDivTypeSimple[i] <- "Diff"
  }
  if (sisterDivType[i] == "Sd-Ap") {
    sisterDivType[i] <- "Ap-Sd"
    sisterDivTypeSimple[i] <- "Diff"
  }
  if (sisterDivType[i] == "Si-Ap") {
    sisterDivType[i] <- "Ap-Si"
    sisterDivTypeSimple[i] <- "Diff"
  }
  if (sisterDivType[i] == "Si-Sd") {
    sisterDivType[i] <- "Sd-Si"
    sisterDivTypeSimple[i] <- "Diff"
  }
  if (sisterDivType[i] == "Si-Si" | sisterDivType[i] == "Sd-Sd" | sisterDivType[i] == "Ad-Ad"| sisterDivType[i] == "Ap-Ap") {sisterDivTypeSimple[i] <- "Same"}
}

for (i in c(1:nrow(TcSisters ))){
  TcSisters [i,1:2] <- TcSisters [i,order(TcSisters [i,1:2], decreasing = F)]
}


### Kruskal-Wallis on the ratio of the Tc
# ratioTc <- TcSisters[,1]/TcSisters[,2]
# test <- data.frame(ratioTc, motherTypes = as.factor(motherTypes), sisterDivTypeSimple, sisterDivType)
# test.aov <- with(test, aov(ratioTc ~ motherTypes * sisterDivTypeSimple))

myCol <- c("purple","orange", "cyan", "red")
mylegend <- c("Ad","Ap", "Sd","Si")
myPch <- c(0:10)

myCor <- cor.test(log(TcSisters[,1]), log(TcSisters[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSisters[,2]) ~ log(TcSisters[,1]))


plot(jitter(log(TcSisters[,1]),factor = 2), jitter(log(TcSisters[,2]), factor = 2), col = myCol[as.factor(sisterDivTypeSimple)], pch = myPch[as.factor(sisterDivType)], 
     main = paste("Correlation Tc of sister cells / div type \n n = ", nrow(TcSisters), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = c("Different division types", "Same division Type"), fill = myCol[c(1,2)], cex = 0.8, bty = 'n')


TcSisterTcSameDivType <- TcSisters[sisterDivTypeSimple == "Same",]
sisterDivTypeSame <- as.factor(sisterDivType[sisterDivTypeSimple == "Same"])
sisterDivTypeSame <- droplevels(sisterDivTypeSame)
myCor <- cor.test(log(TcSisterTcSameDivType[,1]), log(TcSisterTcSameDivType[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSisterTcSameDivType [,2]) ~ log(TcSisterTcSameDivType [,1]))
plot(jitter(log(TcSisterTcSameDivType[,1]),factor = 2), jitter(log(TcSisterTcSameDivType[,2]), factor = 2), pch = myPch[sisterDivTypeSame], 
     main = paste("Correlation Tc of sister cells / same div type \n n = ", nrow(TcSisterTcSameDivType), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = levels(sisterDivTypeSame), pch = myPch[1:3], cex = 0.8, bty = 'n')

TcSisterTcDiffDivType <- TcSisters[sisterDivTypeSimple == "Diff",]
sisterDivTypeDiff <- as.factor(sisterDivType[sisterDivTypeSimple == "Diff"])
sisterDivTypeDiff <- droplevels(sisterDivTypeDiff)
myCor <- cor.test(log(TcSisterTcDiffDivType[,1]), log(TcSisterTcDiffDivType[,2]), method = "spearman", exact = F)
myLm <- lm(log(TcSisterTcDiffDivType [,2]) ~ log(TcSisterTcDiffDivType [,1]))
plot(jitter(log(TcSisterTcDiffDivType[,1]),factor = 2), jitter(log(TcSisterTcDiffDivType[,2]), factor = 2), pch = myPch[sisterDivTypeDiff], 
     main = paste("Correlation Tc of sister cells / Different div types \n n = ", nrow(TcSisterTcDiffDivType), " pairs, Spearman correlation rho = ", 
     round(myCor$estimate,3),"\n p = ", round(myCor$p.value,5)), ylab = "longer sister log(Tc)", xlab = "shorter sister log(Tc)", cex.main = 0.8)
abline(myLm, col = "red", lty = 2)
legend("bottomright", legend = levels(sisterDivTypeDiff), pch = myPch[1:length(levels(sisterDivTypeDiff))], cex = 0.8, bty = 'n')


##########################
## Minimum cell cycle length for root or quiescent cells
## Minimum length of cell cycle, starting from induction time
##########################

# Select root cells and quiescent R cells (rank = 1, celltype2 = 1), remove the putative root cells (time point 0)
rank1Cells <- dataMinMax[dataMinMax$Rank == 1 & dataMinMax$DPI > 0 & dataMinMax$CellType2 == 1 ,] 
## MODIFIED SARA##
# select root cells and quiescent R cells (rank = 1, celltype2 = 1), KEEP the putative root cells (time point 0)
# rank1Cells <- dataMinMax[dataMinMax$Rank == 1 & dataMinMax$CellType2 == 1 ,] 


minTc <- with(rank1Cells, tapply(X = DPI, INDEX = Cell_uniqueID, FUN = max))

par(mfrow= c(1,2))
hist(minTc, main = paste("Histogram of min Tc root / quiescent cells \n mean =", round( mean(minTc),1), "days, n =", length(minTc)), 
     xlab = "minimum Tc [Days]", ylab = "# of Tc", cex.main = 0.9)


# Only quiescent cells
quiescentCells <- rank1Cells[rank1Cells$Clone_uniqueID %in% quiescentClonesID,]
minTcQuiescent <- with(quiescentCells, tapply(X = DPI, INDEX = Cell_uniqueID, FUN = max))


# Only root cells
rootCells <- rank1Cells[!rank1Cells$Clone_uniqueID %in% quiescentClonesID,]
minTcRoot <- with(rootCells, tapply(X = DPI, INDEX = Cell_uniqueID, FUN = max))


# hist(minTcQuiescent, main = paste("Histogram of minimum Tc quiescent cells \n mean =", round( mean(minTcQuiescent),1), "days, n =", 
#      length(minTcQuiescent),"\n root cells \n mean =", round( mean(minTcRoot),1), "days, n =", length(minTcRoot)), 
#      xlab = "minimum Tc [Days]", ylab = "# of Tc", cex.main = 0.9, ylim = c(0,20), col=rgb(0,0,1,0.5), xlim = c(0,70), breaks = seq(from=0, to=70, by=5))
# hist(minTcRoot, main = paste("Histogram of min Tc root cells \n mean =", round( mean(minTcRoot),1), "days, n =", length(minTcRoot)), 
#     xlab = "minimum Tc [Days]", ylab = "# of Tc", cex.main = 0.9, add = T, ylim = c(0,20), col=rgb(1,0,0,0.5), xlim = c(0,70), breaks = seq(from=0, to=70, by=5))
# legend("topright", legend = c("Quiescent cells", "Root cells"), col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), pt.cex=2, pch=15)



# Only R root cells
# rootCellsR <- rank1Cells[!rank1Cells$Clone_uniqueID %in% quiescentClonesID & rank1Cells$CellType2 == 1 & rank1Cells$UncertaintyType2 ==1,]
# minTcRootR <- with(rootCellsR, tapply(X = DPI, INDEX = Cell_uniqueID, FUN = max))

## modified sara##
# It shows also the uncertain once (the putative)
rootCellsR <- rank1Cells[!rank1Cells$Clone_uniqueID %in% quiescentClonesID & rank1Cells$CellType2 == 1 ,] 
minTcRootR <- with(rootCellsR, tapply(X = DPI, INDEX = Cell_uniqueID, FUN = max))

hist(minTcRootR, xlab = "minimum Tc for R root cells", ylab = "Number of root cells", main = paste("Histogramm of minimum Tc root cell \n mean = ", 
    round(mean(minTcRootR),1), " sem = ", round(sd(minTcRootR)/sqrt(length(minTcRootR)),1), ", n = ", length(minTcRootR)))
abline(v = mean(minTcRootR), col = "red", lwd = 2)

# Statistical test
tcRRank1 <- tabTcR$meanTc[tabTcR$DivisionRank == 1]
tcRRank2 <- tabTcR$meanTc[tabTcR$DivisionRank == 2]

tabRoot <- data.frame(Tc =minTcRootR, Rank = "Root")
tabTcRRank1 <- data.frame(Tc = tcRRank1, Rank = "1")
tabTcRRank2 <- data.frame(Tc = tcRRank2, Rank = "2")

tcRRankRootvs1 <- rbind(tabRoot, tabTcRRank1)
tcRRank1vs2 <- rbind(tabTcRRank1, tabTcRRank2)

kruskalRootvsRank1 <- with(tcRRankRootvs1 , kruskal.test(Tc ~ Rank))
kruskalRank1vsRank2 <-with(tcRRank1vs2, kruskal.test(Tc ~ Rank))


# Plot of R cell Tc with minimum Tc for root cell
tcRRank[[1]] <- c(mean(minTcRootR),tcRRank[[1]]) # mean
tcRRank[[2]] <- c(sd(minTcRootR)/sqrt(length(minTcRootR)),tcRRank[[2]]) # sem
tcRRank[[3]] <- c(length(minTcRootR), tcRRank[[3]]) # n

names(tcRRank[[1]])[1] <- "root_Cell"
names(tcRRank[[2]])[1] <- "root_Cell"
myplot <- barplot(tcRRank[[1]], xlab = "Rank of division", ylab = "Tc [Days]", main = paste("Tc of R cells according to the Rank of Division \n root cell vs Rank 1 pvalue = ", 
                  round(kruskalRootvsRank1$p.value,5),"\n Rank 1 vs Rank 2 pvalue =", round(kruskalRank1vsRank2$p.value,3)), ylim = c(0,15), cex.main = 0.8, axis.lty = 1, col = c("gray","red","orange","yellow")) #ylim = c(0,max(tcRRank[[1]])+max(tcRRank[[2]]))
error.bar(myplot, y = tcRRank[[1]], upper = tcRRank[[2]])
text(x = myplot, y = 0.5, labels = tcRRank[[3]], cex = 0.8)





