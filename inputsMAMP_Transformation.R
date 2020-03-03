library(readr)
threatDistribution <- read_csv("data/data_Big/z_anonymous_threat_distribution.csv")
threatDistribution <- threatDistribution[,1:5]

threatDistribution_t <- t(threatDistribution)
threatDistribution_t <- threatDistribution_t[2:5,] #Rows = Threats / Columns = Planning units!
  
threatsDistribution_Tuple_c1 = c()
threatsDistribution_Tuple_c2 = c()
threatsDistribution_Tuple_c3 = c()

numberUnits   = ncol(threatDistribution_t)
numberThreats = nrow(threatDistribution_t) 
it = 0

for(i in 1:numberUnits){
  for(j in 1:numberThreats){
    if(threatDistribution_t[j,i] == 1){
      it    =it +1
      value = threatDistribution_t[j,i]
      print(it)
      threatsDistribution_Tuple_c1 = c(threatsDistribution_Tuple_c1, i)
      threatsDistribution_Tuple_c2 = c(threatsDistribution_Tuple_c2, j)
      threatsDistribution_Tuple_c3 = c(threatsDistribution_Tuple_c3, value)
    }
  }
}

threatsDistribution_Tuple = cbind(threatsDistribution_Tuple_c1, threatsDistribution_Tuple_c2, threatsDistribution_Tuple_c3 )
exportThreatsDistribution = matrix(data = threatsDistribution_Tuple, nrow = nrow(threatsDistribution_Tuple), ncol = ncol(threatsDistribution_Tuple))
dimnames(exportThreatsDistribution) <- list( c(1:nrow(exportThreatsDistribution) ), c("pu", "threats", "amount"))

head(exportThreatsDistribution)

#write.csv(exportThreatsDistribution, file="data05.csv", row.names = F)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
Amenazas_y_sensibilidades <- read_csv("data/data_Big/Amenazas y sensibilidades.csv")
numberSpecies = nrow(Amenazas_y_sensibilidades)
Amenazas_y_sensibilidades[,1] <- c(1:numberSpecies)
head(Amenazas_y_sensibilidades)

sensibilityData <- as.matrix(x = Amenazas_y_sensibilidades[,2:5])
str(sensibilityData)

data_x = c()
data_y = c()
data_z = c()

for(i in 1:numberSpecies){
  for(j in 1:numberThreats){
    if(sensibilityData[i,j] == 1){
      value = sensibilityData[i,j] 
      data_x = c(data_x, i)
      data_y = c(data_y, j)
      data_z = c(data_z, value)
    }
  }
}

sensibilityData_Tuple <- cbind(data_x, data_y, data_z)
head(sensibilityData_Tuple)
exportSensibility <- matrix(data = sensibilityData_Tuple, nrow = nrow(sensibilityData_Tuple), ncol = ncol(sensibilityData_Tuple))
dimnames(exportSensibility) <- list( c(1:nrow(exportSensibility) ), c("species", "threats", "sensibility"))
head(exportSensibility)
str(exportSensibility)

# write.csv(exportSensibility, file="data06.csv", row.names = F)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
actionCost <- read_csv("data/data_Big/data07.csv")
actionCost <- as.matrix(actionCost)
head(actionCost)
str(actionCost)
dimnames(actionCost) <- list( c(1:nrow(actionCost) ), c(1:ncol(actionCost) ))

numberUnits   = ncol(actionCost)
numberThreats = nrow(actionCost)

a_1 = c()
a_2 = c()
a_3 = c()

for(i in 1:numberUnits){
  for(j in 1:numberThreats){
    if(actionCost[j,i] != 0){
      value = actionCost[j,i]
      a_1 = c(a_1, i)
      a_2 = c(a_2, j)
      a_3 = c(a_3, value)
    }
  }
}

actionCost_Tuple <- cbind(a_1, a_2, a_3)
exportActionCost <- matrix(data = actionCost_Tuple, nrow = nrow(actionCost_Tuple), ncol = ncol(actionCost_Tuple))
dimnames(exportActionCost) <- list( c(1:nrow(exportActionCost) ), c("pu", "threats", "cost"))
head(exportActionCost)

write.csv(exportActionCost, file="data07_tuple.csv", row.names = F)




