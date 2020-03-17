# Activación de paquetes
#------------------------------------------------------------------------------------------
sys.source(file = "Package_RLibraries.R") 


# Lectura de datos (inputs) de MARXAN - Importación de Clases/Functiones de C++ (Rcpp)
#------------------------------------------------------------------------------------------
setwd("C:/Users/Irlanda Ceballos/Documents/GitHub/Conservation_Priorization")
sourceCpp(file = "src/RcppClass_MARXANData.cpp")
sourceCpp(file = "src/RcppFunction_GlobalFunctions.cpp")

#Lectura de datos (inputs) de MARXAN (inputs: 10 units/ 5 species/ 3 threats)
dataTarget       <- as.data.frame(read_delim("data/data_ExtremelySmall/target_ExtremelySmall.csv", ";", trim_ws = TRUE)%>%select(1:5))
dataCost         <- as.data.frame(read_delim("data/data_ExtremelySmall/unitCost_ExtremelySmall.csv", ";", trim_ws = TRUE))
dataBoundary     <- as.data.frame(read_delim("data/data_ExtremelySmall/boundary_ExtremelySmall.csv", ";",
                                             trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
dataDistribution <- as.data.frame(read_delim("data/data_ExtremelySmall/speciesDistribution_ExtremelySmall.csv", ";", trim_ws = TRUE))
dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
colnames(dataDistribution) <- c("pu","species","amount")


# # Lectura de datos (inputs) de MARXAN (inputs: 240 units/ 10 species/ 4 threats)
# dataTarget       <- as.data.frame(read_delim("data/data_Small/target_Small.csv", ";", trim_ws = TRUE)%>%select(1:5))
# dataCost         <- as.data.frame(read_delim("data/data_Small/unitCost_Small.csv", ";", trim_ws = TRUE))
# dataBoundary     <- as.data.frame(read_delim("data/data_Small/boundary_Small.csv", ";",
#                                              trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
# dataDistribution <- as.data.frame(read_delim("data/data_Small/speciesDistribution_Small.csv", ";", trim_ws = TRUE))
# dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
# colnames(dataDistribution) <- c("pu","species","amount")
# 
# 
# #Lectura de datos (inputs) de MARXAN (inputs: 2316 units/ 45 species/ 4 threats)
# dataTarget       <- as.data.frame(read_delim("data/data_Big/target_Big.csv", ";", trim_ws = TRUE)%>%select(1:5))
# dataCost         <- as.data.frame(read_delim("data/data_Big/unitCost_Big.csv", ";", trim_ws = TRUE))
# dataBoundary     <- as.data.frame(read_delim("data/data_Big/boundary_Big.csv", ";",
#                                              trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
# dataDistribution <- as.data.frame(read_delim("data/data_Big/speciesDistribution_Big.csv", ";", trim_ws = TRUE))
# dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
# colnames(dataDistribution) <- c("pu","species","amount")


# MARXAN con Penalización, con Rsymphony Solver (de la API de R)
#------------------------------------------------------------------------------------------
y             <- new(MARXANData, dataTarget, dataCost, dataBoundary, dataDistribution)
numberUnits   <- y$getUnits()
numberSpecies <- y$getSpecies()
bFactor       <- y$getB() #bFactor set in 1 for default! 

vectorC_UnitCost  <- as.vector( y$getUnitCost() )
vectorb_Target    <- as.vector( y$getTarget() )


boundarySize       <- nrow(dataBoundary)
indexI             <- 1 
indexJ             <- 1
connectivityCoeff  <- 0
connectivityVector <-  vector(mode = "numeric", length = numberUnits)
auxCoeff           <- 0
#matrixVij corresponds to the original quadrate matrix "V[i1,i2]", the symmetric one!
matrixVij = matrix(data = 0, nrow = numberUnits, ncol = numberUnits)


for(l in 1:boundarySize){
  indexI            <- dataBoundary[l, 1]
  indexJ            <- dataBoundary[l, 2]
  connectivityCoeff <- bFactor*getValueTuple_2(indexI, indexJ, dataBoundary)
  
  matrixVij[indexI, indexJ] = dataBoundary[l, 3]
  matrixVij[indexJ, indexI] = dataBoundary[l, 3]
  
  for(i in 1:numberUnits){
    if(indexI == i){
      auxCoeff              = connectivityVector[i]
      connectivityVector[i] = auxCoeff + connectivityCoeff
      #This allows you to specify in a set, the values associated with W[i] that are linked to Y[i1][i2]
      #Here, we suppose that the original quadrate matrix V[i1,i2], is symmetric!
    }
    if(indexJ == i){
      auxCoeff              = connectivityVector[i]
      connectivityVector[i] = auxCoeff + connectivityCoeff
    }
    
  }#END internal 'for'
}#END external 'for'


# Observation: vector with the associated coefficients of each variable W[i] 
#(number of positions equivalent to the number of planning units)
vectorC_PlanningCost     <- vectorC_UnitCost + connectivityVector

# Observation: vector with the associated coefficients of each variable Y[i1,i2] 
#(number of positions equivalent to the number of planning units)
vectorC_ConnectivityCost <- c();

#Creation of vectors containing the name of the decision variables (variables W[i] and Y[i1,i2]).
vectorVariablesW = c();
vectorVariablesY = c();
vectorIndexI_Yij = c();
vectorIndexJ_Yij = c();

for(i in 1:numberUnits){
  #Observation: "variableW()" Function was written in C++ language.
  varW = variableW(i) 
  vectorVariablesW = c(vectorVariablesW, varW)
  for(j in 1:numberUnits){
    if(i!=j && matrixVij[i,j] != 0 && bFactor > 0){
      #Observation: "variableY()" Function was written in C++ language.
      varY = variableY(i,j)
      
      vectorIndexI_Yij = c(vectorIndexI_Yij, i)
      vectorIndexJ_Yij = c(vectorIndexJ_Yij, j)
      vectorVariablesY = c(vectorVariablesY, varY)
      
      connectivityCoeff = -bFactor*matrixVij[i,j]
      vectorC_ConnectivityCost = c(vectorC_ConnectivityCost, connectivityCoeff)
    }
  }
}


#Number of variables W is equal to #N, where N is the set of {1,..., numberUnits}.
#Number of variables Y can only be calculated after creating the vector that contains 
#all the names associated with the valid variables Y[i1,i2]. Otherwise, the number of variables Y 
#is equal to #N*#N - #N, where N is the set of {1,..., numberUnits}.
numberVariablesW = numberUnits
numberVariablesY = length(vectorVariablesY)
numberVariables  = numberVariablesW + numberVariablesY
sizeRHS          = numberSpecies + (3*numberVariablesY)


vectorVariables = vector(mode = "character", length = numberVariables)
vectorC_Final   = vector(mode = "numeric", length = numberVariables)

auxIndex = numberVariablesW+1
vectorVariables[1:numberVariablesW]        = vectorVariablesW
vectorVariables[auxIndex:numberVariables]  = vectorVariablesY

vectorC_Final[1:numberVariablesW]          = vectorC_PlanningCost
vectorC_Final[auxIndex:numberVariables]    = vectorC_ConnectivityCost


matrixA_Target    <- matrix(NA, nrow = numberSpecies , ncol = numberUnits)
vectorA_Target    = c()

#Create matrix A for restriction of targets (encapsulate later in a generic function)
for(i in 1:numberSpecies){
  for(j in 1:numberUnits){
    matrixA_Target[i,j] = getValueTuple(j,i,dataDistribution) #Here, I used the original DataFrame named "dataDistribution".
    aux                 = as.integer( (matrixA_Target[i,j])*1 )
    vectorA_Target      = c(vectorA_Target, aux)
    #Observation: "getValueTuple()" Function was written in C++ language.
  }
}
matrixA_Target <- matrix(vectorA_Target, nrow = numberSpecies , ncol = numberUnits, byrow = TRUE)


matrixA_Final  <- matrix(data = 0, nrow = sizeRHS, ncol = numberVariables)
matrixA_Final[1:numberSpecies , 1:numberUnits] = matrixA_Target


#Linearity constraints (3 restrictions for each variable Y[i1,i2])
auxNumberUnits   = numberUnits + 1
#posY_Yij is a vector that stores the positions (indexes) of the variables Y[i1,i2] found in vectorVariables.
vectorPosY_Yij   = match(vectorVariablesY , vectorVariables)
vectorIndexI_Yij = vectorIndexI_Yij;
vectorIndexJ_Yij = vectorIndexJ_Yij;

for(k in 1:numberVariablesY){
  #varY         = variableY(i,j)
  posY_Yij     = vectorPosY_Yij[k]
  posX_Yij_ct1 = numberSpecies + (posY_Yij - auxNumberUnits)*3 + 1
  posX_Yij_ct2 = numberSpecies + (posY_Yij - auxNumberUnits)*3 + 2
  posX_Yij_ct3 = numberSpecies + (posY_Yij - auxNumberUnits)*3 + 3
  index_I      = vectorIndexI_Yij[k]
  index_J      = vectorIndexJ_Yij[k]
  
  #Constraint number 1 (Y[i1,i2] - X[i1] <= 0)
  matrixA_Final[ posX_Yij_ct1 , posY_Yij] =  1
  matrixA_Final[ posX_Yij_ct1 , index_I ] = -1
  
  #Constraint number 2 (Y[i1,i2] - X[i2] <= 0)
  matrixA_Final[ posX_Yij_ct2 , posY_Yij] =  1
  matrixA_Final[ posX_Yij_ct2 , index_J ] = -1  
  
  #Constraint number 3 (Y[i1,i2] - X[i1] - X[i2] => -1)
  matrixA_Final[ posX_Yij_ct3 , posY_Yij] =  1
  matrixA_Final[ posX_Yij_ct3 , index_I ] = -1
  matrixA_Final[ posX_Yij_ct3 , index_J ] = -1 
}



vectorb_LinearityConstraints   = rep( c(0,0,-1), numberVariablesY)
vectorb_Final                  = vector(mode = "numeric", length = sizeRHS)
vectorb_Final[1 : numberSpecies] = vectorb_Target
vectorb_Final[(numberSpecies+1):sizeRHS] = vectorb_LinearityConstraints


vectorSense_Target               = rep(">=",numberSpecies)
vectorSense_LinearityConstraints = rep( c("<=","<=",">="), numberVariablesY)
vectorSense                      = vector(mode = "character", length = sizeRHS)
vectorSense[1 : numberSpecies]   = vectorSense_Target
vectorSense[(numberSpecies+1):sizeRHS] = vectorSense_LinearityConstraints

vectorTypesVariables = rep("B", numberVariables)


# length(vectorC_UnitCost) #Size as to be equal to 10 (units).
# length(vectorb_Target)   #Size as to be equal to 5 (species).
# dim(matrixA_Target)      #Size as to be equal to 5*10 (species*units)


#Inputs parameters for Rsymphony's solver
obj      <- vectorC_Final
mat      <- matrixA_Final
dir      <- vectorSense
rhs      <- vectorb_Final

bounds   <- list(lower = list(ind = seq(1:numberVariables), val = rep(0,numberVariables)),
                 upper = list(ind = seq(1:numberVariables), val = rep(1,numberVariables)))

types    <- vectorTypesVariables #"C", "I",and "B" corresponding to continuous, integer, and binary variables.
max      <- FALSE                 #max <- TRUE means F.O = maximize / max <- FALSE means F.O = minimize
write_lp <- TRUE                  #Optional!
write_lp <- TRUE                  #Optional!


View(obj)
View(mat)
View(dir)
View(rhs)
View(bounds)
View(types)


modelSolver <- Rsymphony::Rsymphony_solve_LP(obj, mat, dir, rhs, bounds = bounds, 
               types = types, max=max, write_lp = write_lp, write_mps = write_lp)
print(modelSolver)



