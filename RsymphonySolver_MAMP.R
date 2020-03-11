# Activación de paquetes - Importación de Clases/Functiones de C++ (Rcpp)
#------------------------------------------------------------------------------------------
sys.source(file = "Package_RLibraries.R") 
setwd("C:/Users/Irlanda Ceballos/Documents/GitHub/Conservation_Priorization")
sourceCpp(file = "src/RcppClass_MAMPData.cpp")
sourceCpp(file = "src/RcppFunction_GlobalFunctions.cpp")

# problemData   <- new(MAMPData, target_Data, unitCost_Data, boundary_Data, speciesDistribution_Data, threatsDistribution_Data, sensibility_Data, actionCost_Data)
# subset_Ki_4 <- Ki[[1]]
# subset_Ki_3 <- Ki[[3]]
# intersect(subset_Ki_4,subset_Ki_3)


preprocessingTime_1 = Sys.time()

# Lectura de datos (inputs) de MARXAN/MAMP (inputs: 10 units/ 5 species/ 3 threats)
#------------------------------------------------------------------------------------------
target_Data       <- as.data.frame(read_delim("data/data_ExtremelySmall/spec_ExtremelySmall.csv", ";", trim_ws = TRUE)%>%select(1:5))
unitCost_Data     <- as.data.frame(read_delim("data/data_ExtremelySmall/PU_ExtremelySmall.csv", ";", trim_ws = TRUE))
boundary_Data     <- as.data.frame(read_delim("data/data_ExtremelySmall/bound_ExtremelySmall.csv", ";",
                                              trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
speciesDistribution_Data <- as.data.frame(read_delim("data/data_ExtremelySmall/puvspr2_ExtremelySmall.csv", ";", trim_ws = TRUE))
speciesDistribution_Data <- cbind.data.frame(speciesDistribution_Data$pu, speciesDistribution_Data$species, speciesDistribution_Data$amount);
colnames(speciesDistribution_Data) <- c("pu","species","amount")

#New data required by the MAMP problem (Threats distribution, species-threats sensibility and action cost)
threatsDistribution_Data <- as.data.frame( read_delim("data/data_ExtremelySmall/threatsDistribution_ExtremelySmall.csv",
                                                       ";", trim_ws = TRUE, col_types = cols(pu = col_integer(), threats = col_integer(), amount = col_integer() )) )
sensibility_Data <- as.data.frame( read_delim("data/data_ExtremelySmall/speciesThreatsSensibility_ExtremelySmall.csv",
                                               ";", trim_ws = TRUE, col_types = cols(species = col_integer(), threats = col_integer(), sensibility = col_integer() )) )
actionCost_Data <- as.data.frame( read_delim("data/data_ExtremelySmall/actionCost_ExtremelySmall.csv",
                                              ";", col_names = TRUE, trim_ws = TRUE, col_types = NULL ))


# # Lectura de datos (inputs) de MARXAN/MAMP (inputs: 2316 units/ 45 species/ 4 threats)
# #------------------------------------------------------------------------------------------
# target_Data       <- as.data.frame(read_delim("data/data_Big/spec_Big.csv", ";", trim_ws = TRUE)%>%select(1:5))
# unitCost_Data     <- as.data.frame(read_delim("data/data_Big/PU_Big.csv", ";", trim_ws = TRUE))
# boundary_Data     <- as.data.frame(read_delim("data/data_Big/bound_Big.csv", ";",
#                                              trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
# speciesDistribution_Data <- as.data.frame(read_delim("data/data_Big/puvspr2_Big.csv", ";", trim_ws = TRUE))
# speciesDistribution_Data <- cbind.data.frame(speciesDistribution_Data$pu, speciesDistribution_Data$species, speciesDistribution_Data$amount);
# colnames(speciesDistribution_Data) <- c("pu","species","amount")
# 
# #New data required by the MAMP problem (Threats distribution, species-threats sensibility and action cost)
# threatsDistribution_Data <- as.data.frame( read_delim("data/data_Big/threatsDistribution_Big.csv",
#                                                       ";", trim_ws = TRUE, col_types = cols(pu = col_integer(), threats = col_integer(), amount = col_integer() )) )
# sensibility_Data <- as.data.frame( read_delim("data/data_Big/speciesThreatsSensibility_Big.csv",
#                                               ";", trim_ws = TRUE, col_types = cols(species = col_integer(), threats = col_integer(), sensibility = col_integer() )) )
# actionCost_Data <- as.data.frame( read_delim("data/data_Big/actionCost_Big.csv",
#                                              ";", col_names = TRUE, trim_ws = TRUE, col_types = NULL ))


# # Lectura de datos (inputs) de MARXAN/MAMP (inputs: 240 units/ 10 species/ 4 threats)
# #------------------------------------------------------------------------------------------
# target_Data       <- as.data.frame(read_delim("data/data_Small/spec_Small.csv", ";", trim_ws = TRUE)%>%select(1:5))
# unitCost_Data     <- as.data.frame(read_delim("data/data_Small/PU_Small.csv", ";", trim_ws = TRUE))
# boundary_Data     <- as.data.frame(read_delim("data/data_Small/bound_Small.csv", ";",
#                                              trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
# speciesDistribution_Data <- as.data.frame(read_delim("data/data_Small/puvspr2_Small.csv", ";", trim_ws = TRUE))
# speciesDistribution_Data <- cbind.data.frame(speciesDistribution_Data$pu, speciesDistribution_Data$species, speciesDistribution_Data$amount);
# colnames(speciesDistribution_Data) <- c("pu","species","amount")
# 
# #New data required by the MAMP problem (Threats distribution, species-threats sensibility and action cost)
# threatsDistribution_Data <- as.data.frame( read_delim("data/data_Small/threatsDistribution_Small.csv",
#                                                       ";", trim_ws = TRUE, col_types = cols(pu = col_integer(), threats = col_integer(), amount = col_integer() )) )
# sensibility_Data <- as.data.frame( read_delim("data/data_Small/speciesThreatsSensibility_Small.csv",
#                                               ";", trim_ws = TRUE, col_types = cols(species = col_integer(), threats = col_integer(), sensibility = col_integer() )) )
# actionCost_Data <- as.data.frame( read_delim("data/data_Small/actionCost_Small.csv",
#                                              ";", col_names = TRUE, trim_ws = TRUE, col_types = NULL ))

#------------------------------------------------------------------------------------------
#------------- # MAMP (with gamma = 1), con Rsymphony Solver (de la API de R) -------------
#------------------------------------------------------------------------------------------
#Observation: "MAMPData" Class was writen in the C++ language.
problemData   = new(MAMPData, target_Data, unitCost_Data, boundary_Data, speciesDistribution_Data, threatsDistribution_Data, sensibility_Data, actionCost_Data)

numberUnits   = problemData$getUnits()
numberSpecies = problemData$getSpecies()
numberThreats = problemData$getThreats()
actionCost    = problemData$getActionCost()
# bFactor       = problemData$getB() #bFactor set in 1 for default! 
bFactor       = 1

#Si = problemData$getSet("Si") #It corresponds to the set Si (Si belong to set S).
Ki = problemData$getSet("Ki") #It corresponds to the set Ki (Ki belong to set K).
#Ks = problemData$getSet("Ks") #It corresponds to the set Ks (Ks belong to set K).
#Is = problemData$getSet("Is") #It corresponds to the set Is (Is belong to set I).




#matrix_cv corresponds to the original quadrate matrix "cv[i1,i2]", the symmetric one!
#matrix_speciesDistribution, matrix_threatsDistribution and matrix_c corresponds to the original matrices
#(that is, not the dense matrix in its tuple version).
#Observation: "originalMatrix_cv()", "originalMatrix_Distribution()" and "createMatrix_c()" functions were written in the C++ language.
matrix_cv <- originalMatrix_cv(type="symmetric", data= boundary_Data, units= numberUnits) 
matrix_speciesDistribution <- originalMatrix_Distribution(data= speciesDistribution_Data, units= numberUnits, species= numberSpecies)
matrix_threatsDistribution <- originalMatrix_Distribution(data= threatsDistribution_Data, units= numberUnits, species= numberThreats)
matrix_c <- createMatrix_c(dataDistribution = matrix_threatsDistribution, dataCost = as.vector(x = unitCost_Data[,2] , mode = "numeric"))

# write.csv(matrix_cv, file="NEW_cv_Matrix.csv", row.names = F)
# write.csv(matrix_speciesDistribution, file="NEW_speciesDistribution_Matrix.csv", row.names = F)
# write.csv(matrix_threatsDistribution, file="NEW_threatsDistribution_Matrix.csv", row.names = F)
# write.csv(matrix_c, file="NEW_actionCost_Matrix.csv", row.names = F)

#------------------------------------------------------------------------------------------
#---------------- Vector C - Coefficients associated with the planning cost ---------------
#--------------------- (coefficients associated with W[i] variables) ----------------------
#------------------------------------------------------------------------------------------
boundarySize       <- nrow(boundary_Data)
index_i1           <- 1 
index_i2           <- 1
connectivityCoeff  <- 0
connectivityVector <- vector(mode = "numeric", length = numberUnits)
auxCoeff           <- 0
for(l in 1:boundarySize){#OJO!:VER COMO OPTIMIZAR ESTA PARTE CUANDO bFactor == 0
  index_i1            <- boundary_Data[l, 1]
  index_i2            <- boundary_Data[l, 2]
  connectivityCoeff   <- bFactor*getValueTuple_2(index_i1, index_i2, boundary_Data)
  
  for(i in 1:numberUnits){
    #Here we assume that the "cv[i1,i2]" matrix is symmetric!
    if(index_i1 == i){
      auxCoeff              = connectivityVector[i]
      connectivityVector[i] = auxCoeff + connectivityCoeff
      #This allows you to specify in a set, the values associated with W[i] that are linked to Y[i1][i2]
      #Here, we suppose that the original quadrate matrix cm[i1,i2], is symmetric!
    }
    if(index_i2 == i){
      auxCoeff              = connectivityVector[i]
      connectivityVector[i] = auxCoeff + connectivityCoeff
    }
  }#END internal 'for'
}#END external 'for'

#Observation: "vectorC_PlanningCost", vector with the associated coefficients of each variable W[i], 
#(number of positions equivalent to the number of planning units)
vectorC_UnitCost      <- as.vector( x = unitCost_Data[,2] , mode = "numeric" )
vectorC_PlanningCost  <- vectorC_UnitCost + connectivityVector


#------------------------------------------------------------------------------------------
#-------------- Vector C - Coefficients associated with the connectivity cost -------------
#-------------------- (coefficients associated with Y[i1,i2] variables) -------------------
#------------------------------------------------------------------------------------------
# Observation: vector with the associated coefficients of each variable Y[i1,i2] 
#(number of positions equivalent to the number of planning units)
vectorC_ConnectivityCost <- c()

#Creation of vectors containing the name of the decision variables (variables W[i] and Y[i1,i2]).
#Observation: "createVectorVarW()" Function was written in C++ language.
vectorVariablesW   = createVectorVarW(units = numberUnits);
vectorVariablesY   = c()
vectorIndex_i1_Yij = c()
vectorIndex_i2_Yij = c()

#Important observation!
#If bFactor = 0, the Y[i1,i2] variables and the linearization restrictions associated with these variables are not required!
if(bFactor != 0){
  for(i in 1:numberUnits){
    for(j in 1:numberUnits){
      if(i!=j && matrix_cv[i,j] != 0){
        #Observation: "variableY()" Function was written in C++ language.
        varY = variableY(i,j)
        
        vectorIndex_i1_Yij = c(vectorIndex_i1_Yij, i)
        vectorIndex_i2_Yij = c(vectorIndex_i2_Yij, j)
        vectorVariablesY   = c(vectorVariablesY, varY)
        
        connectivityCoeff = -bFactor*matrix_cv[i,j]
        vectorC_ConnectivityCost = c(vectorC_ConnectivityCost, connectivityCoeff)
      }
    }
  }
}#END If


#------------------------------------------------------------------------------------------
#----------------- Vector C - Coefficients associated with the action cost ----------------
#--------------------- (coefficients associated with X[i,k] variables) --------------------
#------------------------------------------------------------------------------------------
vectorC_ActionCost = c()
vectorVariablesX   = c()
vectorIndexI_Xik   = c()
vectorIndexK_Xik   = c()
actionCostCoeff    = 0

for(i in 1:numberUnits){ 
  for(k in 1:numberThreats){
    if(matrix_c[k,i] > 0){
      varX = variableX(i, k)
      
      vectorVariablesX = c(vectorVariablesX, varX)
      vectorIndexI_Xik = c(vectorIndexI_Xik, i)
      vectorIndexK_Xik = c(vectorIndexK_Xik, k)
      
      actionCostCoeff  = matrix_c[k,i]
      vectorC_ActionCost = c(vectorC_ActionCost, actionCostCoeff)
    } 
  }
}

#------------------------------------------------------------------------------------------
#------------ Vector C - Coefficients associated with the auxiliary variable Z ------------
#--------------------- (coefficients associated with Z[i,s] variables) --------------------
#------------------------------------------------------------------------------------------
Si <- problemData$getSpeciesDistribution() #It corresponds to the set Si (Si belong to set S).
#Ki <- problemData$getThreatsDistribution() #It corresponds to the set Ki (Ki belong to set K).
Ks <- problemData$getSensibility()         #It corresponds to the set Ks (Ks belong to set K).

vectorVariablesZ  = c()
vectorIndex_i_Zis = c()
vectorIndex_s_Zis = c()


for(i in 1:length(Si)){
  i_var <- as.integer(names(Si[i]))
  subset_Si <- as.vector(names(Si[[i]]), mode = "integer")
  if(i!=i_var){
    i=i_var
  }
  if(numberUnits!=length(Ki)){
    i=1
     k_var <- match(i,names(Ki))
     subset_Ki    <- as.vector(names(Ki[[k_var]]), mode = "integer")
  }
  else{
    subset_Ki    <- as.vector(names(Ki[[i]]), mode = "integer")
  }
  
  for(s in subset_Si){
    subset_Ks    <- as.vector(names(Ks[[s]]), mode = "integer")
    subset_Ki_Ks <- intersect(subset_Ki,subset_Ks)
    if(length(subset_Ki_Ks) == 0){
      varZ = variableZ(i,s)
      vectorVariablesZ  = c(vectorVariablesZ, varZ)
      vectorIndex_i_Zis = c(vectorIndex_i_Zis, i) 
      vectorIndex_s_Zis = c(vectorIndex_s_Zis, s)
    }
  }
}


numberVariablesZ = length(vectorVariablesZ)
vectorC_auxVarZ  = vector(mode = "numeric", length = numberVariablesZ)

#Number of variables W is equal to #N, where N is the set of {1,..., numberUnits}.
#Number of variables Y can only be calculated after creating the vector that contains 
#all the names associated with the valid variables Y[i1,i2]. Otherwise, the number of variables Y 
#is equal to #N*#N - #N, where N is the set of {1,..., numberUnits}.
#Number of variables X can only be calculated after creating the vector that contains
#all the names associated with the valid variables X[i,k].
numberVariablesW = numberUnits
numberVariablesY = length(vectorVariablesY)
numberVariablesX = length(vectorVariablesX)
numberVariables  = numberVariablesW + numberVariablesY + numberVariablesX + numberVariablesZ
#sizeRHS          = numberSpecies + (3*numberVariablesY)


vectorVariables = c(vectorVariablesW, vectorVariablesY, vectorVariablesX, vectorVariablesZ)
vectorC_Final   = as.vector(x= c(vectorC_PlanningCost, vectorC_ConnectivityCost, vectorC_ActionCost, vectorC_auxVarZ), mode = "numeric")	

#vectorC_Final_Sparse = Matrix(data = vectorC_Final, nrow=1, sparse = TRUE)
vectorC_Final_Sparse = as(vectorC_Final, "sparseVector")

object.size(vectorC_Final)
object.size(vectorC_Final_Sparse)

#View(vectorC_Final)  

#------------------------------------------------------------------------------------------
#----- Matrix A - Coefficients associated with the variables of the first restriction -----
#--------------------------------------- (MAMP. 2) ----------------------------------------
#------------------------------------------------------------------------------------------
sizeRHS_MAMP2 = numberSpecies
#matrixA_MAMP2 <- matrix(data = 0, nrow = sizeRHS_MAMP2, ncol = numberVariables)
matrixA_MAMP2 <- Matrix(data = 0, nrow = sizeRHS_MAMP2, ncol = numberVariables, sparse = TRUE)

# vectorIndexI_Xik  #Subindex of "i" variable X (creo que no se usa! no es necesario!)
# vectorIndexK_Xik  #Subindex of "k" variable X 
# vectorIndex_i_Zis #Subindex of "i" variable Z 
# vectorIndex_s_Zis #Subindex of "s" variable Z 

Is = problemData$getSpeciesDistribution_t()
posX_Xik = 0
posY_Xik = 0
posX_Zis = 0
posY_Zis = 0

#Break points (x axis) of vector C, which indicate the separations of each set of variables that are within that vector.
break_1 = 0
break_2 = numberVariablesW
break_3 = numberVariablesW + numberVariablesY
break_4 = numberVariablesW + numberVariablesY + numberVariablesX
break_5 = numberVariables

for(s in 1:numberSpecies){
  subset_Is <- as.vector(names(Is[[s]]), mode = "integer")
  
  for(i in subset_Is){
      if(is.na(match(i,as.vector(names(Ki), mode = "integer")))){
        intersectionSize <- 0
      }
      else{
        subset_Ks        <- as.vector(names(Ks[[s]]), mode = "integer")
        subset_Ki        <- as.vector(names(Ki[[i]]), mode = "integer")
        subset_Ki_Ks     <- intersect(subset_Ki,subset_Ks)
        intersectionSize <- length(subset_Ki_Ks)
        print(paste("unit ->" , i))
        print(paste("subset_Ki ->" , subset_Ki) )
      }
  
      if(intersectionSize != 0){
        coeff_Xik = 1/intersectionSize
        for(k in subset_Ki_Ks){
          varX     = variableX(i,k)
          posX_Xik = break_3 + match(varX, vectorVariablesX) 
          posY_Xik = s
          matrixA_MAMP2[ posY_Xik, posX_Xik] =  coeff_Xik
        }#END internal for (about Threats!)
        
      }else{#When |Ki intersect Ks| = 0
        coeff_Zis = 1
        varZ      = variableZ(i,s)
        posX_Zis  = break_4 + match(varZ, vectorVariablesZ)
        posY_Zis  = s
        matrixA_MAMP2[ posY_Zis, posX_Zis] = coeff_Zis
      }#END internal if-else 
  }#END external_1 for (about Planning Units!)
}#END external_2 for (about Species!)

#View(matrixA_MAMP2[ , (break_3+1):numberVariables])

#------------------------------------------------------------------------------------------
#---- Matrix A - Coefficients associated with the variables of the second restriction -----
#--------------------------------------- (MAMP. 3) ----------------------------------------
#------------------------------------------------------------------------------------------
sizeRHS_MAMP3 = length(Ki)
#matrixA_MAMP3 <- matrix(data = 0, nrow = sizeRHS_MAMP3, ncol = numberVariables)
matrixA_MAMP3 <- Matrix(data = 0, nrow = sizeRHS_MAMP3, ncol = numberVariables, sparse = TRUE)

posX_Wi  = 0
posY_Wi  = 0
posX_Xik = 0
posY_Xik = 0


for(i in 1:length(Ki)){
  k_var <- as.integer(names(Ki[i]))
  subset_Ki <- as.vector(names(Ki[[i]]), mode = "integer") 
  coeff_Wi  <- -1*length(subset_Ki)
  if(i!=k_var){
    i=k_var
  }
  for(k in subset_Ki){
    posX_Wi  = i
    posY_Wi  = i
    matrixA_MAMP3[posY_Wi, posX_Wi] = coeff_Wi
    
    coeff_Xik = 1
    varX      = variableX(i,k)
    posX_Xik  = break_3 + match(varX, vectorVariablesX) 
    posY_Xik  = i
    matrixA_MAMP3[posY_Xik, posX_Xik] = coeff_Xik
  }
}

#view(matrixA_MAMP3)

#------------------------------------------------------------------------------------------
#---- Matrix A - Coefficients associated with the variables of the third restriction ------
#--------------------------------------- (MAMP. 4) ----------------------------------------
#------------------------------------------------------------------------------------------
if(numberVariablesZ != 0){
  sizeRHS_MAMP4 = length(unique(vectorIndex_i_Zis))
  #matrixA_MAMP4 <- matrix(data = 0, nrow = sizeRHS_MAMP4, ncol = numberVariables)
  matrixA_MAMP4 <- Matrix(data = 0, nrow = sizeRHS_MAMP4, ncol = numberVariables, sparse = TRUE)
  
  posX_Wi  = 0
  posY_Wi  = 0
  posX_Zis = 0
  posY_Zis = 0
  
  for(i in unique(vectorIndex_i_Zis)){
    subset_Si <-which(vectorIndex_i_Zis %in% i)
    
    coeff_Wi <- -1*length(subset_Si)
    posX_Wi  = i
    posY_Wi  = which(unique(vectorIndex_i_Zis) %in% i)
    matrixA_MAMP4[posY_Wi, posX_Wi] = coeff_Wi
    for(s in subset_Si){
      coeff_Zis = 1
      posX_Zis  = break_4 + s
      posY_Zis  = posY_Wi
      matrixA_MAMP4[posY_Zis, posX_Zis] = coeff_Zis
    }
  }
  
  
  #View(matrixA_MAMP4)
}else{
  sizeRHS_MAMP4 = 0
  matrixA_MAMP4 <- Matrix(data = 0, nrow = sizeRHS_MAMP4, ncol = numberVariables, sparse = TRUE)
} 



#------------------------------------------------------------------------------------------
# Matrix A - Coefficients associated the linearisation of the MAMP model objective function 
#-------------- (3 restrictions for each variable Y[i1,i2]) -------------------------------
#------------------------------------------------------------------------------------------
if(bFactor != 0){
  #Linearity constraints for the fragmentation of planning units.
  sizeRHS_MAMP6 = 3*numberVariablesY 
  #matrixA_MAMP6 <- matrix(data = 0, nrow = sizeRHS_MAMP6, ncol = numberVariables)
  matrixA_MAMP6 <- Matrix(data = 0, nrow = sizeRHS_MAMP6, ncol = numberVariables, sparse = TRUE)
  
  auxNumberUnits = numberUnits + 1
  posX_Yij      = 0 #I used the name "Yij" instead of the original one, "Yi1i2".
  posY_Yij_ct1  = 0
  posY_Yij_ct2  = 0
  posY_Yij_ct3  = 0
  
  for(l in 1:numberVariablesY){
    index_i1 = vectorIndex_i1_Yij[l]
    index_i2 = vectorIndex_i2_Yij[l]
    varY     = variableY(index_i1, index_i2)
    #print(varY)
    
    posX_Yij     = break_2 + match(varY, vectorVariablesY)
    #print(posX_Yij)
    posY_Yij_ct1 = (posX_Yij - auxNumberUnits)*3 + 1
    posY_Yij_ct2 = (posX_Yij - auxNumberUnits)*3 + 2
    posY_Yij_ct3 = (posX_Yij - auxNumberUnits)*3 + 3
    
    #Constraint number 1 (Y[i1,i2] - X[i1] <= 0)
    matrixA_MAMP6[ posY_Yij_ct1 , posX_Yij ] =  1
    matrixA_MAMP6[ posY_Yij_ct1 , index_i1 ] = -1
    
    #Constraint number 2 (Y[i1,i2] - X[i2] <= 0)
    matrixA_MAMP6[ posY_Yij_ct2 , posX_Yij ] =  1
    matrixA_MAMP6[ posY_Yij_ct2 , index_i2 ] = -1  
    
    #Constraint number 3 (Y[i1,i2] - X[i1] - X[i2] => -1)
    matrixA_MAMP6[ posY_Yij_ct3 , posX_Yij ] =  1
    matrixA_MAMP6[ posY_Yij_ct3 , index_i1 ] = -1
    matrixA_MAMP6[ posY_Yij_ct3 , index_i2 ] = -1 
  }
  #View(matrixA_MAMP6)
  matrixA_Final = rbind(matrixA_MAMP2, matrixA_MAMP3, matrixA_MAMP4, matrixA_MAMP6)
}else{
  matrixA_Final = rbind(matrixA_MAMP2, matrixA_MAMP3, matrixA_MAMP4)
} 

#View(matrixA_Final)
matrixA_Final_Sparse = Matrix(data = matrixA_Final, sparse = TRUE)

object.size(matrixA_Final)/(1024*1024)
object.size(matrixA_Final_Sparse)/(1024*1024)

#------------------------------------------------------------------------------------------
#----------------- Vector b - Parameters associated with all restrictions -----------------
#-------------------------------- (RHS, right-hand side) ----------------------------------
#------------------------------------------------------------------------------------------
sizeRHS = sizeRHS_MAMP2 + sizeRHS_MAMP3 + sizeRHS_MAMP4
ts      = problemData$getTarget() #Conservation target of each species!

vectorb_MAMP2 = as.vector(x = ts, mode = "numeric")
vectorb_MAMP3 = vector(mode = "numeric", length = sizeRHS_MAMP3)
vectorb_MAMP4 = vector(mode = "numeric", length = sizeRHS_MAMP4)
vectorb_MAMP6 = rep(c(0,0,-1), numberVariablesY)

vectorb_Final = c(vectorb_MAMP2, vectorb_MAMP3, vectorb_MAMP4, vectorb_MAMP6)

#------------------------------------------------------------------------------------------
#------------------------ Others inputs parameters for the solver -------------------------
#------------ (Vector of sense, Vector of bounds, and Vector of variables type ) ----------
#------------------------------------------------------------------------------------------
vectorSense_MAMP2 = rep(">=", sizeRHS_MAMP2)
vectorSense_MAMP3 = rep("<=", sizeRHS_MAMP3)
vectorSense_MAMP4 = rep("<=", sizeRHS_MAMP4)
vectorSense_MAMP6 = rep( c("<=","<=",">="), numberVariablesY)

vectorSense_Final = c(vectorSense_MAMP2, vectorSense_MAMP3, vectorSense_MAMP4, vectorSense_MAMP6)

vectorBounds <- list(lower = list(ind = seq(1:numberVariables), val = rep(0,numberVariables)),
                     upper = list(ind = seq(1:numberVariables), val = rep(1,numberVariables)))

vectorVariablestype = rep("B", numberVariables) 


preprocessingTime_2 = Sys.time() ; preprocessingTime = preprocessingTime_2 - preprocessingTime_1

#------------------------------------------------------------------------------------------
#------------------------------- Solvers (Rsymphony - GLPK) -------------------------------
#------------------------------------------------------------------------------------------
processingTime_1 = Sys.time()

#Inputs parameters for solvers (matrix notation of the MAMP mathematical model)
#obj      <- vectorC_Final
#mat      <- matrixA_Final
obj      <- vectorC_Final_Sparse
mat      <- matrixA_Final_Sparse
dir      <- vectorSense_Final
rhs      <- vectorb_Final
bounds   <- vectorBounds
types    <- vectorVariablestype #"C", "I",and "B" corresponding to continuous, integer, and binary variables.
max      <- FALSE               #max <- TRUE means F.O = maximize / max <- FALSE means F.O = minimize
write_lp <- TRUE                #Optional!
write_lp <- TRUE                #Optional!



#Rsymphony Solver!
modelSolver_Rsymphony <- Rsymphony::Rsymphony_solve_LP(obj, mat, dir, rhs, bounds = bounds, 
                                             types = types, max=max, verbosity = 100, write_lp = write_lp, write_mps = write_lp)


print(modelSolver_Rsymphony)


processingTime_2 = Sys.time() ; processingTime = processingTime_2 - processingTime_1
print("Preprocessing Time (matrix construction) = ")
print(preprocessingTime)
print("Processing Time (solver execution) = ")
print(processingTime)


# #GLPK Solver!
# modelSolver_GLPK <- Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, bounds = bounds,
#                                           types = types, max=max)
# print(modelSolver_GLPK)





