# Activación de paquetes
#------------------------------------------------------------------------------------------
sys.source(file = "Package_RLibraries.R") 


# Lectura de datos (inputs) de MARXAN - Importación de Clases/Functiones de C++ (Rcpp)
#------------------------------------------------------------------------------------------
setwd("C:/Users/Irlanda Ceballos/Documents/GitHub/Conservation_Priorization")
sourceCpp(file = "src/RcppClass_MARXANData.cpp")
sourceCpp(file = "src/RcppFunction_GlobalFunctions.cpp")

# #Lectura de datos (inputs) de MARXAN (inputs: 10 units/ 5 species/ 3 threats)
# dataTarget       <- as.data.frame(read_delim("data/data_ExtremelySmall/target_ExtremelySmall.csv", ";", trim_ws = TRUE)%>%select(1:5))
# dataCost         <- as.data.frame(read_delim("data/data_ExtremelySmall/unitCost_ExtremelySmall.csv", ";", trim_ws = TRUE))
# dataBoundary     <- as.data.frame(read_delim("data/data_ExtremelySmall/boundary_ExtremelySmall.csv", ";",
#                     trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
# dataDistribution <- as.data.frame(read_delim("data/data_ExtremelySmall/speciesDistribution_ExtremelySmall.csv", ";", trim_ws = TRUE))
# dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
# colnames(dataDistribution) <- c("pu","species","amount")


# # Lectura de datos (inputs) de MARXAN (inputs: 240 units/ 10 species/ 4 threats)
# dataTarget       <- as.data.frame(read_delim("data/data_Small/target_Small.csv", ";", trim_ws = TRUE)%>%select(1:5))
# dataCost         <- as.data.frame(read_delim("data/data_Small/unitCost_Small.csv", ";", trim_ws = TRUE))
# dataBoundary     <- as.data.frame(read_delim("data/data_Small/boundary_Small.csv", ";",
#                     trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
# dataDistribution <- as.data.frame(read_delim("data/data_Small/speciesDistribution_Small.csv", ";", trim_ws = TRUE))
# dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
# colnames(dataDistribution) <- c("pu","species","amount")
# 
# 
#Lectura de datos (inputs) de MARXAN (inputs: 2316 units/ 45 species/ 4 threats)
dataTarget       <- as.data.frame(read_delim("data/data_Big/target_Big.csv", ";", trim_ws = TRUE)%>%select(1:5))
dataCost         <- as.data.frame(read_delim("data/data_Big/unitCost_Big.csv", ";", trim_ws = TRUE))
dataBoundary     <- as.data.frame(read_delim("data/data_Big/boundary_Big.csv", ";",
                    trim_ws = TRUE, col_types = cols(id1 = col_integer(), id2 = col_integer(), boundary = col_double()  )))
dataDistribution <- as.data.frame(read_delim("data/data_Big/speciesDistribution_Big.csv", ";", trim_ws = TRUE))
dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
colnames(dataDistribution) <- c("pu","species","amount")


# MARXAN Básico (sin penalización) con Rsymphony Solver (de la API de R)
#------------------------------------------------------------------------------------------
y             <- new(MARXANData, dataTarget, dataCost, dataBoundary, dataDistribution)
numberUnits   <- y$getUnits()
numberSpecies <- y$getSpecies()
vectorC       <- as.vector( y$getUnitCost() )
vectorb       <- as.vector( y$getTarget() )
#matrixA      <- y$getSpeciesDistribution()
  
matrixA    <- matrix(NA, nrow = numberSpecies , ncol = numberUnits)
vectorA = c()

for(i in 1:numberSpecies){
  for(j in 1:numberUnits){
    matrixA[i,j] = getValueTuple(j,i,dataDistribution) #Here, I used the original DataFrame named "dataDistribution".
    aux          = as.integer( (matrixA[i,j])*1 )
    vectorA      = c(vectorA, aux)
    #Observation: "getValueTuple()" Function was written in C++ language.
  }
}

length(vectorC) #Size as to be equal to 2316 (units).
length(vectorb) #Size as to be equal to 45 (species).
dim(matrixA)    #Size as to be equal to 45*2316 (species*units)

matrixDEF <- matrix(vectorA, nrow = numberSpecies , ncol = numberUnits, byrow = TRUE)
matrixQLA <- matrix(matrixDEF , nrow = numberSpecies , ncol = numberUnits)

#Inputs parameters for Rsymphony's solver
obj      <- vectorC
mat      <- matrixQLA
dir      <- rep(">=",numberSpecies)
rhs      <- vectorb
#rhs      <- c(8,2,6,4,3) #Force the targets!
bounds   <- list(lower = list(ind = seq(1:numberUnits), val = rep(0,numberUnits)),
                 upper = list(ind = seq(1:numberUnits), val = rep(1,numberUnits)))
types    <- rep("B", numberUnits) #"C", "I",and "B" corresponding to continuous, integer, and binary variables.
max      <- FALSE                 #max <- TRUE means F.O = maximize / max <- FALSE means F.O = minimize
write_lp <- TRUE                  #Optional!
write_lp <- TRUE                  #Optional!

modelSolver <- Rsymphony::Rsymphony_solve_LP(obj, mat, dir, rhs, bounds = bounds, 
                              types = types, max=max, write_lp = write_lp, write_mps = write_lp)
print(modelSolver)

# modelParameters <- list(obj, mat, dir, rhs, bounds, types, max)
# modelSolver     <- do.call(Rsymphony::Rsymphony_solve_LP, modelParameters)
# print(modelSolver)

