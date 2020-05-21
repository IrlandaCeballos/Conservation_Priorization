# Instalación de paquetes
#------------------------------------------------------------
install.packages("purrr")
install.packages("Rcpp")
install.packages("RcppExamples")
install.packages("prioritizr")
install.packages("tidyverse")
#install.packages("prioritizrdata", repos = "https://cran.rstudio.com/")
#install.packages("BH")
#install.packages("Rsymphony")

# Activación de paquetes
#------------------------------------------------------------
library("purrr")
library("Rcpp")
library("RcppExamples")
#library("prioritizr")
library("tidyverse")
library("readr")
library("BH")
library("Rsymphony")

# Conexión con código C++
#------------------------------------------------------------
require( Rcpp )
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(file = "src/C_Script_Example.cpp")

# Modelo Matemático
#------------------------------------------------------------
setwd("C:/Users/Irlanda Ceballos/Documents/GitHub/Conservation_Priorization")
dataTarget       <- as.data.frame(read_delim("data/spec.csv", ";", trim_ws = TRUE)%>%select(1:5))
dataCost         <- as.data.frame(read_delim("data/PU.csv", ";", trim_ws = TRUE))
dataBoundary     <- as.data.frame(read_delim("data/bound.dat", "\t", escape_double = FALSE, trim_ws = TRUE))
dataDistribution <- as.data.frame(read_delim("data/puvspr2.dat", "\t", escape_double = FALSE, trim_ws = TRUE))
dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
colnames(dataDistribution) <- c("pu","species","amount")

#sólo para test! (funciona getValueTuple function development in C++)
#dataDistribution[2, 3] = 0;
#
# save(dataBoundary, file = "data/dataTarget.rda")
# save(dataBoundary, file = "data/dataCost.rda")
# save(dataBoundary, file = "data/dataBoundary.rda")
# save(dataBoundary, file = "data/dataDistribution.rda")
# 

Load(file = "data/dataTarget.rda")
# Load(file = "data/dataCost.rda")
# Load(file = "data/dataBoundary.rda")
# Load(file = "data/dataDistribution.rda")


sourceCpp(file = "src/RcppClass_DataInstance.cpp")
y <- new(DataInstance, dataTarget, dataCost, dataBoundary, dataDistribution)
y$print()
y$data01
y$data02
y$data03
y$data04

y$getB()
y$getSpecies()
y$getUnits()
y$getTarget()
y$getCost()
y$getDistribution()
y$getBoundary()

# test_1 <- y$getDistribution()
# test_2 <- y$getBoundary()
# test_1$`1`
# test_1[[1]]
# names(test_1)
# test_1 <- c(1,0,1,1,1,0)
# test_2 <- convertToBoolean(test_1)



# MARXAN con Rsymphony Solver (de la API de R)
#------------------------------------------------------------
library("Rcpp")
library("BH")
library("Matrix")
library("Rsymphony")

setwd("C:/Users/Irlanda Ceballos/Documents/GitHub/Conservation_Priorization")
dataTarget       <- as.data.frame(read_delim("data/data_Small/spec_Small.csv", ";", trim_ws = TRUE)%>%select(1:5))
dataCost         <- as.data.frame(read_delim("data/data_Small/PU_Small.csv", ";", trim_ws = TRUE))
dataBoundary     <- as.data.frame(read_delim("data/data_Small/bound_Small.csv", ";", trim_ws = TRUE))
dataDistribution <- as.data.frame(read_delim("data/data_Small/puvspr2_Small.csv", ";", trim_ws = TRUE))
dataDistribution <- cbind.data.frame(dataDistribution$pu, dataDistribution$species, dataDistribution$amount);
colnames(dataDistribution) <- c("pu","species","amount")

sourceCpp(file = "src/RcppClass_DataInstance.cpp")
sourceCpp(file = "src/RcppFunction_GlobalFunctions.cpp")
y <- new(DataInstance, dataTarget, dataCost, dataBoundary, dataDistribution)
numberUnits   <- y$getUnits()
numberSpecies <- y$getSpecies()
vectorC       <- as.vector( y$getCost() )
vectorb       <- as.vector( y$getTarget() )
#matrixA      <- y$getDistribution()

matrixA      <- matrix(NA, nrow = numberSpecies , ncol = numberUnits)
vectorA <- vector("numeric", numberSpecies*numberUnits)

#matrixA[1,1] <- 20 #indexes in R start from 1 (in C++, start from 0)
for(i in 1:numberSpecies){
  for(j in 1:numberUnits){
    #matrixA[i,j] <- getValueTuple(j,i,dataDistribution) #Here, I used the original DataFrame named "dataDistribution".
    #Observation: "getValueTuple()" Function was written in C++ language.
    vectorA[i,j] <- getValueTuple(j,i,dataDistribution) #Here, I used the original DataFrame named "dataDistribution".
    
  }
}

matrixA <- matrix(1*matrixA, nrow = numberSpecies , ncol = numberUnits)

length(vectorC) #Size as to be equal to 2316 (units).
length(vectorb) #Size as to be equal to 45 (species).
dim(matrixA)    #Size as to be equal to 45*2316 (species*units)

#Inputs parameters for Rsymphony's solver
obj      <- vectorC
mat      <- matrixA 
dir      <- rep(">=",numberSpecies)
rhs      <- vectorb
# bounds <- list(lower = list(ind = c(1L, 3L), val = c(-Inf, 2)),
#                upper = list(ind = c(1L, 2L), val = c(4, 100)))
bounds <- list(lower = list(ind = seq(1:numberUnits), val = rep(0,numberUnits)),
                upper = list(ind = seq(1:numberUnits), val = rep(1,numberUnits)))

types <- rep("B", numberUnits) #"C", "I",and "B" corresponding to continuous, integer, and binary variables.
max      <- FALSE #max <- TRUE means F.O = maximize / max <- FALSE means F.O = minimize
write_lp <- TRUE  #Optional!

Rsymphony::Rsymphony_solve_LP(obj, mat, dir, rhs, bounds = bounds, types = types, max=max, write_lp = write_lp)
#Rsymphony::Rsymphony_solve_LP(obj, mat, dir, rhs, types = types, max=max, write_lp = write_lp)



modelParameters <- list(obj, mat, dir, rhs, bounds, types, max)
# modelSolver     <- do.call(Rsymphony::Rsymphony_solve_LP, modelParameters)
# print(modelSolver)








# Ejemplos de conexión/interacción R con funciones escritas en C++
#------------------------------------------------------------
# Obteniendo constantes
evalCpp("PI")
evalCpp("std::numeric_limits<double>::max()")
evalCpp("__cplusplus")
#------------------------------------------------------------
# Datos
x   <- seq(10)
x_2 <- rep(4, 10)
x_3 <- as.list(x_2)
x_4 <- c(1, 1, 4, 5, 4, 6, 5, 6, 7)
y   <- matrix(sample(100),10)
z   <- 20 
t   <- c(FALSE,FALSE,FALSE)
mod <- lm(mpg~wt, data = mtcars)
mod
mod["residuals"]
mod["fitted.values"]
attributes(mod)
xx <- 2:18
v <- c(5, 10, 15) # create two bins [5,10) and [10,15)
#------------------------------------------------------------
# Ejecución de funciones en C++
test_copy()
foo(x)
bar(x)
rowSumsC(y)
pdistC(z, x)
pdistC2(z, x)
f1(x)
f2(x)
f3(t)
f5(x, x_2)
attribs() #Definiendo atributos a objetos
mpe(mod)  #Mean percentage error of a linear model  
callWithOne(function(x) x + 1)
callWithOne(paste)
hola <-lapply1(x_3, function(x) x + 1)
hola
is_naC(c(NA, 5.4, 3.2, NA))  #Retorna un vector logico con TRUE cuando encuentra elementos de un vector como NA.-
is_naC2(c(NA, 5.4, 3.2, NA)) #Retorna un vector logico con TRUE cuando encuentra elementos de un vector como NA.-
sum3(x)
sum4(x)
cbind(xx, findInterval(xx, v)) #Función "findInterval()" base de R
cbind(xx, findInterval2(xx,v)) #Función "findInterval()" desarrollada en C++
rleC(c(2,2,2,3,4,2,7,8,9,4))
duplicatedC(x_4) #no entiendo esta función como funciona :s
tableC(x_4)
#------------------------------------------------------------
# Ejecución de clases en C++, creación del objeto y acceso a sus métodos.
sourceCpp(file = "src/C_Script_ClassExample.cpp")
x <- new(Location, 10, 20)
x$print()

#------------------------------------------------------------
# Ejecución de funciones en C++, que utilizan boost libraries
sourceCpp(file = "src/RcppFunction_TestFunction.cpp")
computeGCD(10,100)
computeLCM(100,50)
