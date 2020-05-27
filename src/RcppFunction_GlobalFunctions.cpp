#include "Package.h"
/* 
 * Function that converts any vector of binary integers (only 1 and 0) into 
 * a logical vector (only true and false).
 * Function only needs the vector you want to transform.
 */
// [[Rcpp::export(name = "convertToBoolean")]]
LogicalVector convertToBoolean(IntegerVector x){
  int dataSize = x.length();
  LogicalVector y(dataSize); 
  
  for(int i=0; i<dataSize; i++){
    if(x[i] == 1){
      y[i] = true;
    }else{
      y[i] = false;
    }
  }
  return y; 
} 

/* 
 * If the values of the r-th row of columns 1 and 2 of a 2D array correspond to 
 * the specified coordinate (i, j), the value of the r-th row and 3rd column of said array is obtained.
 */
// [[Rcpp::export(name = "getValueTuple")]]
bool getValueTuple(int indexI, int indexJ, DataFrame arrayData){
  bool value;
  int dataSize = arrayData.nrows();
  IntegerVector vectorAux08 = arrayData[0];
  IntegerVector vectorAux09 = arrayData[1];
  LogicalVector vectorAux10 = arrayData[2]; 
  
  for(int i = 0; i < dataSize; i++){
    if( vectorAux08[i] == indexI  && vectorAux09[i] == indexJ  ){
      value = vectorAux10[i];
      break;
    }else{
      value = false;
    }
  }
  return value;
}


/* 
 * If the values of the r-th row of columns 1 and 2 of a 2D array correspond to 
 * the specified coordinate (i, j), the value of the r-th row and 3rd column of said array is obtained.
 */
// [[Rcpp::export(name = "getValueTuple_2")]]
double getValueTuple_2(int indexI, int indexJ, DataFrame arrayData){
  double value;
  int dataSize = arrayData.nrows();
  IntegerVector vectorAux08 = arrayData[0];
  IntegerVector vectorAux09 = arrayData[1];
  NumericVector vectorAux10 = arrayData[2];
  
  for(int i = 0; i < dataSize; i++){
    if( vectorAux08[i] == indexI  && vectorAux09[i] == indexJ  ){
      value = vectorAux10[i];
      break;
    }else{
      value = 0;
    }
  }
  return value;
}

/* 
 * Function sets the name of a variable "W" with is respective index "i" (i.e, "W[i]").
 */
// [[Rcpp::export(name = "variableW")]]
String variableW(int indexI){
  std::string index_I       = std::to_string(indexI);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string varW          = "W" + parenthesis_1 + index_I + parenthesis_2;
  return varW;
}

/* 
 * Function that creates a vector with the name of all variables W, 
 * (i.e, vectorVariablesW = {"W[1]","W[2]", ..., "W[n]"}).
 */
// [[Rcpp::export(name = "createVectorVarW")]]
CharacterVector createVectorVarW(int units){
  CharacterVector vectorVariablesW;
  std::string varW;
  
  for(int i = 0; i < units; i++){
    varW = variableW(i+1);
    vectorVariablesW.push_back(varW); 
    //I use ".push_back()" because I don't prefix the length of the vector "vectorVariablesW".
}
  return vectorVariablesW;
}

/*
 * Function sets the name of a variable "Y" with is respective index "i1" and "i2" (i.e, "Y[i1,i2]").
 */
// [[Rcpp::export(name = "variableY")]]
String variableY(int indexI, int indexJ){
  std::string index_I       = std::to_string(indexI);
  std::string index_J       = std::to_string(indexJ);
  std::string i_NameIndex   = "i";
  std::string parenthesis_1 = "[";  
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  //std::string varY          = "Y" + parenthesis_1 + i_NameIndex + index_I + comma + i_NameIndex + index_J + parenthesis_2; 
  std::string varY          =  "Y" + parenthesis_1 + index_I + comma + index_J + parenthesis_2; 
  return varY;
}

/* 
 * Function sets the name of a variable "X" with is respective index "i" and "k" (i.e, "X[i,k]").
 */
// [[Rcpp::export(name = "variableX")]]
String variableX(int indexI, int indexK){
  std::string index_I       = std::to_string(indexI);
  std::string index_K       = std::to_string(indexK);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  std::string varX          = "X" + parenthesis_1 + index_I + comma + index_K + parenthesis_2;
  return varX;
}

/* 
 * Function sets the name of a variable "Z" with is respective index "i" and "s" (i.e, "Z[i,s]").
 */
// [[Rcpp::export(name = "variableZ")]]
String variableZ(int indexI, int indexS){
  std::string index_I       = std::to_string(indexI);
  std::string index_S       = std::to_string(indexS);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  std::string varZ          = "Z" + parenthesis_1 + index_I + comma + index_S + parenthesis_2;
  return varZ;
}

/* 
 * Function sets the name of a variable "B" with is respective index "i" and "s" (i.e, "B[i,s]").
 */
// [[Rcpp::export(name = "variableB")]]
String variableB(int indexI, int indexS){
  std::string index_I       = std::to_string(indexI);
  std::string index_S       = std::to_string(indexS);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  std::string varB          = "B" + parenthesis_1 + index_I + comma + index_S + parenthesis_2;
  return varB;
}

/* 
 * Function sets the name of a variable "Lambda" with is respective sub-index "i", "s" and "m" (i.e, "Lambda[i,s,m]").
 */
// [[Rcpp::export(name = "variableLambda")]]
String variableLambda(int indexI, int indexS, int indexM){
  std::string index_I       = std::to_string(indexI);
  std::string index_S       = std::to_string(indexS);
  std::string index_M       = std::to_string(indexM);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  std::string varLambda     = "Lambda" + parenthesis_1 + index_I + comma + index_S + comma + index_M + parenthesis_2;
  return varLambda;
}

/* 
 * Function sets the name of a variable "V" with is respective sub-index "i", "s" and "m" (i.e, "V[i,s,m]").
 */
// [[Rcpp::export(name = "variableV")]]
String variableV(int indexI, int indexS, int indexM){
  std::string index_I       = std::to_string(indexI);
  std::string index_S       = std::to_string(indexS);
  std::string index_M       = std::to_string(indexM);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  std::string varV          = "V" + parenthesis_1 + index_I + comma + index_S + comma + index_M + parenthesis_2;
  return varV;
}

/*
 * Function sets the name of a variable "P" with is respective sub-index "i1", "i2" and "k" (i.e, "P[i1,i2,k]").
 */
// [[Rcpp::export(name = "variableP")]]
String variableP(int indexI1, int indexI2, int indexK){
  std::string index_I1      = std::to_string(indexI1);
  std::string index_I2      = std::to_string(indexI2);
  std::string index_K       = std::to_string(indexK);
  std::string parenthesis_1 = "[";
  std::string parenthesis_2 = "]";
  std::string comma         = ",";
  std::string varP          = "P" + parenthesis_1 + index_I1 + comma + index_I2 + comma + index_K + parenthesis_2;
  return varP;
}


/*
 * Function that creates the original matrix "cv[i1, i2]" from the information of the 
 * dataframe "boundary_Data" (dense matrix with tuple data). 
 * Function only needs the number of planning units, the type the data distribution 
 * that you want (symmetric or asymmetric), and the data (tuple data), obviously.
 */
// [[Rcpp::export(name = "originalMatrix_cv")]]
NumericMatrix originalMatrix_cv(String type, DataFrame data, int units){
  IntegerVector vectorAux01 = data[0]; //id1
  IntegerVector vectorAux02 = data[1]; //id2
  NumericVector vectorAux03 = data[2]; //boundary
  NumericMatrix matrix_cv(units, units); //nrow = units, ncol = units (quadrate matrix)
  
  int dataSize = data.nrows();
  int index_i1 = 0;
  int index_i2 = 0;
  double value;
  
  for(int l = 0; l < dataSize; l++){
    index_i1 = vectorAux01[l] -1;
    index_i2 = vectorAux02[l] -1;
    value    = vectorAux03[l];
    
    // Rcpp::Rcout<<"indexI: "<<index_i1<<std::endl;
    // Rcpp::Rcout<<"indexJ: "<<index_i2<<std::endl;
    // Rcpp::Rcout<<"value: "<<value<<std::endl;

    if(type == "symmetric" ){
      matrix_cv(index_i1, index_i2) = value;
      matrix_cv(index_i2, index_i1) = value;
    }
    if(type == "asymmetric" ){
      matrix_cv(index_i1, index_i2) = value;
    }
  }

  
  
  
    // if(type == "symmetric" ){
    //   matrix_cv(index_i1, index_i2) = value;
    // }}
  
  
  return matrix_cv;
}

/*
 * Function that creates the original matrix (sparse matrix) from the information of the data frames: 
 * "speciesDistribution_Data" or "threatsDistribution_Data" (dense matrices with tuple data),
 * depending on what the user wants.
 * Function only needs the number of planning units, the number of species (or threats) 
 * and the data (tuple data), obviously. 
 */
// [[Rcpp::export(name = "originalMatrix_Distribution")]]
IntegerMatrix originalMatrix_Distribution(DataFrame data, int units, int species){
  IntegerVector vectorAux01 = data[0]; //pu
  IntegerVector vectorAux02 = data[1]; //species (or threats)
  IntegerVector vectorAux03 = data[2]; //amount
  IntegerMatrix matrix_Distribution(species, units); //nrow = species (or threats), ncol = units
  
  int dataSize =  data.nrows();
  int index_i  = 0; //index for units
  int index_s  = 0; //index for species (or threats)
  double value;
  
  for(int l = 0; l < dataSize; l++){
    index_i = vectorAux01[l] -1;
    index_s = vectorAux02[l] -1;
    value   = vectorAux03[l];
    matrix_Distribution(index_s, index_i) = value;
  }
return matrix_Distribution;  
}

/*
 * Function that creates matrix "c[i,k]" from vector "cm[i]", assuming that the cost of applying an action to eliminate 
 * threat k in K from unit i in I is the same for all threats that are in i (set Ki).
 * Function only needs the original matrix of threats distribution and the vector of planning unit cost.
 */
// [[Rcpp::export(name = "createMatrix_c")]]
NumericMatrix createMatrix_c(IntegerMatrix dataDistribution, NumericVector dataCost){
  int numberRows    = dataDistribution.nrow();
  int numberColumns = dataDistribution.ncol();
  NumericMatrix matrix_c(numberRows, numberColumns);
  
  for(int i = 0; i < numberColumns; i++){
    for(int j = 0; j < numberRows; j++){
      if(dataDistribution(j, i) != 0){
        matrix_c(j, i) = dataCost[i];
      } 
    }
  }
  return matrix_c;
}


/*** R
#Write code in R language.
***/