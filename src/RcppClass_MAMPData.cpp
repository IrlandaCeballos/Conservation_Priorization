#include "MAMPData.h"

/* 
 * *************************************************************
 * ********************** CLASS: MAMPData **********************
 * *************************************************************
 */

//CLASS DEFINITION.
//CONSTRUCTORS.
//MAMPData::MAMPData() :x(0), y(0) { } //To set a default constructor!
//(NOTE: it's advisable to change the name of data01 -> _data01)
MAMPData::MAMPData(DataFrame _data01, DataFrame _data02, DataFrame _data03, DataFrame _data04, DataFrame _data05, DataFrame _data06, DataFrame _data07) : data01(_data01), data02(_data02), data03(_data03), data04(_data04), data05(_data05), data06(_data06), data07(_data07) { }


//METHOD: It gives you the message "Hello world!"
void MAMPData::print(){
  Rcpp::Rcout << "Hello world!" << std::endl; //You can use only this: "Rcout".    
  //printf("Hello world!\n");
  //std::cout << "Hello world!" << std::endl;  
}

//METHOD: It gives you an 'INT Number' from: parameter B.
int MAMPData::getB(){
  //'Parameter B' is a boundary length modiﬁer: it can be varied for more or less connected reserve systems. 
  //'Parameter B' = 1 by default!
  int parameterB = 1; 
  return parameterB;
}

//METHOD: It gets you an 'INT Number' from: Planning Units (returns the number of planning units).
int MAMPData::getUnits(){
  int units = data02.nrows();
  return units;
}

//METHOD: It gets you an 'INT Number' from: Species (returns the number of species).
int MAMPData::getSpecies(){
  int species = data01.nrows();
  return species;
}

//METHOD: It gets you an 'INT Number' from: Threats (returns the number of threats).
int MAMPData::getThreats(){
  IntegerVector vectorAux = data05[1];
  int threats = *std::max_element(vectorAux.begin(), vectorAux.end()); //t
  //I've used asterisk (*) to convert an 'iterator' type (default return from max_element() function) in an 'int' type.
  return threats;
}

//METHOD: It gets you a 'MAP<int, double>' from: Species|Target.
std::map<int, double> MAMPData::getTarget(){
  typedef std::pair<int, double> targetDataAux;
  std::map<int, double> targetData;
  std::map<int, double>::iterator it = targetData.begin();
  IntegerVector vectorAux01 = data01[0];
  NumericVector vectorAux02 = data01[1];
  
  for(int i = 0; i < species; i++){
    targetData.insert(it, targetDataAux(vectorAux01[i], vectorAux02[i]));
    it ++;
  }
  
  return targetData;     
}

//METHOD: It gets you a 'MAP<int, double>' from: Planning Units|Unit Cost.
std::map<int, double> MAMPData::getUnitCost(){
  typedef std::pair<int, double> dataUnitCostAux;
  std::map<int, double> dataUnitCost;
  std::map<int, double>::iterator it = dataUnitCost.begin();
  IntegerVector vectorAux03 = data02[0];
  NumericVector vectorAux04 = data02[1];

  for(int i = 0; i < units; i++){
    dataUnitCost.insert(it, dataUnitCostAux(vectorAux03[i], vectorAux04[i]));
    it ++;
    //dataUnitCost[i] = vectorAux04[i]; //Simple way to insert data.
  }
  return dataUnitCost;
}

//METHOD: It gets you a 'MAP<int,std::map<int,int>>' from: Planning Units|Threat|Action Cost.
std::map<int,std::map<int,int>> MAMPData::getActionCost(){
  std::map<int,std::map<int,int>> actionCostData; //It stores the cost of applying an action to eliminate a threat k in planning unit i.
  
  IntegerVector vectorAux05 = data07[0]; //To identify the planning units.
  IntegerVector vectorAux06 = data07[1]; //To identify the threats.
  IntegerVector vectorAux07 = data07[2]; //To identify the action cost.
  int dataSize = data07.nrows();
  int iAux     = 0;
  int jAux     = 0;
  int kAux     = 0;
  
  for(int r = 0; r < dataSize; r++){
    iAux = vectorAux05[r]; //Unit in specific.
    jAux = vectorAux06[r]; //Threat in specific.
    kAux = vectorAux07[r]; //Cost of eliminating threat k in unit i.
    actionCostData[iAux][jAux] = kAux;
  }
  return actionCostData;
}

//METHOD: It gets you a 'MAP<int,std::map<int,bool>>' from: Planning Units|Species|Amount.
//Corresponds to the set S[i] (only create S[i] subsets that are NOT empty!).
std::map<int,std::map<int,bool>> MAMPData::getSpeciesDistribution(){
  std::map<int,std::map<int,bool>> speciesDistributionData; //It stores the value of species and their quantity in the respective planning unit.
  
  IntegerVector vectorAux05 = data04[0]; //To identify the planning units.
  IntegerVector vectorAux06 = data04[1]; //To identify the species.
  LogicalVector vectorAux07 = data04[2]; //To identify the quantity of each species.
  int dataSize = data04.nrows();
  int iAux     = 0;
  int jAux     = 0;
  bool kAux    = 0;
  
  for(int r = 0; r < dataSize; r++){
    iAux = vectorAux05[r]; //Unit in specific.
    jAux = vectorAux06[r]; //Specie in specific.
    kAux = vectorAux07[r]; //Is specie j in unit i? k = TRUE if so, k = FALSE if it isn't.
    speciesDistributionData[iAux][jAux] = kAux;
  }
  
  return speciesDistributionData;
}

//METHOD: It gets you a 'MAP<int,std::map<int,bool>>' from: Species|Planning Units|Amount.
//Corresponds to the "transpose" of the species distribution data (i.e. set I[s]), 
//but only create I[s] subsets that are NOT empty!
std::map<int,std::map<int,bool>> MAMPData::getSpeciesDistribution_t(){
  std::map<int,std::map<int,bool>> speciesDistributionData_t; //It stores the value of species and their quantity in the respective planning unit.
  
  IntegerVector vectorAux05 = data04[1]; //To identify the species.
  IntegerVector vectorAux06 = data04[0]; //To identify the planning units.
  LogicalVector vectorAux07 = data04[2]; //To identify the quantity of each species.
  int dataSize = data04.nrows();
  int iAux     = 0;
  int jAux     = 0;
  bool kAux    = 0;
  
  for(int r = 0; r < dataSize; r++){
    iAux = vectorAux05[r]; //Specie in specific.
    jAux = vectorAux06[r]; //Unit in specific.
    kAux = vectorAux07[r]; //Unit i is part of the habitat of species j? k = TRUE if so, k = FALSE if it isn't.
    speciesDistributionData_t[iAux][jAux] = kAux;
  }
  return speciesDistributionData_t;
}

//METHOD: It gets you a 'MAP<int,std::map<int,bool>>' from: Planning Units|Threats|Amount.
//Corresponds to the set K[i] (only create K[i] subsets that are NOT empty!).
std::map<int,std::map<int,bool>> MAMPData::getThreatsDistribution(){
  std::map<int,std::map<int,bool>> threatsDistributionData; //It stores the value of threats and their quantity in the respective planning unit.
  
  IntegerVector vectorAux05 = data05[0]; //To identify the planning units.
  IntegerVector vectorAux06 = data05[1]; //To identify the threats.
  LogicalVector vectorAux07 = data05[2]; //To identify the quantity of each threats.
  int dataSize = data05.nrows();
  int iAux     = 0;
  int jAux     = 0;
  bool kAux    = 0;
  
  for(int r = 0; r < dataSize; r++){
    iAux = vectorAux05[r]; //Unit in specific.
    jAux = vectorAux06[r]; //Threat in specific.
    kAux = vectorAux07[r]; //Is threat j in unit i? k = TRUE if so, k = FALSE if it isn't.
    threatsDistributionData[iAux][jAux] = kAux;
  }
  
  return threatsDistributionData;
}

//METHOD: It gets you a 'MAP<int,std::map<int,bool>>' from: Species|Threats|Sensibility.
//Corresponds to the set K[s] (only create K[s] subsets that are NOT empty!).
std::map<int,std::map<int,bool>> MAMPData::getSensibility(){
  std::map<int,std::map<int,bool>> sensibilityData; //It stores the sensibility of species to threats.
  
  IntegerVector vectorAux05 = data06[0]; //To identify the species.
  IntegerVector vectorAux06 = data06[1]; //To identify the threats.
  LogicalVector vectorAux07 = data06[2]; //To identify the threats that affect the species s.
  int dataSize = data06.nrows();
  int iAux     = 0;
  int jAux     = 0;
  bool kAux    = 0;
  
  for(int r = 0; r < dataSize; r++){
    iAux = vectorAux05[r]; //Specie in specific.
    jAux = vectorAux06[r]; //Threat in specific.
    kAux = vectorAux07[r]; //Does threat j affect the species i? k = TRUE if so, k = FALSE if it isn't.
    sensibilityData[iAux][jAux] = kAux;
  }
  
  return sensibilityData;
}

//METHOD: It gets you a 'MAP<int,std::map<int,double>>' from: Planning Unit 'i'|Planning Unit 'j'|Boundary.
std::map<int,std::map<int,double>> MAMPData::getBoundary(){
  std::map<int,std::map<int,double>> boundaryData;
  
  IntegerVector vectorAux08 = data03[0]; //To identify the planning unit i.
  IntegerVector vectorAux09 = data03[1]; //To identify the planning unit j.
  NumericVector vectorAux10 = data03[2]; //To identify the boundary between both units i and j.
  int dataSize = data03.nrows();
  int iAux     = 0;
  int jAux     = 0;
  double kAux  = 0;
  
  for(int r = 0; r < dataSize; r++){
    iAux = vectorAux08[r]; //Unit 'i' in specific.
    jAux = vectorAux09[r]; //Unit 'j' in specific.
    kAux = vectorAux10[r]; //Boundary between both units 'i' and 'j', in specific!
    boundaryData[iAux][jAux] = kAux;
  }
  return boundaryData;
}

//METHOD: It gets you a LIST with the set that you required, for example, 
//"Ki" gives you the set K[i] (the options are: "Si", "Is", "Ki", "Ks"). 
//The method only needs the name of the set (Rcpp::String).
List MAMPData::getSet(String setName){
  std::map<int,std::map<int,bool>> setData; //To identify the data input that is needed.
  int setCardinality = 0;    //To identify the expected cardinality of the set.
  int subsetCardinality = 0; //To identify the expected cardinality of the subsets.
  List setRequired;
  StringVector namesList;
  bool filterCondition;
  
  if(setName == "Si"){
    setData = getSpeciesDistribution();
    setCardinality = units;
    subsetCardinality = species;
  }
  if(setName == "Is"){
    setData = getSpeciesDistribution_t();
    setCardinality = species; 
    subsetCardinality = units;
  }
  if(setName == "Ki"){
    setData = getThreatsDistribution();
    setCardinality = units;
    subsetCardinality = threats; 
  }
  if(setName == "Ks"){
    setData = getSensibility();   
    setCardinality = species;
    subsetCardinality = threats;
  }
  
  for(int i = 1; i <= setCardinality; i++){
    namesList.push_back(std::to_string(i));
    IntegerVector listValues;
    for(int j = 1; j<= subsetCardinality; j++){
      filterCondition = setData[i][j];
      if(filterCondition == true){
        listValues.push_back(j);
      }
    }
    setRequired.push_back(listValues);
  }//END external for! 
  
  setRequired.names() = namesList;
    
  return setRequired;
}



//RCPP 'EXPOSURE BLOCK'.
RCPP_MODULE(MAMPDatamodule){
  Rcpp::class_<MAMPData>( "MAMPData" )
  //.constructor("documentation for default constructor")
  .constructor<DataFrame,DataFrame,DataFrame,DataFrame,DataFrame,DataFrame,DataFrame>("documentation for constructor")
  .field( "data01", &MAMPData::data01, "documentation for data01")
  .field( "data02", &MAMPData::data02, "documentation for data02")
  .field( "data03", &MAMPData::data03, "documentation for data03")
  .field( "data04", &MAMPData::data04, "documentation for data04")
  .field( "data05", &MAMPData::data05, "documentation for data05")
  .field( "data06", &MAMPData::data06, "documentation for data06")
  .field( "data07", &MAMPData::data07, "documentation for data07")

  .method( "print", &MAMPData::print,  "documentation for print")
  .method( "getB",  &MAMPData::getB,   "documentation for getB")
  .method( "getUnits",   &MAMPData::getUnits,   "documentation for getUnits")
  .method( "getSpecies", &MAMPData::getSpecies, "documentation for getSpecies")
  .method( "getThreats", &MAMPData::getThreats, "documentation for getThreats")
  .method( "getTarget",  &MAMPData::getTarget,  "documentation for getTarget")
  .method( "getUnitCost", &MAMPData::getUnitCost, "documentation for getUnitCost")
  .method( "getActionCost", &MAMPData::getActionCost, "documentation for getActionCost")
  .method( "getSpeciesDistribution", &MAMPData::getSpeciesDistribution, "documentation for getSpeciesDistribution")
  .method( "getSpeciesDistribution_t", &MAMPData::getSpeciesDistribution_t, "documentation for getSpeciesDistribution_t")
  .method( "getThreatsDistribution", &MAMPData::getThreatsDistribution, "documentation for getThreatsDistribution")
  .method( "getSensibility", &MAMPData::getSensibility, "documentation for getSensibility")
  .method( "getBoundary",    &MAMPData::getBoundary,    "documentation for getBoundary")
  .method( "getSet",         &MAMPData::getSet,         "documentation for getSet");
}

