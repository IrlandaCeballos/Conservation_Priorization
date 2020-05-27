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
MAMPData::MAMPData(DataFrame _data01, DataFrame _data02, DataFrame _data03, DataFrame _data04, DataFrame _data05, DataFrame _data06, List _data07) : data01(_data01), data02(_data02), data03(_data03), data04(_data04), data05(_data05), data06(_data06), data07(_data07) { }

//METHOD: It gives you the message "Hello world!"
void MAMPData::print(){
  Rcpp::Rcout << "Hello world!" << std::endl; //You can use only this: "Rcout".    
  //printf("Hello world!\n");
}

//METHOD: It gives you an 'INT Number' from: parameter "beta1".
double MAMPData::getBeta1(){
  //'beta1' is a boundary length modiﬁer: it can be varied for more or less connected reserve systems. 
  double beta1 = data07(0); //'beta1' = 1 by default!
  return beta1;
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
  
  IntegerVector vectorAux05 = data05[0]; //To identify the planning units.
  IntegerVector vectorAux06 = data05[1]; //To identify the threats.
  IntegerVector vectorAux07 = data05[3]; //To identify the action cost.
  int dataSize = data05.nrows();
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

//New methods for linearization of the measure of the "local benefit of the species" (constraint MAMP.2)
//METHOD:
double MAMPData::getExponent(){
  double exponent = data07(2);
  return exponent;
}

//METHOD:
int MAMPData::getBreakpoints(){
  int breakpoints = data07(3);
  return breakpoints;
}

//METHOD:
int MAMPData::getSegments(){
  int segments = breakpoints - 1;
  return segments;
}

//METHOD:
NumericVector MAMPData::get_bp(){
  //"bp" is a numerical vector of cardinality equal to the number of breakpoints.
  NumericVector bp(breakpoints); 
  double valueX;
  
  for(int m = 1; m <= breakpoints; m++){
    valueX    = (1.0/segments)*(m - 1);
    bp(m - 1) = valueX;
  }
  return bp;     
} 

//METHOD:
NumericVector MAMPData::get_bp3(){
  //"bp3" is a numerical vector of cardinality equal to the number of breakpoints.
  NumericVector bp3(breakpoints); 
  NumericVector bp = get_bp();
  double valueY;
  
  for(int m = 1; m <= breakpoints; m++){
    valueY     = pow( bp(m - 1) , exponent );
    bp3(m - 1) = valueY;
  }
  return bp3;
}

//METHOD:
NumericVector MAMPData::get_slope(){
  //"slope" is a numerical vector of cardinality equal to the number of segments.
  NumericVector slope(segments);
  NumericVector bp  = get_bp();
  NumericVector bp3 = get_bp3();
  double deltaY = 0.0;
  double deltaX = 0.0;
  
  for(int m = 0; m < segments; m++){
    deltaY = bp3(m + 1) - bp3(m);
    deltaX = bp(m + 1) - bp(m); 
    slope(m) = deltaY/deltaX;
  }
  return slope;
}


//METHOD: It gives you an 'INT Number' from: parameter 'beta2'.
double MAMPData::getBeta2(){
  //'beta2' is a penalty factor associated to the spatial fragmentation of actions, which has the same goal than "beta1" in MAMP model.
  double beta2 = data07(1); //'beta2' = 1 by default!
  return beta2;
}


//METHOD: 
List MAMPData::get_UnitStatus(){
  IntegerVector unitsUnrestricted; //Unrestricted planning units.
  IntegerVector unitsLockedIn;     //Pre-included planning units.
  IntegerVector unitsLockedOut;    //Pre-excluded planning units.
  IntegerVector vectorAux01 = data02[0]; //To identify the planning unit i.
  IntegerVector vectorAux02 = data02[2]; //To identify the "status" of planning unit i (variable W[i]).
  int dataSize = data02.nrows();
  
  for(int r = 0; r < dataSize; r++){  
    int filterCondition = vectorAux02[r];
    if(filterCondition == 0){
      int unitIndex = vectorAux01[r];
      unitsUnrestricted.push_back(unitIndex);
    }
    if(filterCondition == 2){ //If the unit status equals 2, it means the unit is "locked-in"
      int unitIndex = vectorAux01[r];
      unitsLockedIn.push_back(unitIndex);
    }
    if(filterCondition == 3){ //If the unit status equals 3, it means the unit is "locked-out"
      int unitIndex = vectorAux01[r];
      unitsLockedOut.push_back(unitIndex);
    }
  }
  
  List unitStatus = List::create(Named("Unrestricted") = Rcpp::sort_unique(unitsUnrestricted), _["LockedIn"] = Rcpp::sort_unique(unitsLockedIn),_["LockedOut"] = Rcpp::sort_unique(unitsLockedOut));
  return unitStatus;
}





//METHOD: 
List MAMPData::get_ActionStatus(){
  IntegerVector vectorAux01 = data05[0]; //To identify the planning unit i.
  IntegerVector vectorAux02 = data05[1]; //To identify the threat k.
  IntegerVector vectorAux03 = data05[4]; //To identify the "status" of a possible action (variable X[i,k]).
  int dataSize = data05.nrows();

  IntegerVector unrestrictedIndexI;
  IntegerVector unrestrictedIndexK;
  
  IntegerVector lockedInIndexI;
  IntegerVector lockedInIndexK;
  
  IntegerVector lockedOutIndexI;
  IntegerVector lockedOutIndexK;
  
  for(int r = 0; r < dataSize; r++){
    int filterCondition = vectorAux03[r];
    if(filterCondition == 0){
      unrestrictedIndexI.push_back(vectorAux01[r]);
      unrestrictedIndexK.push_back(vectorAux02[r]);
    }
    if(filterCondition == 2){ //If the action status equals 2, it means the action is "locked-in"
      lockedInIndexI.push_back(vectorAux01[r]);
      lockedInIndexK.push_back(vectorAux02[r]);
    }
    if(filterCondition == 3){ //If the action status equals 3, it means the action is "locked-out"
      lockedOutIndexI.push_back(vectorAux01[r]);
      lockedOutIndexK.push_back(vectorAux02[r]);
    }
    
  }//END for

  IntegerMatrix actionsUnrestricted = cbind(unrestrictedIndexI, unrestrictedIndexK); //Unrestricted planning units.
  actionsUnrestricted.attr("dim")   = Dimension(actionsUnrestricted.nrow(), 2);
  colnames(actionsUnrestricted)     = CharacterVector::create("i", "k");
  
  IntegerMatrix actionsLockedIn = cbind(lockedInIndexI, lockedInIndexK); //Pre-included planning units.
  actionsLockedIn.attr("dim")   = Dimension(actionsLockedIn.nrow(), 2);
  colnames(actionsLockedIn)     = CharacterVector::create("i", "k");
  
  IntegerMatrix actionsLockedOut = cbind(lockedOutIndexI, lockedOutIndexK); //Pre-excluded planning units.
  actionsLockedOut.attr("dim")   = Dimension(actionsLockedOut.nrow(), 2);
  colnames(actionsLockedOut)     = CharacterVector::create("i", "k");

  
  List actionStatus = List::create(Named("Unrestricted") = actionsUnrestricted, _["LockedIn"] = actionsLockedIn,_["LockedOut"] = actionsLockedOut);
  return actionStatus;
}

//RCPP 'EXPOSURE BLOCK'.
RCPP_MODULE(MAMPDatamodule){
  Rcpp::class_<MAMPData>( "MAMPData" )
  //.constructor("documentation for default constructor")
  .constructor<DataFrame,DataFrame,DataFrame,DataFrame,DataFrame,DataFrame,List>("documentation for constructor")
  .field( "data01", &MAMPData::data01, "documentation for data01")
  .field( "data02", &MAMPData::data02, "documentation for data02")
  .field( "data03", &MAMPData::data03, "documentation for data03")
  .field( "data04", &MAMPData::data04, "documentation for data04")
  .field( "data05", &MAMPData::data05, "documentation for data05")
  .field( "data06", &MAMPData::data06, "documentation for data06")
  .field( "data07", &MAMPData::data07, "documentation for data07")
  //New input data for linearization of the measure of the "local benefit of the species" (constraint MAMP.2)
  //.field( "exponent", &MAMPData::exponent   , "documentation for exponent")
  //.field( "breakpoints", &MAMPData::breakpoints, "documentation for breakpoints")

  .method( "print", &MAMPData::print,  "documentation for print")
  .method( "getBeta1",  &MAMPData::getBeta1,   "documentation for getBeta1")
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
  .method( "getSet",         &MAMPData::getSet,         "documentation for getSet")
  //New methods for linearization of the measure of the "local benefit of the species" (constraint MAMP.2)
  .method( "getExponent",    &  MAMPData::getExponent,  "documentation for getExponent")
  .method( "getBreakpoints", &MAMPData::getBreakpoints, "documentation for getBreakpoints")
  .method( "getSegments",    &MAMPData::getSegments,    "documentation for getSegments")
  .method( "get_bp",         &MAMPData::get_bp,         "documentation for get_bp") 
  .method( "get_bp3",        &MAMPData::get_bp3,        "documentation for get_bp3") 
  .method( "get_slope",      &MAMPData::get_slope,      "documentation for get_slope")   
  .method( "getBeta2",       &MAMPData::getBeta2,       "documentation for getBeta2") 
  //New methods for "status column" (in "unitCost_Data" and "threatsDistribution_Data") 
  .method( "get_UnitStatus", &MAMPData::get_UnitStatus, "documentation for get_UnitStatus")
  .method( "get_ActionStatus", &MAMPData::get_ActionStatus, "documentation for get_ActionStatus")
  
  ; //ATENTION! With ";" 
}

