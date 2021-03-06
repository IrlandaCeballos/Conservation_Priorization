#ifndef GUARD_MAMPDataClass_h
#define GUARD_MAMPDataClass_h
#include "Package.h"

/* 
 * *************************************************************
 * ********************** CLASS: MAMPData **********************
 * *************************************************************
 */

//CLASS DECLARATION.
class MAMPData{
public:
  //MAMPData();
  MAMPData(DataFrame, DataFrame, DataFrame, DataFrame, DataFrame, DataFrame, List); //Constructor
  DataFrame data01; //Equivalent to target_Data
  DataFrame data02; //Equivalent to unitCost_Data
  DataFrame data03; //Equivalent to boundary_Data
  DataFrame data04; //Equivalent to speciesDistribution_Data
  //New input data for the MAMP problem
  DataFrame data05; //Equivalent to threatsDistribution_Data
  DataFrame data06; //Equivalent to sensibility_Data
  List data07;      //Equivalent to settings_Data 

  
  void print();
  double getBeta1();
  int getUnits();
  int getSpecies();
  int getThreats();
  std::map<int, double> getTarget(); 
  std::map<int, double> getUnitCost(); 
  std::map<int,std::map<int,int>>  getActionCost(); 
  std::map<int,std::map<int,bool>> getSpeciesDistribution();
  std::map<int,std::map<int,bool>> getSpeciesDistribution_t();
  std::map<int,std::map<int,bool>> getThreatsDistribution();
  std::map<int,std::map<int,bool>> getSensibility();
  std::map<int,std::map<int,double>> getBoundary();
  List getSet(String setName);
  //New methods for linearization of the measure of the "local benefit of the species" (constraint MAMP.2)
  double getExponent();
  int getBreakpoints();
  int getSegments();
  NumericVector get_bp();
  NumericVector get_bp3();
  NumericVector get_slope();
  //New methods for the MAMP-E model
  double getBeta2();
  //New methods for "status column" (in "unitCost_Data" and "threatsDistribution_Data")
  List get_UnitStatus();
  List get_ActionStatus(); 

  
private:
  int units    = getUnits();
  int species  = getSpecies();
  int threats  = getThreats();
  double beta1 = getBeta1();
  //New input data for linearization of the measure of the "local benefit of the species" (constraint MAMP.2)
  double exponent = getExponent();    //Equivalent to exponent
  int breakpoints = getBreakpoints(); //Equivalent to numberBreakpoints
  int segments    = getSegments();    //NEW! Equivalent to numberSegments
  //New input data for the MAMP-E model
  double beta2    = getBeta2();
    
};


#endif
