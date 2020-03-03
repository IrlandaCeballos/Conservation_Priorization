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
  MAMPData(DataFrame, DataFrame, DataFrame, DataFrame, DataFrame, DataFrame, DataFrame); //Constructor
  DataFrame data01; //Equivalent to target_Data
  DataFrame data02; //Equivalent to unitCost_Data
  DataFrame data03; //Equivalent to boundary_Data
  DataFrame data04; //Equivalent to speciesDistribution_Data
  //New input data for the MAMP problem
  DataFrame data05; //Equivalent to threatsDistribution_Data
  DataFrame data06; //Equivalent to sensibility_Data
  DataFrame data07; //Equivalent to actionCost_Data
  
  void print();
  int getB();
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
  
private:
  int units   = getUnits();
  int species = getSpecies();
  int threats = getThreats();
  
};


#endif