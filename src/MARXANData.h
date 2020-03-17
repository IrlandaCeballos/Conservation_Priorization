#ifndef GUARD_MARXANDataClass_h
#define GUARD_MARXANDataClass_h
#include "Package.h"

/* 
 * *************************************************************
 * ********************* CLASS: MARXANData *********************
 * *************************************************************
 */

//CLASS DECLARATION.
class MARXANData{
public:
  //MARXANData();
  MARXANData(DataFrame, DataFrame, DataFrame, DataFrame); //Constructor
  DataFrame data01; //Equivalent to target_Data
  DataFrame data02; //Equivalent to unitCost_Data
  DataFrame data03; //Equivalent to boundary_Data
  DataFrame data04; //Equivalent to speciesDistribution_Data
  void print();
  int getB();
  int getUnits();
  int getSpecies();
  std::map<int, double> getTarget(); 
  std::map<int, double> getUnitCost(); 
  std::map<int,std::map<int,bool>> getSpeciesDistribution();
  std::map<int,std::map<int,bool>> getSpeciesDistribution_t();
  std::map<int,std::map<int,double>> getBoundary();
  List getSet(String setName);
  
private:
  int units   = getUnits();
  int species = getSpecies();

};


#endif
