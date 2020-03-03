#include <Rcpp.h>
#include <vector>
#include <map>
#include <unordered_set> 
#include <utility>        //"Utility" library to use the "pair" template.
#include <iostream> 
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

/* 
 * <<Function that not belong to the Class "DataInstance">>
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
 * <<Function that not belong to the class "DataInstance">>
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

/*** R
#Write code in R language.
test <- getValueTuple(1, 7, dataDistribution)
print(test)
rm(test)

***/



/* 
 * *************************************************************
 * ******************** CLASS: DataInstance ********************
 * *************************************************************
 */

//CLASS DECLARATION.
class DataInstance{
public:
  //DataInstance();
  DataInstance(DataFrame, DataFrame, DataFrame, DataFrame); //Constructor
  DataFrame data01; //Equivalent to dataTarget
  DataFrame data02; //Equivalent to dataCost
  DataFrame data03; //Equivalent to dataBoundary
  DataFrame data04; //Equivalent to dataDistribution
  void print();
  int getB();
  int getUnits();
  int getSpecies();
  std::map<int, double> getTarget(); 
  std::map<int, double> getCost(); 
  std::map<int,std::map<int,bool>> getDistribution();
  std::map<int,std::map<int,double>> getBoundary();
  
private:
  int units = data02.nrows();
  int species = data01.nrows();
  //int threats;

};

//CLASS DEFINITION.
//CONSTRUCTORS.
//DataInstance::DataInstance() :x(0), y(0) { } //To set a default constructor!
//(NOTE: it's advisable to change the name of data01 -> _data01)
DataInstance::DataInstance(DataFrame _data01, DataFrame _data02, DataFrame _data03, DataFrame _data04) : data01(_data01), data02(_data02), data03(_data03), data04(_data04) { }


//METHOD: It gives you the message "Hello world!"
void DataInstance::print(){
  Rcpp::Rcout << "Hello world!" << std::endl; //You can use only this: "Rcout".    
  //printf("Hello world!\n");
  //std::cout << "Hello world!" << std::endl;  
}

//METHOD: It gives you an 'INT Number' from: parameter B.
int DataInstance::getB(){
  //'Parameter B' is a boundary length modiﬁer: it can be varied for more or less connected reserve systems. 
  //'Parameter B' = 1 by default!
  int parameterB = 1; 
  return parameterB;
}

//METHOD: It gets you an 'INT Number' from: Planning Units (returns the number of planning units).
int DataInstance::getUnits(){
  return units;
}

//METHOD: It gets you an 'INT Number' from: Species (returns the number of species).
int DataInstance::getSpecies(){
  return species;
}

//METHOD: It gets you a 'MAP<int, double>' from: Species|Target.
std::map<int, double> DataInstance::getTarget(){
  typedef std::pair<int, double> dataTargetAux;
  std::map<int, double> dataTarget;
  std::map<int, double>::iterator it = dataTarget.begin();
  IntegerVector vectorAux01 = data01[0];
  NumericVector vectorAux02 = data01[1];
  
  for(int i = 0; i < species; i++){
    dataTarget.insert(it, dataTargetAux(vectorAux01[i], vectorAux02[i]));
    it ++;
  }
  
  return dataTarget;     
}

//METHOD: It gets you a 'MAP<int, double>' from: Planning Units|Cost.
std::map<int, double> DataInstance::getCost(){
  typedef std::pair<int, double> dataCostAux;
  std::map<int, double> dataCost;
  std::map<int, double>::iterator it = dataCost.begin();
  IntegerVector vectorAux03 = data02[0];
  NumericVector vectorAux04 = data02[1];

  for(int i = 0; i < units; i++){
    dataCost.insert(it, dataCostAux(vectorAux03[i], vectorAux04[i]));
    it ++;
    //dataCost[i] = vectorAux04[i]; //Simple way to insert data.
  }
  return dataCost;
}

//METHOD: It gets you a 'MAP<int,std::map<int,bool>>' from: Planning Units|Species|Amount.
std::map<int,std::map<int,bool>> DataInstance::getDistribution(){
  std::map<int,std::map<int,bool>> dataDistribution; //It stores the value of species and their quantity in the respective planning unit.
  
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
    dataDistribution[iAux][jAux] = kAux;
  }
  
  return dataDistribution;
}

//METHOD: It gets you a 'MAP<int,std::map<int,double>>' from: Planning Unit 'i'|Planning Unit 'j'|Boundary.
std::map<int,std::map<int,double>> DataInstance::getBoundary(){
  std::map<int,std::map<int,double>> dataBoundary;
  
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
    dataBoundary[iAux][jAux] = kAux;
  }

  return dataBoundary;
}

//RCPP 'EXPOSURE BLOCK'.
RCPP_MODULE(DataInstancemodule){
  Rcpp::class_<DataInstance>( "DataInstance" )
  //.constructor("documentation for default constructor")
  .constructor<DataFrame,DataFrame,DataFrame,DataFrame>("documentation for constructor")
  .field( "data01", &DataInstance::data01, "documentation for data01")
  .field( "data02", &DataInstance::data02, "documentation for data02")
  .field( "data03", &DataInstance::data03, "documentation for data03")
  .field( "data04", &DataInstance::data04, "documentation for data04")
  .method( "print", &DataInstance::print,  "documentation for print")
  .method( "getB",  &DataInstance::getB,   "documentation for getB")
  .method( "getUnits",   &DataInstance::getUnits,   "documentation for getUnits")
  .method( "getSpecies", &DataInstance::getSpecies, "documentation for getSpecies")
  .method( "getTarget",  &DataInstance::getTarget,  "documentation for getTarget")
  .method( "getCost",    &DataInstance::getCost,    "documentation for getCost")
  .method( "getDistribution", &DataInstance::getDistribution, "documentation for getDistribution")
  .method( "getBoundary",     &DataInstance::getBoundary,     "documentation for getBoundary");
}

