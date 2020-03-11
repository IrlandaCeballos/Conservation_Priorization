#include <Rcpp.h>
#include <iostream>
#include <math.h>        //Librería "math.h" es la librería matemática y trigonométrica de C++
#include <numeric>       //Librería "numeric" que permite usar funciones como: accumulate(); adjacent_difference(); inner_product(); partial_sum().-
#include <unordered_set> //Librería para usar "sets" desordenados
#include <vector>
#include <algorithm>     //Librería "algorithm" que provee un largo num. de algoritmos que trabajan con "iterators".-
#include <string>
#include <array>


// #include <utility> //Librería "utility" para poder usar la plantilla pair


// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// [[Rcpp::export]]
void test_copy(){
  NumericVector A = NumericVector::create(1, 2, 3);
  NumericVector B = A;
  Rcout << "Before: " << std::endl << "A: " << A << std::endl << "B: " << B << std::endl; 
  A[1] = 5; // 2 -> 5
  Rcout << "After: " << std::endl << "A: " << A << std::endl << "B: " << B << std::endl; 
}

// [[Rcpp::export]]
Rcpp::NumericVector foo(const Rcpp::NumericVector& x) {
  Rcpp::NumericVector tmp(x.length());
  for (int i = 0; i < x.length(); i++)
    tmp[i] = x[i] + 1.0;
  return tmp;
}

// [[Rcpp::export]]
std::vector<double> bar(const std::vector<double>& x) {
  std::vector<double> tmp(x.size());
  for (int i = 0; i < x.size(); i++)
    tmp[i] = x[i] + 1.0;
  return tmp;
}

// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix x){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  for(int i=0;i<nrow;i++){
    double total=0;
    for(int j=0;j<ncol;j++){
      total+=x(i,j);
    }
    out[i]=total;
  }
  return out;
}

// Función que obtiene la distancia euclidiana entre un escalar "x" y un escalar ys[i] del vector ys.
// [[Rcpp::export]]
NumericVector pdistC(double x, NumericVector ys){
  int n = ys.size();
  NumericVector out(n);
  for(int i=0; i<n; ++i){
    out[i]=sqrt(pow(ys[i]-x, 2.0));
  }
  return out;
}

// Función que obtiene la distancia euclidiana entre un escalar "x" y un escalar ys[i] del vector ys.
// [[Rcpp::export]]
NumericVector pdistC2(double x, NumericVector ys){
  return sqrt(pow((x-ys), 2)); //Con funciones de "sugar"
}

// Equivale a la función mean() en R.
// [[Rcpp::export]]
double meanC(NumericVector x){
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i){
    total += x[i];
    }
  return total/n;
}

// Equivale a la función all() en R.
// [[Rcpp::export]]
double f1(NumericVector x){
  int n = x.size();
  double y = 0;
  
  for(int i = 0; i<n; ++i){
    y += x[i]/n;  
  }
  return y;
}

// [[Rcpp::export]]
NumericVector f2(NumericVector x){
  int n = x.size();
  NumericVector out(n);
  
  out[0] = x[0];
  for(int i = 1; i<n; ++i){
    out[i] = out[i - 1] + x[i];
    }
  return out;
}

// [[Rcpp::export]]
bool f3(LogicalVector x){
  int n = x.size();
  
  for(int i = 0; i < n; ++i){
    if(x[i]) return true;
  }
  return false;
}

// [[Rcpp::export]]
NumericVector f5(NumericVector x, NumericVector y){
  int n = std::max(x.size(), y.size());
  NumericVector x1 = rep_len(x, n);
  NumericVector y1 = rep_len(y, n);
  
  NumericVector out(n);
  for(int i = 0; i < n; ++i){
    out[i] = std::min(x1[i], y1[i]);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector attribs(){
  NumericVector out = NumericVector::create(1,2,3);
  out.names() = CharacterVector::create("a","b","c");
  out.attr("my-attr") = "my-value";
  out.attr("class") = "my-class";
  
  return out;
  }

// [[Rcpp::export]]
double mpe(List mod){
  if(!mod.inherits("lm")) stop("Input must be a linear model");
  
  NumericVector resid  = as<NumericVector>(mod["residuals"]);  
  NumericVector fitted = as<NumericVector>(mod["fitted.values"]);
  
  int n = resid.size();
  double err = 0;
  for(int i = 0; i<n; ++i){
    err += resid[i] / (fitted[i] + resid[i]);
  }
  return err/n;
}

// [[Rcpp::export]]
RObject callWithOne(Function f){
  return f(1);
}

//Equivalente a la función lapply de R
// [[Rcpp::export]]
List lapply1(List input, Function f){
  int n = input.size();
  List out(n);
  
  for(int i = 0; i < n; i++){
    out[i] = f(input[i]);
  }
  return out;
}

//Checkea si un valor en un vector es NA
// [[Rcpp::export]]
LogicalVector is_naC(NumericVector x){
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i){
    out[i] = NumericVector::is_na(x[i]);
  }
  return out;
}

//Checkea si un valor en un vector es NA
// [[Rcpp::export]]
LogicalVector is_naC2(NumericVector x){
  return is_na(x); //usando el la función de "sugar"
}  

// Suma los elementos de un vector numérico X dado.
// [[Rcpp::export]]  
double sum3(NumericVector x){
  double total = 0;
  
  NumericVector::iterator it;
  for(it = x.begin(); it != x.end(); ++it){
    total += *it;
  }
  return total;
}

// Suma los elementos de un vector numérico X dado.
// [[Rcpp::export]]  
double sum4(NumericVector x){
  return std::accumulate(x.begin(), x.end(), 0.0);
}

//Versión en C++ de la función findInterval()  
// [[Rcpp::export]]   
IntegerVector findInterval2(NumericVector x, NumericVector breaks){
  IntegerVector out(x.size());
  
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  for(it = x.begin(), out_it = out.begin(); it!= x.end();++it, ++out_it){
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }
  return out;
}

//Versión en C++ de una "tabla de frecuencias", pero solo contabiliza la frecuencia de un número esn sus posiciones adyacentes 
// [[Rcpp::export]]   
List rleC(NumericVector x){
  std::vector<int> lengths;
  std::vector<double> values;
  
  //Initialise first value
  int i = 0;
  double prev = x[0];
  
  values.push_back(prev);
  lengths.push_back(1);
  
  NumericVector::iterator it;
  for(it = x.begin() + 1; it!=x.end(); ++it){
    if(prev == *it){
      lengths[i]++;
    }else{
      values.push_back(*it);
      lengths.push_back(1);
      i++;
      prev = *it;
    }
  }
  return List::create(
    _["Lengths"] = lengths,
    _["values"]  = values
    );
}

//Versión en C++ de duplicated()(en R) para vectores de enteros.-
// [[Rcpp::export]]  
LogicalVector duplicatedC(IntegerVector x){
  std::unordered_set<int> seen;
  int n = x.size();
  LogicalVector out(n);
  
  for(int i = 0; i<n; ++i){
    out[i] = !seen.insert(x[i]).second; 
    //mucho razonamiento aquí :s, insert() devuelve dos parámetros: 
    //.first y .second (este es TRUE si el valor es una nueva adición al set)
  return out;
  }
}

//Versión en C++ de una tabla de frecuencias, dado un vector numérico x.-
// [[Rcpp::export]]  
std::map<double, int> tableC (NumericVector x){
  std::map<double, int> counts;
  int n = x.size();
    
  for(int i = 0; i < n; i++){
    counts[x[i]]++;
  } 
  return counts;
}

  
  
/*** R
#Escribir codigo en lenguaje R.-

#library(microbenchmark)
#w <- runif(1e5)
#microbenchmark(
#  mean(w),
#  meanC(w)
#)

***/

