// We can now use the BH package
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/common_factor.hpp>  

using namespace Rcpp;

// [[Rcpp::export]]
int computeGCD(int a, int b) {
  return boost::math::gcd(a, b);
}

// [[Rcpp::export]]
int computeLCM(int a, int b) {
  return boost::math::lcm(a, b);
}

// [[Rcpp::export]]
IntegerVector rcpp_match( CharacterVector x, CharacterVector table){
  return match( x, table ) ;
}

// [[Rcpp::export]]
IntegerVector rcpp_match_2( StringVector x, StringVector table){
  return match( x, table ) ;
}

// [[Rcpp::export]]
IntegerVector rcpp_hola(){
  Rcpp::CharacterVector vector_01 = Rcpp::CharacterVector::create("A", "B", "C", "D");
  Rcpp::StringVector vector_02    = Rcpp::StringVector::create("A", "B", "C", "D");
  
  Rcpp::CharacterVector vector_03 = Rcpp::CharacterVector::create("B", "C");
  Rcpp::StringVector vector_04 = Rcpp::StringVector::create("B", "C");
  
  
  std::string find_string = "C";
  return rcpp_match(vector_04, vector_02 );
}


