#include<Rcpp.h>
using namespace Rcpp;

class Location{
public:
  Location();
  //Location(int, int, int, int, int, int, int, int, int, int);
  Location(int, int, int, int, int, int, int);
  int a;
  int b;
  int c;
  int d;
  int e;
  int f;
  int g;
  //int h;
  // int i;
  // int j;
  void print();
  
private:

};

//constructors
// Location::Location() :a(0), b(0), c(0), d(0), e(0), f(0), g(0), h(0), i(0), j(0) { }
// Location::Location(int ai, int bi, int ci, int di, int ei, int fi, int gi, int hi, int ii, int ji) :a(ai), b(bi), c(ci), d(di), e(ei), f(fi), g(gi), h(hi), i(ii), j(ji){ }

Location::Location() :a(0), b(0), c(0), d(0), e(0), f(0), g(0) { }
Location::Location(int ai, int bi, int ci, int di, int ei, int fi, int gi) :a(ai), b(bi), c(ci), d(di), e(ei), f(fi), g(gi){ }


//print function
void Location::print(){
  Rcpp::Rcout << "a = " << a << std::endl;
  Rcpp::Rcout << "b = " << b << std::endl;
}

RCPP_MODULE(locationmodule){
  Rcpp::class_<Location>( "Location" )
  .constructor("documentation for default constructor")
  //.constructor<int, int, int, int, int, int, int, int, int, int>("documentation for constructor")
  .constructor<int, int, int, int, int, int, int>("documentation for constructor")
  .field_readonly( "a", &Location::a, "documentation for a")
  .field( "b", &Location::b, "documentation for b")
  .field( "c", &Location::c, "documentation for c")
  .field( "d", &Location::d, "documentation for d")
  .field( "e", &Location::e, "documentation for e")
  .field( "f", &Location::f, "documentation for f")
  .field( "g", &Location::g, "documentation for g")
  //.field( "h", &Location::h, "documentation for h")
  // .field( "i", &Location::i, "documentation for i")
  // .field( "j", &Location::j, "documentation for j")
  .method( "print", &Location::print, "documentation for print")
  ;
}





