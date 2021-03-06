#ifndef GUARD_Package_h
#define GUARD_Package_h


/* Load header files, set plugins, load Rcpp namespace */

#include <Rcpp.h>            //Librería "Rcpp.h" para conectar c?digo C++ con c?digo R.
// [[Rcpp::depends(BH)]]

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

#include <algorithm>        //Librería "algorithm" que provee un largo num. de algoritmos que trabajan con "iterators".
#include <array>
#include <iostream>
#include <map>
#include <math.h>           //Librería "math.h" es la Librería matem?tica y trigonom?trica de C++.
#include <numeric>          //Librería "numeric" que permite usar funciones como: accumulate(); adjacent_difference(); inner_product(); partial_sum().    
#include <string>
#include <unordered_set>    //Librería para usar "sets" desordenados.
#include <utility>          //Librería "utility" para poder usar la plantilla pair.
#include <vector>
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>


#endif
