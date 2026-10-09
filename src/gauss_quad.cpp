#include <Rcpp.h>
#include <string>

using namespace Rcpp;

//' Gauss-Laguerre quadrature
//'
//' Computes the node and the weights for the Gauss-Laguerre quadrature (integral on the whole real line)
//'
//' @name gauss_laguerre
//' @param N the number of evaluations
//' @return a list containing two numeric vectors of length N, the first one containing the nodes and the second one the weights
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
List gauss_laguerre(int N) {
  double eps = 1E-14;
  double z, z1, p1, p2, p3, pp, ai;
  NumericVector x(N);
  NumericVector w(N);
  for (int i = 0; i < N; i++){
    if (i == 0){
      z = 3 / (1 + 2.4 * N);
    } else if (i == 1){
      z = z + 15 / (1 + 2.5 * N);
    } else {
      ai = i - 1;
      z = z + ( (1 + 2.55 * ai) / (1.9 * ai)) * (z - x[i - 2]);
    }
    z1 = 1E03;
    while ( sqrt(pow(z1 - z, 2)) > eps){
      p1 = 1;
      p2 = 0;
      for (int j = 0; j < N; j++){
	p3 = p2;
	p2 = p1;
	p1 = ( (2 * j + 1 - z) * p2 - j * p3) / (j + 1);
      }
      pp = (N * p1 - N * p2) / z;
      z1 = z;
      z = z1 - p1 / pp;
    }
    x[i] = z;
    w[i] = - 1 / (pp * N * p2);
  }
  return(List::create(_["nodes"] = x, _["weights"] = w));
}

//' Gauss-Hermitte quadrature
//'
//' Computes the node and the weights for the Gauss-Hermite quadrature (integral on the whole real line)
//'
//' @name gaussian_quad
//' @param N the number of evaluations
//' @return a list containing two numeric vectors of length N, the first one containing the nodes and the second one the weights
//' @export
// [[Rcpp::export]]
List gauss_hermite(int N) {
  double eps = 1E-14;
  double z, z1, p1, p2, p3, pp;
  NumericVector x(N);
  NumericVector w(N);
  int m = (N + 1) / 2;
  for (int i = 0; i < m; i++){
    if (i == 0){
      z = sqrt(2 * N + 1) - 1.85575 * pow(2 * N + 1, - 0.1667);
    } else if (i == 1){
      z -= 1.14 * pow(N, 0.426) / z;
    } else if (i == 2) {
      z =  1.86 * z - 0.86 * x[0];
    } else if (i == 3) {
      z = 1.91 * z - 0.91 * x[1];
    } else {
      z = 2 * z - x[i - 2];
    }
    z1 = 1E03;
    while ( sqrt(pow(z1 - z, 2)) > eps){
      p1 = pow(3.141593, - 0.25);
      p2 = 0;
      for (int j = 0; j < N; j++){
	p3 = p2;
	p2 = p1;
	p1 = z * sqrt(2 / (1.0 * j + 1)) * p2 - sqrt(1.0 * j / (1.0 * j + 1)) * p3;
      }
      pp = sqrt(2 * N) * p2;
      z1 = z;
      z = z1 - p1 / pp;
    }
    x[i] = z;
    x[N - 1 - i] = - z;
    w[i] = 2 / pow(pp, 2);
    w[N - 1 - i] = w[i];
  }
 return(List::create(_["nodes"] = x, _["weights"] = w));
}


// //' Gaussian quadrature
// //'
// //' Computes the node and the weights for the Gauss-Hermite quadrature (integral on the whole real line)
// //'
// //' @name gaussian_quad
// //' @param N the number of evaluations
// //' @param type one of "hermite" or "laguerre"
// //' @return a list containing two numeric vectors of length N, the first one containing the nodes and the second one the weights
// //' @export
// // [[Rcpp::export]]
// List gaussian_quad(int N, std::string kind){
//   if (kind == "hermite"){
//     return(gauss_hermite(N));
//   } else {
//     if (kind == "laguerre"){
//       return(gauss_laguerre(N));
//     }
//   }
// }
