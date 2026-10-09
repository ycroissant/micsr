#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double punorm0(double x) {
  double xabs, cumnorm, expon, build;
  xabs = sqrt(pow(x, 2));
  if (xabs > 37){
    cumnorm = 0;
  } else {
    expon = exp(- pow(xabs, 2) / 2);
    if (xabs < 7.07106781186547) {
      build = 3.52624965998911E-2 * xabs + 0.700383064443688;
      build = build * xabs + 6.37396220353165;
      build = build * xabs + 33.912866078383;
      build = build * xabs + 112.079291497871;
      build = build * xabs + 221.213596169931;
      build = build * xabs + 220.206867912376;
      cumnorm = expon * build;
      build = 8.83883476483184E-02 * xabs + 1.75566716318264;
      build = build * xabs + 16.064177579207;
      build = build * xabs + 86.7807322029461;
      build = build * xabs + 296.564248779674;
      build = build * xabs + 637.333633378831;
      build = build * xabs + 793.826512519948;
      build = build * xabs + 440.413735824752;
      cumnorm = cumnorm / build;
    } else{
      build = xabs + 0.65;
      build = xabs + 4 / build;
      build = xabs + 3 / build;
      build = xabs + 2 / build;
      build = xabs + 1 / build;
      cumnorm = expon / build / 2.506628274631;
    }
  }
  if (x > 0){
    cumnorm = 1 - cumnorm;
  }
  return(cumnorm);
}

// [[Rcpp::export]]
double pbnorm0(double h1, double h2, double rho) {
  NumericVector x(5), w(5);
  x = {0.04691008,  0.23076534,  0.5,          0.76923466,  0.95308992};
  w = {0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042};
  double lh, h12, aa, ab, rr, bcum;
  double r1, r2, r3;
  double h3, h5, h6, h7, h8;
  int j;
  lh = 0.0;
  h12 = (pow(h1, 2) + pow(h2, 2)) / 2;
  if (sqrt(pow(rho, 2)) >= 0.7){
    r2 = 1 - pow(rho, 2);
    r3 = sqrt(r2);
    if (rho < 0){
      h2 = - h2;
    }
    h3 = h1 * h2;
    h7 = exp(- h3 / 2);
    if (sqrt(pow(rho, 2)) < 1){
      h6 = sqrt(pow(h1 - h2, 2));
      h5 = pow(h6, 2) / 2;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8;
      ab = 3 - 2 * aa * h5;
      lh = 0.13298076 * h6 * ab * (1 - punorm0(h6)) - exp(- h5 / r2) * (ab + aa * r2) * 0.053051647;
      for (j = 0; j < 5; j++){
	r1 = r3 * x[j];
	rr = pow(r1, 2);
	r2 = sqrt(1 - rr);
	if (h7 == 0) {
	  h8 = 0;
	} else {
	  h8 = exp(- h3 / (1 + r2)) / r2 / h7;
	}
	lh = lh - w[j] * exp(- h5 / rr) * (h8 - 1 - aa * rr);
      }
    }
    if (h1 >= h2){
      bcum = lh * r3 * h7 + punorm0(h2);
    } else {
      bcum = lh * r3 * h7 + punorm0(h1);
    }      
    if (rho < 0){
      bcum = punorm0(h1) - bcum;
    }
  } else {
    h3 = h1 * h2;
    if (rho != 0){
      for (j = 0; j < 5; j++){
	r1 = rho * x[j];
	r2 = 1 - pow(r1, 2);
	lh = lh + w[j] * exp( (r1 * h3 - h12) / r2) / sqrt(r2);
      }
    }
    bcum = punorm0(h1) * punorm0(h2) + rho * lh;
  }
  return(bcum);
}

// [[Rcpp::export]]
double ptnorm0(double z1, double z2, double z3, double rho12, double rho13, double rho23){
  NumericVector x(5), w(5), z(3), rho(3);
  x = {0.04691008,  0.23076534,  0.5,          0.76923466,  0.95308992};
  w = {0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042};
  z = {z1, z2, z3};
  rho = {rho12, rho13, rho23};
  double r12, r13, r23;
  double xi, xis;
  double tnorm, del, fac;
  double hp1, hp2, z12, z13;
  int j;
  if ((pow(rho12, 2) >= pow(rho13, 2)) && (pow(rho12, 2) >= pow(rho23, 2))){
    z1 = z[0]; z2 = z[1]; z3 = z[2];
    r12 = rho[0]; r13 = rho[1]; r23 = rho[2];
  }
  if ((pow(rho13, 2) >= pow(rho12, 2)) && (pow(rho13, 2) >= pow(rho23, 2))){
    z1 = z[1]; z2 = z[0]; z3 = z[3];
    r12 = rho[0]; r13 = rho[2]; r23 = rho[1];
  }
  if ((pow(rho23, 2) >= pow(rho12, 2)) && (pow(rho23, 2) >= pow(rho13, 2))){
    z1 = z[1]; z2 = z[2]; z3 = z[0];
    r12 = rho[2]; r13 = rho[0]; r23 = rho[1];
  }
  tnorm = 0.0;
  for (j = 0; j < 5; j++){
    xi = x[j];
    xis = pow(xi, 2);
    z12 = exp(- 0.5 * (pow(z1, 2) + pow(z2, 2) - 2 * xi * r12 * z1 * z2) /
	      (1.0 - xis * pow(r12, 2))) /
      sqrt(1.0 - xis * pow(r12, 2));
    z13 = exp(- 0.5 * (pow(z1, 2) + pow(z3, 2) - 2 * xi * r13 * z1 * z3) /
	      (1.0 - xis * pow(r13, 2))) /
      sqrt(1.0 - xis * pow(r13, 2));
    del = 1.0 - xis * pow(r12, 2) - xis * pow(r13, 2) - pow(r23, 2) + 2 * xis * r12 * r13 * r23;
    fac = sqrt(del);
    hp1 = (z3 * (1.0 - xis * pow(r12, 2)) - z1 * (xi * r13 - xi * r12 * r23) -
	   z2 * (r23 - xis * r12 * r13)) /
      fac / sqrt(1.0 - xis * pow(r12, 2));
    hp2 = (z2 * (1.0 - xis * pow(r13, 2)) - z1 * (xi * r12 - xi * r13 * r23) -
	   z3 * (r23 - xis * r12 * r13)) /
      fac / sqrt(1.0 - xis * pow(r13, 2));
    tnorm = tnorm + w[j] * z12 * punorm0(hp1) * r12;
    tnorm = tnorm + w[j] * z13 * punorm0(hp2) * r13;
  }
  tnorm = tnorm + punorm0(z1) * pbnorm0(z2, z3, r23);
  return(tnorm);
}

//' Compute the probability for the univariate normal function
//'
//' @param z a numeric vector
//' @return a numeric vector
//' @export
// [[Rcpp::export]]
NumericVector punorm(NumericVector z) {
  int N = z.size();
  NumericVector result(N);
  for (int i = 0; i < N; i++){
    result[i] = punorm0(z[i]);
  }
  return(result);
}

//' Compute the probability for the bivariate normal function
//'
//' @param z1,z2 two numeric vectors
//' @param rho a numeric vector
//' @return a numeric vector
//' @export
// [[Rcpp::export]]
NumericVector pbnorm(NumericVector z1, NumericVector z2, NumericVector rho) {
  int N = z1.size();
  int i;
  double arho;
  NumericVector result(N);
  int Nrho = rho.size();
  for (i = 0; i < N; i++){
    if (Nrho == 1){
      arho = rho[0];
    } else {
      arho = rho[i];
    }
    result[i] = pbnorm0(z1[i], z2[i], arho);
  }
  return(result);
}
    
// [[Rcpp::export]]
NumericVector ptnormv(NumericVector z1, NumericVector z2, NumericVector z3, NumericVector rho12, NumericVector rho13, NumericVector rho23) {
  int N = z1.size();
  int i, j;
  NumericVector result(N);
  int Nrho = rho12.size();
  for (i = 0; i < N; i++){
    if (Nrho == 1){
      j = 0;
    } else {
      j = i;
    }
    result[i] = ptnorm0(z1[i], z2[i], z3[i], rho12[j], rho13[j], rho23[j]);
  }
  return(result);
}

//' Compute the probability for the trivariate normal function
//'
//' @param z a matrix with three columns
//' @param rho a matrix with three columns
//' @return a numeric vector
//' @export
// [[Rcpp::export]]
NumericVector ptnorm(NumericMatrix z, NumericMatrix rho) {
  int N = z.nrow();
  int i, j;
  NumericVector result(N);
  int Nrho = rho.nrow();
  for (i = 0; i < N; i++){
    if (Nrho == 1){
      j = 0;
    } else {
      j = i;
    }
    result[i] = ptnorm0(z(i, 0), z(i, 1), z(i, 2), rho(j, 0), rho(j, 1), rho(j, 2));
  }
  return(result);
}
