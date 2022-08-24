
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double unorm(double x){
  double xabs;
  double cumnorm;
  double expon;
  double build;
  xabs = std::abs(x);
  if (xabs > 37) {
    cumnorm = 0;
  } else {
    expon = exp(- exp(2 * log(xabs)) / 2);
  }
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
  } else {
    build = xabs + 0.65;
    build = xabs + 4 / build;
    build = xabs + 3 / build;
    build = xabs + 2 / build;
    build = xabs + 1 / build;
    cumnorm = expon / build / 2.506628274631;
  }
  if (x > 0){
    cumnorm = 1 - cumnorm;
  }
  return cumnorm;
}

// [[Rcpp::export]]
double bnorm(double h1, double h2, double r){
  double x[] = { 0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992 };
  double w[] = { 0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042 };
  double lh, h12, h3, h5, h6, h7, h8, aa, ab, r1, r2, r3, bcum, rr;
  h12 = (exp(2 * log(std::abs(h1))) + exp(2 * log(std::abs(h2)))) / 2;
  lh = 0.;
  if (std::abs(r) > 0.7){
    r2 = 1 - exp(2 * log(std::abs(r)));
    r3 = sqrt(r2);
    if (r < 0) h2 = - h2;
    h3 = h1 * h2;
    h7 = exp(- h3 / 2);
    if (std::abs(r) < 1.0){
      h6 = std::abs(h1 - h2);
      h5 = exp(2 * log(std::abs(h6))) / 2;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8;
      ab = 3 - 2 * aa * h5;
      lh = 0.13298076 * h6 * ab * (1 - unorm(h6)) -
	exp(- h5 / r2) * (ab + aa * r2) * 0.053051647;
      for (int i = 0; i <= 4; i++){
	r1 = r3 * x[i];
	rr = exp(2 * log(std::abs(r1)));
	r2 = sqrt(1 - rr);
	if (h7 == 0){
	  h8 = 0;
	} else {
	  h8 = exp(- h3 / (1 + r2)) / r2 / h7;
	}
      lh = lh - w[i] * exp(- h5 / rr) * (h8 - 1.0 - aa * rr);
      }
    }
    bcum = lh * r3 * h7 + unorm(std::min(h1, h2));
    if (r < 0){
      bcum = unorm(h1) - bcum;
    }
  } else {
    h3 = h1 * h2;
    if (r != 0){
      for (int i = 0; i <= 4; i++){
	r1 = r * x[i];
	r2 = 1 - exp(2 * log(std::abs(r1)));
	lh = lh + w[i] * exp((r1 * h3 - h12) / r2) / sqrt(r2);
      }
    }
    bcum = unorm(h1) * unorm(h2) + r * lh;
  }
  return bcum;
}

// [[Rcpp::export]]
NumericVector punorm(NumericVector x){
  int n = x.size();
  NumericVector probs(n);
  for (int i = 0; i < n; i++){
    probs[i] = unorm(x[i]);
  }
  return probs;
}

//' Bivariate normal probability
//'
//' @param z1,z2 two vectors of normal variates
//' @param rho the coefficient of correlation
//' @return the probability
// [[Rcpp::export]]
NumericVector pbnorm(NumericVector z1, NumericVector z2, NumericVector rho){
  int n = z1.size();
  NumericVector probs(n);
  for (int i = 0; i < n; i++){
    probs[i] = bnorm(z1[i], z2[i], rho[i]);
  }
  return probs;
}
