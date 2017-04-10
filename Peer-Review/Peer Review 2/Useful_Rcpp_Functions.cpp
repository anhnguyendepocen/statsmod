// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

/*-------------------------------------------------------
#		Note about sampling from gammas:
#		R::rgamma(sh,rt) seems to work fine.  It is not vectorized; to sample multiples, use a for loop.  Still works quickly.
#-------------------------------------------------------*/	

/*-------------------------------------------------------
#	Replicates R's seq(from,to,by) function.
#-------------------------------------------------------*/	
	struct add_multiple {
	  int incr;
	  int count;
	  add_multiple(int incr)
	    : incr(incr), count(0)
	    {}
	  inline int operator()(int d) {
	    return d + incr * count++;
	  }
	};
	
	// [[Rcpp::export]]
	Rcpp::NumericVector rcpp_seq(double from_, double to_, double by_ = 1.0) {
	  int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
	  int from = adjust * from_;
	  int to = adjust * to_;
	  int by = adjust * by_;

	  std::size_t n = ((to - from) / by) + 1;
	  Rcpp::IntegerVector res = Rcpp::rep(from, n);
	  add_multiple ftor(by);

	  std::transform(res.begin(), res.end(), res.begin(), ftor);
	  return Rcpp::NumericVector(res) / adjust;
	}
		

// [[Rcpp::export]]
/*-------------------------------------------------------
# Generate draws from a multivariate normal distribution
#--------------------------------------------------------
#  n        number of samples
#  mu       mean vector
#  Sigma    covariance matrix
#-------------------------------------------------------*/
mat mvrnorm(int n, vec mu, mat sigma) {
	   int ncols = sigma.n_cols;
	   mat Y = randn(n, ncols);
		 mat result = (repmat(mu, 1, n).t() + Y * chol(sigma)).t();
		 return result;
}

/*-------------------------------------------------------
# Generate Draws from an Inverse Wishart Distribution
# via the Bartlett Decomposition. Replicates riwish from MCMCPack.
#--------------------------------------------------------
#   n     number of samples
#   S     scale matrix 
#   v     degrees of freedom
#-------------------------------------------------------*/
// [[Rcpp::export]]
cube rinvwish(int n, int v, mat S){
  RNGScope scope;
  int p = S.n_rows;
  mat L = chol(inv_sympd(S), "lower");
  cube sims(p, p, n, fill::zeros);
  for(int j = 0; j < n; j++){
  mat A(p,p, fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df)); 
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
  sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}