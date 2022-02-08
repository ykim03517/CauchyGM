#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gridG(NumericVector g, NumericVector seqx) {
    int n = g.length();
    int p = seqx.length();
    NumericMatrix X(n,p);
    
    for (int j = 0; j < p; j++) {
        for (int i = 0; i < n; i++) {
            if ( g(i) == 0 ) {
                X(i,j) = 0;
            } else if ( g(i) == 1) {
                X(i,j) = seqx(j);
            } else {
                X(i,j) = 1;
            }
        }
    }
    return X;
    
}

