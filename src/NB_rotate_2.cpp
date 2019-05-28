#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


double dbl_abs(double x){
  if(x>0) return x;
  else return -x;
}

// [[Rcpp::export]]

NumericVector dksum(NumericVector x, NumericVector p, NumericVector y, int n, int m, double h){
  int ord = 1;
  NumericVector output(m);
  double denom;
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  for(int i=0; i<=ord; i++){
    Ly(i,0) = pow(-x[0], i)*y[0];
  }
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      Ly(j,i) = pow(-x[i],j)*y[i] + exp((x[i-1]-x[i])/h)*Ly(j,i-1);
      Ry(j,n-i-1) = exp((x[n-i-1]-x[n-i])/h)*(pow(x[n-i],j)*y[n-i]+Ry(j,n-i));
    }
  }
  int count;
  NumericVector counts(m);
  count = 0;
  for(int i=0; i<m; i++){
    if(p[i]>=x[n-1]){
      for(int j=i; j<m; j++) counts[j] = n;
      break;
    }
    else{
      while(x[count]<=p[i]) count += 1;
      counts[i] = count;
    } 
  }
  for(int orddo=1; orddo<=ord; orddo++){
    NumericVector coefs(orddo + 1);
    coefs[0] = coefs[orddo] = 1;
    if(orddo>1){
      double num = 1;
      for(int j=2; j<=orddo; j++) num *= j;
      double denom1 = 1;
      double denom2 = num/orddo;
      for(int i=2; i<=orddo; i++){
        coefs[i-1] = num/denom1/denom2;
        denom1 *= i;
        denom2 /= (orddo-i+1);
      }
    }
    denom = pow(h, orddo);
    double exp_mult;
    for(int i=0; i<m; i++){
      if(counts[i]==0){
        exp_mult = exp((p[i]-x[0])/h);
        output[i] -= (x[0]-p[i])*y[0]/denom/4.0*exp_mult;
        for(int j=0; j<=1; j++) output[i] += coefs[j]*pow(-p[i],orddo-j)*Ry(j,0)/denom/4.0*exp_mult;
      }
      else{
        exp_mult = exp((x[counts[i]-1]-p[i])/h);
        for(int j=0; j<=1; j++) output[i] -= coefs[j]*(pow(p[i],orddo-j)*Ly(j,counts[i]-1)*exp_mult-pow(-p[i],orddo-j)*Ry(j,counts[i]-1)/exp_mult)/denom/4.0;
      }
    }
  }
  return output;
}




// [[Rcpp::export]]

NumericVector den_fast(NumericVector x, NumericVector p, int n, int m, double h){
  int ord = 1;
  NumericVector output(m);
  double denom;
  NumericMatrix L(ord + 1, n);
  NumericMatrix R(ord + 1, n);
  for(int i=0; i<=ord; i++){
    L(i,0) = pow(-x[0], i);
  }
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      L(j,i) = pow(-x[i],j) + exp((x[i-1]-x[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((x[n-i-1]-x[n-i])/h)*(pow(x[n-i],j)+R(j,n-i));
    }
  }
  int count;
  NumericVector counts(m);
  count = 0;
  for(int i=0; i<m; i++){
    if(p[i]>=x[n-1]){
      for(int j=i; j<m; j++) counts[j] = n;
      break;
    }
    else{
      while(x[count]<=p[i]) count += 1;
      counts[i] = count;
    } 
  }
  for(int orddo=0; orddo<=ord; orddo++){
    NumericVector coefs(orddo + 1);
    coefs[0] = coefs[orddo] = 1;
    if(orddo>1){
      double num = 1;
      for(int j=2; j<=orddo; j++) num *= j;
      double denom1 = 1;
      double denom2 = num/orddo;
      for(int i=2; i<=orddo; i++){
        coefs[i-1] = num/denom1/denom2;
        denom1 *= i;
        denom2 /= (orddo-i+1);
      }
    }
    denom = pow(h, orddo+1)*n;
    double exp_mult;
    for(int i=0; i<m; i++){
      if(counts[i]==0){
        exp_mult = exp((p[i]-x[0])/h);
        output[i] += .25*pow(x[0]-p[i], orddo)/denom*exp_mult;
        for(int j=0; j<=orddo; j++) output[i] += .25*coefs[j]*pow(-p[i],orddo-j)*R(j,0)/denom*exp_mult;
      }
      else{
        exp_mult = exp((x[counts[i]-1]-p[i])/h);
        for(int j=0; j<=orddo; j++) output[i] += .25*coefs[j]*(pow(p[i],orddo-j)*L(j,counts[i]-1)*exp_mult+pow(-p[i],orddo-j)*R(j,counts[i]-1)/exp_mult)/denom;
      }
    }
  }
  for(int i=0; i<m; i++){
    if(counts[i] > 0){
      if(dbl_abs(p[i]-x[counts[i]-1])<.0000000000000000000000000001) output[i] = (n+0.0)/(n-1.0)*(output[i]-.25/n/h);
    }
  }
  return output;
  //NumericVector out(n);
  //for(int i=0; i<n; i++) out[i] = L(0,i);
  //return out;
}




// [[Rcpp::export]]

NumericVector den_fast_all(NumericVector x, IntegerVector y, int n, IntegerVector ns, int nc, int maxn, NumericVector hs){
  NumericMatrix L0(maxn,nc);
  NumericMatrix L1(maxn,nc);
  NumericMatrix R0(maxn,nc);
  NumericMatrix R1(maxn,nc);
  IntegerVector counts(nc);
  IntegerVector counts2(nc);
  for(int i=0; i<nc; i++){
    counts[i] = 0;
    counts2[i] = 0;
  }
  IntegerVector previous(nc);
  IntegerVector next(nc);
  NumericVector first(nc);
  int cl;
  for(int i=0; i<n; i++){
    cl = y[i];
    if(counts[cl]==0){
      L0(0,cl) = 1.0;
      L1(0,cl) = -x[i];
      first[cl] = x[i];
    }
    else{
      L0(counts[cl],cl) = 1.0 + exp((x[previous[cl]]-x[i])/hs[cl])*L0(counts[cl]-1,cl);
      L1(counts[cl],cl) = -x[i] + exp((x[previous[cl]]-x[i])/hs[cl])*L1(counts[cl]-1,cl);
    }
    previous[cl] = i;
    counts[cl] += 1;
    cl = y[n-i-1];
    if(counts2[cl]>0){
      R0(ns[cl]-counts2[cl]-1,cl) = exp((x[n-i-1]-x[next[cl]])/hs[cl])*(1.0+R0(ns[cl]-counts2[cl],cl));
      R1(ns[cl]-counts2[cl]-1,cl) = exp((x[n-i-1]-x[next[cl]])/hs[cl])*(x[next[cl]]+R1(ns[cl]-counts2[cl],cl));
    }
    next[cl] = n-i-1;
    counts2[cl] += 1;
  }
  for(int i=0; i<nc; i++){
    counts[i] = 0;
    previous[i] = 0;
  }
  NumericMatrix output(n,nc);
  double exp_mult;
  for(int i=0; i<n; i++){
    cl = y[i];
    for(int c=0; c<nc; c++){
      if(c==cl){
        output(i,c) += .25/ns[c]/hs[c]*(L0(counts[c],c)*(1.0+x[i]/hs[c])+L1(counts[c],c)/hs[c] + R0(counts[c],c)*(1.0-x[i]/hs[c]) + R1(counts[c],c)/hs[c]);
      }
      else if(counts[c]==0){
        exp_mult = exp((x[i]-first[c])/hs[c]);
        output(i,c) += .25/ns[c]/hs[c]*exp_mult*(1.0+(first[c]-x[i])/hs[c]);
        output(i,c) += .25/ns[c]/hs[c]*exp_mult*(R0(0,c)*(1.0-x[i]/hs[c]) + R1(0,c)/hs[c]);
      }
      else{
        exp_mult = exp((x[previous[c]]-x[i])/hs[c]);
        output(i,c) += .25/ns[c]/hs[c]*((L0(counts[c]-1,c)*(1.0+x[i]/hs[c])+L1(counts[c]-1,c)/hs[c])*exp_mult + (R0(counts[c]-1,c)*(1.0-x[i]/hs[c]) + R1(counts[c]-1,c)/hs[c])/exp_mult);
      }
    }
    counts[cl] += 1;
    previous[cl] = i;
  }
  for(int i=0; i<n; i++){
    output(i,y[i]) = ns[y[i]]/(ns[y[i]]-1.0)*(output(i,y[i])-.25/ns[y[i]]/hs[y[i]]);
  }
  return output;
  //return L0;
}



// [[Rcpp::export]]

NumericMatrix dksum_cross(NumericVector x, NumericMatrix w, IntegerVector y, int n, IntegerVector ns, int nc, int maxn, NumericVector hs){
  NumericMatrix L0(maxn,nc);
  NumericMatrix L1(maxn,nc);
  NumericMatrix R0(maxn,nc);
  NumericMatrix R1(maxn,nc);
  IntegerVector counts(nc);
  IntegerVector counts2(nc);
  for(int i=0; i<nc; i++){
    counts[i] = 0;
    counts2[i] = 0;
  }
  IntegerVector previous(nc);
  IntegerVector next(nc);
  NumericVector first(nc);
  int cl;
  for(int i=0; i<n; i++){
    cl = y[i];
    if(counts[cl]==0){
      L0(0,cl) = w(i,cl);
      L1(0,cl) = -x[i]*w(i,cl);
      first[cl] = x[i];
    }
    else{
      L0(counts[cl],cl) = w(i,cl) + exp((x[previous[cl]]-x[i])/hs[cl])*L0(counts[cl]-1,cl);
      L1(counts[cl],cl) = -x[i]*w(i,cl) + exp((x[previous[cl]]-x[i])/hs[cl])*L1(counts[cl]-1,cl);
    }
    previous[cl] = i;
    counts[cl] += 1;
    cl = y[n-i-1];
    if(counts2[cl]>0){
      R0(ns[cl]-counts2[cl]-1,cl) = exp((x[n-i-1]-x[next[cl]])/hs[cl])*(w(next[cl],cl)+R0(ns[cl]-counts2[cl],cl));
      R1(ns[cl]-counts2[cl]-1,cl) = exp((x[n-i-1]-x[next[cl]])/hs[cl])*(x[next[cl]]*w(next[cl],cl)+R1(ns[cl]-counts2[cl],cl));
    }
    next[cl] = n-i-1;
    counts2[cl] += 1;
  }
  for(int i=0; i<nc; i++){
    counts[i] = 0;
    previous[i] = 0;
  }
  NumericMatrix output(n,nc);
  double exp_mult;
  for(int i=0; i<n; i++){
    cl = y[i];
    for(int c=0; c<nc; c++){
      if(c==cl){
        output(i,c) -= .25/hs[c]*(L0(counts[c],c)*x[i]+L1(counts[c],c) + R0(counts[c],c)*x[i] - R1(counts[c],c));
      }
      else if(counts[c]==0){
        exp_mult = exp((x[i]-first[c])/hs[c]);
        output(i,c) -= .25/hs[c]*exp_mult*(first[c]-x[i])*w(first[c],c);
        output(i,c) += .25/hs[c]*exp_mult*(R1(0,c)-x[i]*R0(0,c));
      }
      else{
        exp_mult = exp((x[previous[c]]-x[i])/hs[c]);
        output(i,c) -= .25/hs[c]*((L0(counts[c]-1,c)*x[i]+L1(counts[c]-1,c))*exp_mult + (R0(counts[c]-1,c)*x[i] - R1(counts[c]-1,c))/exp_mult);
      }
    }
    counts[cl] += 1;
    previous[cl] = i;
  }
  /*for(int i=0; i<n; i++){
    cl = y[i];
    for(int c=0; c<nc; c++){
      if(c==cl){
        output(i,c) -= .25/ns[c]/hs[c]/hs[c]*(L0(counts[c],c)*x[i]+L1(counts[c],c) + R0(counts[c],c)*x[i] - R1(counts[c],c));
      }
      else if(counts[c]==0){
        exp_mult = exp((x[i]-first[c])/hs[c]);
        output(i,c) -= .25/ns[c]/hs[c]*exp_mult*(first[c]-x[i])*w(first[c],c)/hs[c];
        output(i,c) += .25/ns[c]/hs[c]*exp_mult*(R1(0,c)-x[i]*R0(0,c))/hs[c];
      }
      else{
        exp_mult = exp((x[previous[c]]-x[i])/hs[c]);
        output(i,c) -= .25/ns[c]/hs[c]/hs[c]*((L0(counts[c]-1,c)*x[i]+L1(counts[c]-1,c))*exp_mult + (R0(counts[c]-1,c)*x[i] - R1(counts[c]-1,c))/exp_mult);
      }
    }
    counts[cl] += 1;
    previous[cl] = i;
  }*/
  return output;
  //return L0;
}








// [[Rcpp::export]]

NumericVector dksum_all(NumericVector x, NumericMatrix w, IntegerVector y, int n, int nc, NumericVector hs){
  NumericMatrix L0(n,nc);
  NumericMatrix L1(n,nc);
  NumericMatrix R0(n,nc);
  NumericMatrix R1(n,nc);
  for(int c=0; c<nc; c++){
    L0(0,c) = w(0,c);
    L1(0,c) = -x[0]*w(0,c);
  }
  for(int i=1; i<n; i++){
    for(int c=0; c<nc; c++){
      L0(i,c) = w(i,c) + exp((x[i-1]-x[i])/hs[c])*L0(i-1,c);
      L1(i,c) = -x[i]*w(i,c) + exp((x[i-1]-x[i])/hs[c])*L1(i-1,c);
      R0(n-i-1,c) = exp((x[n-i-1]-x[n-i])/hs[c])*(w(n-i,c)+R0(n-i,c));
      R1(n-i-1,c) = exp((x[n-i-1]-x[n-i])/hs[c])*(w(n-i,c)*x[n-i]+R1(n-i,c));
    }
  }
  NumericVector output(n);
  int cl;
  for(int i=0; i<n; i++){
    cl = y[i];
    output[i] -= .25/hs[cl]*(L0(i,cl)*x[i]+L1(i,cl) + R0(i,cl)*x[i] - R1(i,cl));
  }
  return output;
  //return L0;
}





// [[Rcpp::export]]

NumericVector dksum_all_in_class(NumericVector x, NumericMatrix w, IntegerVector y, int n, IntegerVector ns, int nc, int maxn, NumericVector hs){
  NumericMatrix L0(maxn,nc);
  NumericMatrix L1(maxn,nc);
  NumericMatrix R0(maxn,nc);
  NumericMatrix R1(maxn,nc);
  IntegerVector counts(nc);
  IntegerVector counts2(nc);
  for(int i=0; i<nc; i++){
    counts[i] = 0;
    counts2[i] = 0;
  }
  IntegerVector previous(nc);
  IntegerVector next(nc);
  NumericVector first(nc);
  int cl;
  for(int i=0; i<n; i++){
    cl = y[i];
    if(counts[cl]==0){
      L0(0,cl) = w(i,cl);
      L1(0,cl) = -x[i]*w(i,cl);
      first[cl] = x[i];
    }
    else{
      L0(counts[cl],cl) = w(i,cl) + exp((x[previous[cl]]-x[i])/hs[cl])*L0(counts[cl]-1,cl);
      L1(counts[cl],cl) = -x[i]*w(i,cl) + exp((x[previous[cl]]-x[i])/hs[cl])*L1(counts[cl]-1,cl);
    }
    previous[cl] = i;
    counts[cl] += 1;
    cl = y[n-i-1];
    if(counts2[cl]>0){
      R0(ns[cl]-counts2[cl]-1,cl) = exp((x[n-i-1]-x[next[cl]])/hs[cl])*(w(next[cl],cl)+R0(ns[cl]-counts2[cl],cl));
      R1(ns[cl]-counts2[cl]-1,cl) = exp((x[n-i-1]-x[next[cl]])/hs[cl])*(x[next[cl]]*w(next[cl],cl)+R1(ns[cl]-counts2[cl],cl));
    }
    next[cl] = n-i-1;
    counts2[cl] += 1;
  }
  for(int i=0; i<nc; i++){
    counts[i] = 0;
  }
  NumericVector output(n);
  double exp_mult;
  for(int i=0; i<n; i++){
    cl = y[i];
    output[i] -= .25/hs[cl]*(L0(counts[cl],cl)*x[i]+L1(counts[cl],cl) + R0(counts[cl],cl)*x[i] - R1(counts[cl],cl));
    counts[cl] += 1;
  }
  return output;
  //return L0;
}




// [[Rcpp::export]]

NumericVector den_fast_no_loo(NumericVector x, NumericVector p, int n, int m, double h){
  int ord = 1;
  NumericVector output(m);
  double denom;
  NumericMatrix L(ord + 1, n);
  NumericMatrix R(ord + 1, n);
  for(int i=0; i<=ord; i++){
    L(i,0) = pow(-x[0], i);
  }
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      L(j,i) = pow(-x[i],j) + exp((x[i-1]-x[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((x[n-i-1]-x[n-i])/h)*(pow(x[n-i],j)+R(j,n-i));
    }
  }
  int count;
  NumericVector counts(m);
  count = 0;
  for(int i=0; i<m; i++){
    if(p[i]>=x[n-1]){
      for(int j=i; j<m; j++) counts[j] = n;
      break;
    }
    else{
      while(x[count]<=p[i]) count += 1;
      counts[i] = count;
    } 
  }
  for(int orddo=0; orddo<=ord; orddo++){
    NumericVector coefs(orddo + 1);
    coefs[0] = coefs[orddo] = 1;
    if(orddo>1){
      double num = 1;
      for(int j=2; j<=orddo; j++) num *= j;
      double denom1 = 1;
      double denom2 = num/orddo;
      for(int i=2; i<=orddo; i++){
        coefs[i-1] = num/denom1/denom2;
        denom1 *= i;
        denom2 /= (orddo-i+1);
      }
    }
    denom = pow(h, orddo+1)*n;
    double exp_mult;
    for(int i=0; i<m; i++){
      if(counts[i]==0){
        exp_mult = exp((p[i]-x[0])/h);
        output[i] += .25*pow(x[0]-p[i], orddo)/denom*exp_mult;
        for(int j=0; j<=orddo; j++) output[i] += .25*coefs[j]*pow(-p[i],orddo-j)*R(j,0)/denom*exp_mult;
      }
      else{
        exp_mult = exp((x[counts[i]-1]-p[i])/h);
        for(int j=0; j<=orddo; j++) output[i] += .25*coefs[j]*(pow(p[i],orddo-j)*L(j,counts[i]-1)*exp_mult+pow(-p[i],orddo-j)*R(j,counts[i]-1)/exp_mult)/denom;
      }
    }
  }
  return output;
}


