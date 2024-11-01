// simple PO CKMR with exponential pop growth and known ages
// this version outputs adreport on terminal abundance
#include <TMB.hpp>

template<class Type>
Type dbinom_kern_log(Type n, Type x, Type p){
  Type p0 = p==0;
  return x*log(p+p0)+(n-x)*log(1.0-p);
}

template<class Type>
Type get_PO_prob(int psy, int oby, vector<Type> N_yr){
  // parent sampling year, offspring birth year, abundance vector
  return 2.0 / N_yr(oby);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_yrs);
  DATA_INTEGER(data_start);
  DATA_INTEGER(data_end);
  DATA_MATRIX(n_match_PO_dibj);
  DATA_MATRIX(n_comp_PO_dibj);
  PARAMETER(n0_log);
  PARAMETER(lambda_log);

  Type n0 = exp(n0_log);
  Type lambda = exp(lambda_log);

  matrix<Type> P(n_yrs,n_yrs);
  P.setZero();
  
  Type nll = 0.0;
  vector<Type> N_y(n_yrs); 
  Type cur_prob=0;
  N_y(0)=n0;
  for(int iyr=1;iyr<n_yrs;iyr++)N_y(iyr)=N_y(iyr-1)*lambda;

  for(int jyr=0;jyr<n_yrs; jyr++){
    for(int isamp=data_start;isamp<=data_end;isamp++){
      if(n_comp_PO_dibj(isamp,jyr)>0){
        cur_prob = get_PO_prob(isamp,jyr,N_y);
        P(jyr,isamp)=cur_prob;
        nll -= dbinom_kern_log(n_comp_PO_dibj(isamp,jyr),n_match_PO_dibj(isamp,jyr),cur_prob);
      }
    }
  }
  Type N_mean = N_y.mean();
  
  REPORT(N_y);
  REPORT(lambda);
  REPORT(cur_prob);
  REPORT(P);
  REPORT(nll);
  REPORT(N_mean);
  //ADREPORT(N_y);
  ADREPORT(lambda);
  ADREPORT(N_mean);

  return nll;
}
