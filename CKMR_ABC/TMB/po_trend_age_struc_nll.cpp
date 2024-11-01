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
  DATA_INTEGER(n_ages);
  DATA_INTEGER(data_start);
  DATA_INTEGER(data_end);
  DATA_SCALAR(phiA); //adult survival
  DATA_SCALAR(f); //fecundity
  DATA_MATRIX(n_match_PO_dibj);
  DATA_MATRIX(n_comp_PO_dibj);
  PARAMETER(n0_log);
  PARAMETER(phi0_logit);

  Type n0 = exp(n0_log);
  Type phi0 = invlogit(phi0_logit);

  matrix<Type> P(n_yrs,n_yrs);
  P.setZero();
  
  Type nll = 0.0;
  
  //population dynamics - run through twice (first time to get approx 
  // stable age proportions, second time for real using those proportions in year 1)
  matrix<Type> N_ya(n_yrs,n_ages);
  vector<Type> N_y(n_yrs); 
  N_ya(0,0)=n0;
  N_ya(0,1)=n0*phi0;
  for(int iage=2; iage<n_ages; iage++){
    N_ya(0,iage)=N_ya(0,iage-1)*phiA;
  }
  N_y(0)=N_ya(0,3)+N_ya(0,4);
  for(int iyr=1;iyr<n_yrs;iyr++){
    N_ya(iyr,1)=N_ya(iyr-1,0)*phi0;
    for(int iage =2;iage<n_ages;iage++){
      N_ya(iyr,iage)=N_ya(iyr-1,iage-1)*phiA;
    }
    N_y(iyr)=N_ya(iyr,3)+N_ya(iyr,4); //num breeders 
    N_ya(iyr,0)=N_y(iyr)*f;  //recruitment 
  }
  vector<Type> Age_mult(n_ages-1);
  for(int iage=1;iage<n_ages;iage++)Age_mult(iage-1)=N_ya(n_yrs-1,iage)/N_ya(n_yrs-1,iage-1);
  
  for(int iage=1; iage<n_ages; iage++){
    N_ya(0,iage)=N_ya(0,iage-1)*Age_mult(iage-1);
  }
  N_y(0)=N_ya(0,3)+N_ya(0,4);
  for(int iyr=1;iyr<n_yrs;iyr++){
    N_ya(iyr,1)=N_ya(iyr-1,0)*phi0;
    for(int iage =2;iage<n_ages;iage++){
      N_ya(iyr,iage)=N_ya(iyr-1,iage-1)*phiA;
    }
    N_y(iyr)=N_ya(iyr,3)+N_ya(iyr,4); //num breeders 
    N_ya(iyr,0)=N_y(iyr)*f;  //recruitment 
  }
  
  //data likelihood
  Type cur_prob=0;
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
  Type lambda = N_y(n_yrs-1)/N_y(n_yrs-2);
  
  REPORT(phi0);
  REPORT(N_ya);
  REPORT(N_y);
  REPORT(lambda);
  REPORT(cur_prob);
  REPORT(P);
  REPORT(nll);
  REPORT(N_mean);
  //ADREPORT(N_y);
  ADREPORT(lambda);
  ADREPORT(phi0);
  ADREPORT(N_mean);

  return nll;
}
