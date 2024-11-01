#  single sample simulation study
library(spatstat)
library(Matrix)
set.seed(12345)

N_opt = c(100,400)
n_sim = 1
N_sample = c(50,100)

quantile_ABC = 0.01  #take top 1% of sufficient statistic diferences for "simple" ABC algorithm
quantile_ABC2 = 0.05  # take top 1% for KDEs associated with Beaumont regression ABC algorithm

n_sim_ABC = 1000000

Ests = Var = array(NA,dim=c(2,2,5,n_sim))  #N_opt, offspring model, estimator,sim
Cover = array(0,dim=c(2,2,5,n_sim))
Epsilon = array(NA,dim=c(2,2,n_sim))  #suff stat tolerance for ABC-simple
Fec = array(NA,dim=c(2,2,n_sim))

bin_nlogL <- function(logN,n,x){
  p = 1/exp(logN)
  -(x * log(p) + (n-x)*log(1-p))
}
bin_nlogL_prof_ssq <- function(logN,n,x,obj){
  p = 1/exp(logN)
  cur_obj = -(x * log(p) + (n-x)*log(1-p))
  (cur_obj-(obj+1.353))^2  #for 90% CI - X^2_{1,0.9}
}

nlogL_hg <- function(logN,n_0,n_f,x){
  N_f = exp(logN)+n_0+n_f-x
  -(lgamma(N_f-n_0+1)+lgamma(N_f-n_f+1)-lgamma(N_f-n_0-n_f+x+1)-lgamma(N_f+1))
}

nlogL_hg_prof_ssq <- function(logN,n_0,n_f,x,obj){
  N_f = exp(logN)+n_0+n_f-x
  cur_obj =   -(lgamma(N_f-n_0+1)+lgamma(N_f-n_f+1)-lgamma(N_f-n_0-n_f+x+1)-lgamma(N_f+1))
  (cur_obj-(obj+1.353))^2
}  

post_mode <- function(Y,Wts=NULL){
  Dens = density(Y,weights=Wts)
  Dens$x[which(Dens$y==max(Dens$y))[1]]
}


start_time = Sys.time()

#set up all potential Y matrices for ABC sims
Y1_list = Y2_list = vector("list",2000)
for(i in 1:2000){
  Y1_list[[i]]=Diagonal(i)
  Y2_list[[i]]=Matrix(0,i*2,i)
  Y2_list[[i]][cbind(2*c(1:i),c(1:i))]=1
  Y2_list[[i]][cbind(2*c(1:i)-1,c(1:i))]=1
}

for(iN in 1:2){
  N_f = N_opt[iN]
  n_sample = N_sample[iN]
  Y1 = Y1_list[[N_opt[iN]]]
  Y2 = Y2_list[[N_opt[iN]]]
  for(irepro in 1:2){
    if(irepro==1){
      N_o = N_f
      Y = Y1
    }
    for(isim in 1:n_sim){
      if(irepro==2){
        N_o = 2*N_f
        Y=Y2
      }
      Sampled_moms = sample(c(1:N_f),n_sample)
      Sampled_o = sample(c(1:N_o),n_sample)
      Y_samp = Y[Sampled_o,Sampled_moms]
      n_pairs = sum(Y_samp)
      
      #Estimation
      
      #1) binomial pseudo-likelihood
      bin_est = optim(3,bin_nlogL,hessian=FALSE,lower=0,upper=13,n=n_sample^2,x=n_pairs,method="Brent")
      Ests[iN,irepro,1,isim]=exp(bin_est$par)
      min_log_L = bin_est$value
      #Var[iN,irepro,1,isim] = (Ests[iN,irepro,1,isim])^2 / bin_est$hessian  #accounting for transformation from log space

      #come down 1.353 on either side of logL for 90% profile interval
      prof_min = optimize(bin_nlogL_prof_ssq,interval=c(0,bin_est$par),n=n_sample^2,x=n_pairs,obj=min_log_L)
      prof_max = optimize(bin_nlogL_prof_ssq,interval=c(bin_est$par,13),n=n_sample^2,x=n_pairs,obj=min_log_L)
      if(exp(prof_min$minimum)<N_f & exp(prof_max$minimum)>N_f)Cover[iN,irepro,1,isim]=1
      
      #2) hyper-geometric
      hg_est = optim(3,nlogL_hg,hessian=FALSE,lower=0,upper=13,n_0=n_sample,n_f=n_sample,x=n_pairs,method="Brent")
      Ests[iN,irepro,2,isim]=exp(hg_est$par)+2*n_sample-n_pairs
      min_log_L = hg_est$value
      prof_min = optimize(nlogL_hg_prof_ssq,interval=c(0,hg_est$par),n_0=n_sample,n_f=n_sample,x=n_pairs,obj=min_log_L)
      prof_max = optimize(nlogL_hg_prof_ssq,interval=c(hg_est$par,13),n_0=n_sample,n_f=n_sample,x=n_pairs,obj=min_log_L)
      if((exp(prof_min$minimum)+2*n_sample-n_pairs)<N_f & (exp(prof_max$minimum)+2*n_sample-n_pairs)>N_f)Cover[iN,irepro,2,isim]=1
      
      #3) ABC methods
      N_prior = round(runif(n_sim_ABC,n_sample,5*N_f))
      Fec_prior = runif(n_sim_ABC,3,7)
      N_pairs = rep(NA,n_sim_ABC)
      for(iABC in 1:n_sim_ABC){
        if(irepro==1){
          N_o_sim = N_prior[iABC]
          Y_sim = Y1_list[[N_o_sim]]
        }
        
        if(irepro==2){
          N_o_sim = N_prior[iABC]*2
          Y_sim=Y2_list[[N_prior[iABC]]]
        }
        Sampled_moms = sample(c(1:N_prior[iABC]),n_sample)
        Sampled_o = sample(c(1:N_o_sim),n_sample)
        Y_samp = Y_sim[Sampled_o,Sampled_moms]
        N_pairs[iABC] = sum(Y_samp)
      }
      #simple estimator
      Eps_pairs = N_pairs-n_pairs
      SSQ_pairs = (Eps_pairs)^2
      low_SSQ = min(SSQ_pairs)
      quant1 = quantile(SSQ_pairs,quantile_ABC)
      quant2 = quantile(SSQ_pairs,quantile_ABC2)
      
      Which1 = which(N_pairs==n_pairs)  #now doing exact
      Which2 = which(SSQ_pairs<=quant2)
      Epsilon[iN,irepro,isim]=quantile(Eps_pairs,quantile_ABC2)
      
      N_ABC1 = N_prior[Which1]  
      Ests[iN,irepro,3,isim] = post_mode(N_ABC1)
      Var[iN,irepro,3,isim]=var(N_ABC1)
      lower = quantile(N_ABC1,0.05)
      upper = quantile(N_ABC1,0.95)
      if(lower<N_f & upper>N_f)Cover[iN,irepro,3,isim]=1
      
      Fec_ABC1 = Fec_prior[Which1]
      Fec[iN,irepro,isim]=median(Fec_prior)
      
      #simple estimator, 0.05 quantile
      N_ABC2 = N_prior[Which2]
      Ests[iN,irepro,4,isim]=post_mode(N_ABC2)
      lower = quantile(N_ABC2,0.05)
      upper = quantile(N_ABC2,0.95)
      if(lower<N_f & upper>N_f)Cover[iN,irepro,4,isim]=1
      
      
      # regression adjustment estimator
      N_red = N_prior[Which2]
      SSQ_red = SSQ_pairs[Which2]
      eps = max(SSQ_red)
      if(eps>0){  #o/w lm will be singular
        Wts = 1/eps * (1-(SSQ_red/eps)^2)
        lin_reg = lm(N_red ~ N_pairs[Which2],weights=Wts)
        E_N_given_S = predict(lin_reg)
        N_star = N_red - E_N_given_S + predict(lin_reg,new_data=n_pairs)
        ECDF = ewcdf(N_star,Wts/sum(Wts))  #non-exported function from spatstat
        Ests[iN,irepro,5,isim] =post_mode(N_star,Wts/sum(Wts))
        lower = quantile(ECDF,0.05)
        upper = quantile(ECDF,0.95)
        if(lower<N_f & upper>N_f)Cover[iN,irepro,5,isim]=1
        #Var[iN,irepro,4,isim]=var(sample(N_star,10000,replace=TRUE,prob=Wts))
      }
      print(paste0("iN ",iN," irepro ",irepro," isim ",isim,"\n"))
    }
  }
}

fin_time = Sys.time()
print(fin_time - start_time)
