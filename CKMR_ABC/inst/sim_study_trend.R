# Simulation study 2: population trend

#w/ parallel, 1000 sims @ 12 clusters and n_ABC_sims=12000 takes 7.7 days
#on my laptop

library(fishSim)
library(TMB) 
library(spatstat)
library(parallel)
library(abc)

st_time = Sys.time()

set.seed(1111)
n_sim = 100
n_ABC_sims = 24000  
n_clusters=12
cl <- makeCluster(n_clusters)

mat_age = 2 #sexual maturity age

N_init = 200
ageMort = 1-c(0.65,0.9,0.9,0.9,0)
L = matrix(0,5,5)
L[1,]=c(0,0,.9,.9,0)
L[2,1]=0.65
L[3,2]=L[4,3]=L[5,4]=0.9
Age_props = eigen(L)$vectors[,1]/sum(eigen(L)$vectors[,1])  #stable age structure proportions

Repro = c(0,0,0,1,1)  #fishsim adds a year to maturity vector 

#compile TMB estimation files
TmbFile = "c:/users/paul.conn/git/ckmr/ckmr_abc/ckmr_abc/TMB/po_trend_rem_nll.cpp"
compile(TmbFile )
TmbExec="c:/users/paul.conn/git/ckmr/ckmr_abc/ckmr_abc/TMB/po_trend_rem_nll"
dyn.load(dynlib(TmbExec))

TmbFile = "c:/users/paul.conn/git/ckmr/ckmr_abc/ckmr_abc/TMB/po_trend_age_struc_nll.cpp"
compile(TmbFile )
TmbExec="c:/users/paul.conn/git/ckmr/ckmr_abc/ckmr_abc/TMB/po_trend_age_struc_nll"
dyn.load(dynlib(TmbExec))

n_ages=5
n_yrs =20
N_yrs_sample = c(11:20)
n_yrs_sample = length(N_yrs_sample)
N_breed = rep(0,n_yrs)
samples_per_year = 25
N_ta = matrix(0,n_yrs,n_ages-1)
N_ta2 = matrix(0,n_yrs,n_ages)
Yrs = c(1:n_yrs)

N_breed_sim = matrix(NA,n_sim,n_yrs)
Lambda_est = Lambda_SE = N_breed_est = Cover_lambda=Cover_Nbreed= matrix(NA,n_sim,3)  #3 estimators: pseudo-exp, pseudo-age, ABC
Lambda_realized = rep(NA,n_sim)  # realized lambda for each sim

ABC_FT = rep(NA,n_ABC_sims)  
ABC_pars = matrix(0,n_ABC_sims,2) #N_init, phi0
Lambda_ABC = N_breed_ABC = rep(0,n_ABC_sims)

sum_sq_lambda <- function(par,N){  #estimate lambda using realized N for each simulation
  lambda = exp(par[1])
  N0 = exp(par[2])
  n_yrs = length(N)
  Nvec = rep(0,n_yrs)
  Nvec[1]=N0
  for(iyr in 2:n_yrs)Nvec[iyr]=Nvec[iyr-1]*lambda
  sum((Nvec-N)^2)
}

post_mode <- function(Y,Wts=NULL){
  Dens = density(Y,weights=Wts)
  Dens$x[which(Dens$y==max(Dens$y))[1]]
}

sim_ABC <- function(inputs){
  library(fishSim)
  N_init=inputs$N_init
  curMort=inputs$curMort
  phi0=1-inputs$curMort[1]
  L = inputs$L
  n_yrs=inputs$Const$n_yrs
  n_ages=inputs$Const$n_ages
  N_yrs_sample = inputs$Const$N_yrs_sample
  n_yrs_sample = inputs$Const$n_yrs_sample  
  N_breed = inputs$Const$N_breed
  samples_per_year=inputs$Const$samples_per_year
  N_ta=inputs$Const$N_ta
  N_ta2=inputs$Const$N_ta2
  Repro=inputs$Const$Repro
  mat_age=inputs$Const$mat_age
  Yrs=inputs$Const$Yrs
  
  Age_props = eigen(L)$vectors[,1]/sum(eigen(L)$vectors[,1])  #stable age structure proportions
  
  archive <- make_archive()
  
  indiv = makeFounders(pop=N_init,osr=c(0.5,0.5),stocks=1,minAge=1,maxAge=5,survCurv=Age_props)
  flag_ABC=0
  for(iyr in 1:n_yrs){
    if(flag_ABC==0){
      indiv <- mort(indiv = indiv, type = "age", ageMort = c(0,curMort), year = iyr,maxAge=n_ages)
      archive <- archive_dead(indiv = indiv, archive = archive)
      indiv <- remove_dead(indiv = indiv)
      if(sum(indiv$Sex=="F" & indiv$AgeLast%in% c(3,4))==0 | sum(indiv$Sex=="M" & indiv$AgeLast%in% c(3,4))==0)flag_ABC=1 #ran out of mature breeders
      if(flag_ABC==0){
        N_breed[iyr] = sum(indiv$AgeLast %in% c(3,4))
        N_ta[iyr,]=summary(factor(indiv$AgeLast,levels=c('1','2','3','4')))
        indiv <- altMate(indiv = indiv, type='ageSex',maleCurve=Repro,femaleCurve=Repro,osr = c(0.5,0.5),batchSize=2, fecundityDist="poisson",year = iyr)
        indiv <- birthdays(indiv = indiv)
        N_ta2[iyr,]=summary(factor(indiv$AgeLast,levels=c('1','2','3','4','5')))
      }
    }
  }
  # sample from archive (dead only)
  for(isampyr in 1:n_yrs_sample){
    cur_yr = N_yrs_sample[isampyr]
    Cur_which = which(archive$DeathY==cur_yr)
    if(length(Cur_which)>=samples_per_year)
      archive[sample(Cur_which,samples_per_year),"SampY"]=N_yrs_sample[isampyr]
    else flag_ABC=1
  }
  
  if(flag_ABC==0){
    indiv=rbind(archive,indiv)
    
    POPs <- quickin(indiv,max_gen=2)[["POP"]]
    
    ###compile statistics for pseudo-likelihood trend model
    # here, 'si', and 'bj' index sex of older animal, and birth of the younger animal
    # as far as dimensions, there are 5 years of data, but birth dates can be earlier (up one generation time)
    n_comp_PO_sidibj=n_comp_PO_sidibj=n_match_PO_sidibj=n_match_PO_sidibj=array(0,dim=c(2,n_yrs,n_yrs))  #suff stat
    n_comp_PO_bidisibj=n_comp_PO_sidisibj=n_match_PO_bidisibj=n_match_PO_bidisibj=array(0,dim=c(n_yrs,n_yrs,2,n_yrs)) #for debugging
    for(i in 1:nrow(POPs)){
      indiv1 = indiv[which(indiv[,1]==POPs[i,1]),]
      indiv2 = indiv[which(indiv[,1]==POPs[i,2]),]
      
      FIRST = (indiv1$BirthY < indiv2$BirthY)
      # case 1: parent is animal 1
      if(FIRST){
        p_by = indiv1$BirthY
        p_sex = (indiv1$Sex=="M")+1  #1 for female, 2 for males
        p_sy = indiv1$SampY
        o_by = indiv2$BirthY  #year 11 in simulation = year 1 for CKMR inference
      }
      else{
        p_by = indiv2$BirthY
        p_sex = (indiv2$Sex=="M")+1  #1 for female, 2 for males
        p_sy = indiv2$SampY
        o_by = indiv1$BirthY  #year 11 in simulation = year 1 for CKMR inference
      }
      n_match_PO_sidibj[p_sex,p_sy,o_by]=n_match_PO_sidibj[p_sex,p_sy,o_by]+1 #revisit / think about maturity, sperm storage, repro timing
      n_match_PO_bidisibj[p_by,p_sy,p_sex,o_by]=n_match_PO_bidisibj[p_by,p_sy,p_sex,o_by]+1 
    }
    
    Sampled = indiv[is.na(indiv$SampY)==FALSE,]
    n_indiv = nrow(Sampled)
    for(i1 in 1:(n_indiv-1)){
      for(i2 in (i1+1):n_indiv){
        born_diff = Sampled$BirthY[i1] - Sampled$BirthY[i2]
        if(born_diff<0){ #first animal is older
          p_by = Sampled$BirthY[i1]
          p_sy = Sampled$SampY[i1]
          p_sex = (Sampled$Sex[i1]=="M")+1  #1 for female, 2 for males
          o_by = Sampled$BirthY[i2]    
        }
        else{
          if(born_diff>0){ #note no comparisons for animals born in same year
            p_by = Sampled$BirthY[i2]
            p_sy = Sampled$SampY[i2]
            p_sex = (Sampled$Sex[i2]=="M")+1  #1 for female, 2 for males
            o_by = Sampled$BirthY[i1]    
          }
        }
        if(o_by>(p_by+mat_age) & o_by<p_sy){
          n_comp_PO_sidibj[p_sex,p_sy,o_by]=n_comp_PO_sidibj[p_sex,p_sy,o_by]+1 #revisit / think about maturity, sperm storage, repro timing
          n_comp_PO_bidisibj[p_by,p_sy,p_sex,o_by]=n_comp_PO_bidisibj[p_by,p_sy,p_sex,o_by]+1 
        }
      }
    }
    
    #collapse over sex since dynamics the same
    n_comp_PO_dibj = apply(n_comp_PO_sidibj,c(2,3),'sum')
    n_match_PO_dibj = apply(n_match_PO_sidibj,c(2,3),'sum')
    
    Obs_prop = n_match_PO_dibj/n_comp_PO_dibj
    mean_prop = mean(Obs_prop,na.rm=TRUE)
    Prop_bj = colMeans(Obs_prop,na.rm=TRUE)
    Which_yrs = which(!is.na(Prop_bj))
    slope = as.numeric(lm(Prop_bj[Which_yrs]~Yrs[Which_yrs])$coefficients[2])
  }
  if(flag_ABC==1)return(NULL)
  else return(list("n_match_PO_dibj"=n_match_PO_dibj,"n_comp_PO_dibj"=n_comp_PO_dibj,"slope"=slope,"N_mean"=mean(N_breed)))
}


for(isim in 1:n_sim){
  cat(paste("isim",isim,"\n"))
  flag=0
  archive <- make_archive()
  indiv = makeFounders(pop=N_init,osr=c(0.5,0.5),stocks=1,minAge=1,maxAge=5,survCurv=Age_props)
  for(iyr in 1:n_yrs){
    indiv <- mort(indiv = indiv, type = "age", ageMort = c(0,ageMort), year = iyr,maxAge=n_ages)
    archive <- archive_dead(indiv = indiv, archive = archive)
    indiv <- remove_dead(indiv = indiv)
    N_breed[iyr] = sum(indiv$AgeLast %in% c(3,4))
    N_ta[iyr,]=summary(factor(indiv$AgeLast,levels=c('1','2','3','4')))
    indiv <- altMate(indiv = indiv, type='ageSex',maleCurve=Repro,femaleCurve=Repro,osr = c(0.5,0.5),batchSize=2, fecundityDist="poisson",year = iyr)
    indiv <- birthdays(indiv = indiv)
    N_ta2[iyr,]=summary(factor(indiv$AgeLast,levels=c('1','2','3','4','5')))
    
  }
  N_breed_sim[isim,] = N_breed
  lam_mod <- nlminb(c(0,log(100)),sum_sq_lambda,N=N_breed)
  Lambda_realized[isim] = exp(lam_mod$par[1])
  Nbreed_mean = mean(N_breed)
  
  
  # sample from archive (dead only)
  for(isampyr in 1:n_yrs_sample){
    cur_yr = N_yrs_sample[isampyr]
    Cur_which = which(archive$DeathY==cur_yr)
    if(length(Cur_which)>=samples_per_year){
      archive[sample(Cur_which,samples_per_year),"SampY"]=N_yrs_sample[isampyr]
    }
    else flag=1
  }
  if(flag==0){
    indiv=rbind(archive,indiv)
    
    POPs <- quickin(indiv,max_gen=2)[["POP"]]
    
    ###compile statistics for pseudo-likelihood trend model
    # here, 'si', and 'bj' index sex of older animal, and birth of the younger animal
    # as far as dimensions, there are 5 years of data, but birth dates can be earlier (up one generation time)
    n_comp_PO_sidibj=n_comp_PO_sidibj=n_match_PO_sidibj=n_match_PO_sidibj=array(0,dim=c(2,n_yrs,n_yrs))  #suff stat
    n_comp_PO_bidisibj=n_comp_PO_sidisibj=n_match_PO_bidisibj=n_match_PO_bidisibj=array(0,dim=c(n_yrs,n_yrs,2,n_yrs)) #for debugging
    for(i in 1:nrow(POPs)){
      indiv1 = indiv[which(indiv[,1]==POPs[i,1]),]
      indiv2 = indiv[which(indiv[,1]==POPs[i,2]),]
      
      FIRST = (indiv1$BirthY < indiv2$BirthY)
      # case 1: parent is animal 1
      if(FIRST){
        p_by = indiv1$BirthY
        p_sex = (indiv1$Sex=="M")+1  #1 for female, 2 for males
        p_sy = indiv1$SampY
        o_by = indiv2$BirthY  #year 11 in simulation = year 1 for CKMR inference
      }
      else{
        p_by = indiv2$BirthY
        p_sex = (indiv2$Sex=="M")+1  #1 for female, 2 for males
        p_sy = indiv2$SampY
        o_by = indiv1$BirthY  #year 11 in simulation = year 1 for CKMR inference
      }
      n_match_PO_sidibj[p_sex,p_sy,o_by]=n_match_PO_sidibj[p_sex,p_sy,o_by]+1 #revisit / think about maturity, sperm storage, repro timing
      n_match_PO_bidisibj[p_by,p_sy,p_sex,o_by]=n_match_PO_bidisibj[p_by,p_sy,p_sex,o_by]+1 
    }
    
    Sampled = indiv[is.na(indiv$SampY)==FALSE,]
    n_indiv = nrow(Sampled)
    for(i1 in 1:(n_indiv-1)){
      for(i2 in (i1+1):n_indiv){
        born_diff = Sampled$BirthY[i1] - Sampled$BirthY[i2]
        if(born_diff<0){ #first animal is older
          p_by = Sampled$BirthY[i1]
          p_sy = Sampled$SampY[i1]
          p_sex = (Sampled$Sex[i1]=="M")+1  #1 for female, 2 for males
          o_by = Sampled$BirthY[i2]    
        }
        else{
          if(born_diff>0){ #note no comparisons for animals born in same year
            p_by = Sampled$BirthY[i2]
            p_sy = Sampled$SampY[i2]
            p_sex = (Sampled$Sex[i2]=="M")+1  #1 for female, 2 for males
            o_by = Sampled$BirthY[i1]    
          }
        }
        if(o_by>(p_by+mat_age) & o_by<p_sy){
          n_comp_PO_sidibj[p_sex,p_sy,o_by]=n_comp_PO_sidibj[p_sex,p_sy,o_by]+1 #revisit / think about maturity, sperm storage, repro timing
          n_comp_PO_bidisibj[p_by,p_sy,p_sex,o_by]=n_comp_PO_bidisibj[p_by,p_sy,p_sex,o_by]+1 
        }
      }
    }
    
    #collapse over sex since dynamics the same
    n_comp_PO_dibj = apply(n_comp_PO_sidibj,c(2,3),'sum')
    n_match_PO_dibj = apply(n_match_PO_sidibj,c(2,3),'sum')
    
    #trend model, no age stucture (just model number of breeders)
    Data=list("n_yrs"=n_yrs,"data_start"=9,"data_end"=19,  #note C convention in start, end
              "n_match_PO_dibj"=n_match_PO_dibj,"n_comp_PO_dibj"=n_comp_PO_dibj)
    Params = list("n0_log"=log(mean(N_breed)),"lambda_log"=0) #intial param values
    
    Map = list()  #specify fixed parameter values
    #Random= c("log_eta1","log_eta2","log_eta3")
    Random=NULL
    
    Obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, DLL="po_trend_rem_nll")
    
    Obj$fn( Obj$par )
    Report = Obj$report()  
    
    #Minimize negative log likelihood and time it
    Start_time = Sys.time()
    Opt = nlminb(start=Params, objective=Obj$fn, gradient=Obj$gr)
    End_time = Sys.time()
    
    Report=Obj$report()
    SD_report=sdreport(Obj)
    N_est = SD_report$value[which(names(SD_report$value)=="N_y")]
    lambda_est=SD_report$value[which(names(SD_report$value)=="lambda")]
    N_breed_est[isim,1]=Report$N_mean
    Lambda_est[isim,1]=lambda_est
    Lambda_SE[isim,1]=SD_report$sd[which(names(SD_report$value)=="lambda")]
    cv_lambda = Lambda_SE[isim,1]/lambda_est
    c_lambda = exp(1.645*sqrt(log(1+cv_lambda^2)))
    CI_lambda = c(lambda_est/c_lambda,lambda_est*c_lambda)
    Cover_lambda[isim,1]= 1*(CI_lambda[1]<Lambda_realized[isim] & CI_lambda[2]>Lambda_realized[isim])
    Nbreed_SE = SD_report$sd[which(names(SD_report$value)=="N_mean")]
    cv_Nbreed = Nbreed_SE/Report$N_mean
    c_Nbreed = exp(1.645*sqrt(log(1+cv_Nbreed^2)))
    CI_Nbreed = c(Report$N_mean/c_Nbreed,Report$N_mean*c_Nbreed)
    Cover_Nbreed[isim,1]= 1*(CI_Nbreed[1]<Nbreed_mean & CI_Nbreed[2]>Nbreed_mean)
    
    
    #trend model, age stucture (model pup survival instead of lambda)
    Data=list("n_yrs"=n_yrs,"n_ages"=5,"data_start"=9,"data_end"=19,  #note C convention in start, end
              "phiA"=0.9,"f"=1,"n_match_PO_dibj"=n_match_PO_dibj,"n_comp_PO_dibj"=n_comp_PO_dibj)
    Params = list("n0_log"=log(mean(N_breed)),"phi0_logit"=log(0.65/0.35)) #intial param values
    
    Map = list()  #specify fixed parameter values
    #Random= c("log_eta1","log_eta2","log_eta3")
    Random=NULL
    
    Obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, DLL="po_trend_age_struc_nll")
    
    Obj$fn( Obj$par )
    
    #Minimize negative log likelihood and time it
    Start_time = Sys.time()
    Opt = nlminb(start=Params, objective=Obj$fn, gradient=Obj$gr)
    End_time = Sys.time()
    
    Report=Obj$report()
    SD_report=sdreport(Obj)
    N_est = SD_report$value[which(names(SD_report$value)=="N_y")]
    lambda_est=SD_report$value[which(names(SD_report$value)=="lambda")]
    N_breed_est[isim,2]=Report$N_mean
    Lambda_est[isim,2]=lambda_est
    Lambda_SE[isim,2]=SD_report$sd[which(names(SD_report$value)=="lambda")]
    cv_lambda = Lambda_SE[isim,2]/lambda_est
    c_lambda = exp(1.645*sqrt(log(1+cv_lambda^2)))
    CI_lambda = c(lambda_est/c_lambda,lambda_est*c_lambda)
    Cover_lambda[isim,2]= 1*(CI_lambda[1]<Lambda_realized[isim] & CI_lambda[2]>Lambda_realized[isim])
    Nbreed_SE = SD_report$sd[which(names(SD_report$value)=="N_mean")]
    cv_Nbreed = Nbreed_SE/Report$N_mean
    c_Nbreed = exp(1.645*sqrt(log(1+cv_Nbreed^2)))
    CI_Nbreed = c(Report$N_mean/c_Nbreed,Report$N_mean*c_Nbreed)
    Cover_Nbreed[isim,2]= 1*(CI_Nbreed[1]<Nbreed_mean & CI_Nbreed[2]>Nbreed_mean)
    
    #now the ABC part
    #summary stats for real data
    N_prop = n_match_PO_dibj/n_comp_PO_dibj
    Which_comp = which(n_comp_PO_dibj>0)

    ABC_inputs_list = vector("list",n_ABC_sims)
    ABC_pars = matrix(0,n_ABC_sims,2)
    #ABC_pars[,1]=runif(n_ABC_sims,50,800)
    for(i in 1:n_ABC_sims){  #N~1/N Jeffreys prior
      done = 0
      while(done==0){
        N_try = runif(1,50,800)
        if(runif(1)<(50/N_try)){
          ABC_pars[i,1]=N_try
          done=1
        }
      }
    }
    ABC_pars[,2]=runif(n_ABC_sims,0.4,0.9)
    Const = list("n_yrs"=n_yrs,"n_ages"=n_ages,"N_yrs_sample"=N_yrs_sample,
                 "n_yrs_sample"=n_yrs_sample,"N_breed"=N_breed,
                 "samples_per_year"=samples_per_year,"N_ta"=N_ta,"N_ta2"=N_ta2,
                 "Repro"=Repro,"mat_age"=mat_age,"Yrs"=Yrs)

    for(i in 1:n_ABC_sims){
      L[2,1]=ABC_pars[i,2]
      Lambda_ABC[i]=eigen(L)$values[1]
      ABC_inputs_list[[i]]=list(N_inits=ABC_pars[i,1],
                                curMort=c(1-ABC_pars[i,2],ageMort[-1]),
                                L=L,
                                Const=Const)
    }
      
    sim_out = parSapply(cl,ABC_inputs_list,sim_ABC)
    
    Stat_mat = matrix(0,n_ABC_sims,200)
    No_est = rep(0,n_ABC_sims)
    
    for(i in 1:n_ABC_sims){
      if(is.null(sim_out[[i]])==FALSE){
        N_breed_ABC[i]=sim_out[[i]]$N_mean
        Stat_mat[i,]=as.vector(sim_out[[i]]$n_match_PO_dibj[11:20,])
      }
      else No_est[i]=1
    }
    
    Which_no_est = which(No_est==1)
    Par_dm = data.frame("Lambda"=Re(Lambda_ABC),"Nbreed"=N_breed_ABC)
    Par_dm = Par_dm[-Which_no_est,]
    Stat_mat= Stat_mat[-Which_no_est,]
    abc_out <- abc(target=as.vector(n_match_PO_dibj[11:20,]),param=Par_dm,tol=0.05,
                   sumstat=Stat_mat,transf=c("log","log"),method="neuralnet",
                   sizenet=4)
    summary_abc = summary(abc_out,print=FALSE,intvl=0.9)
    
    
    Lambda_est[isim,3]=summary_abc[5,1]
    N_breed_est[isim,3]=summary_abc[5,2]
    lower = summary_abc[2,1]
    upper = summary_abc[6,1]
    if(Lambda_est[isim,3]>lower & Lambda_est[isim,3]<upper)Cover_lambda[isim,3]=1
    lower = summary_abc[2,2]
    upper = summary_abc[6,2]
    if(N_breed_est[isim,3]>lower & N_breed_est[isim,3]<upper)Cover_Nbreed[isim,3]=1
  }
  
  if((isim %% 10)==0)save.image("Sim_trend_nn.RData")
}

Sys.time()-st_time


Bias_N = (N_breed_est-rowMeans(N_breed_sim))/rowMeans(N_breed_sim)
Bias_Lambda=Lambda_est-Lambda_realized/Lambda_realized
colMeans(Bias_N,na.rm=TRUE)
colMeans(Bias_Lambda,na.rm=TRUE)
colMeans(Cover_Nbreed,na.rm=TRUE)
colMeans(Cover_lambda,na.rm=TRUE)



