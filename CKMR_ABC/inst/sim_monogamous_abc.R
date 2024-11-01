# test ABC approaches to CKMR on monogamous beaver-like population
# this version includes semi-automated ABC (Fernhead & Prangle 2012) as well as the 
# neural network approach of Blum and Francois

library(parallel)
library(abc)

logit <- function(x){
  log(x/(1-x))
}

birth.interval = function(Sampled,Kin){
  n_pairs = nrow(Kin)
  Byear = matrix(0,n_pairs,2)
  for(ipair in 1:n_pairs){
    for(ikin in 1:2){
      Byear[ipair,ikin] = Sampled[which(Sampled$Me==Kin[ipair,ikin]),"BirthY"]
    }
  }
  return(abs(Byear[,1]-Byear[,2]))
}

sum.sq.lambda <- function(par,N){  #for estimating "realized" lambda for each simulation using least squares
  lambda = exp(par[1])
  N0 = exp(par[2])
  n_yrs = length(N)
  Nvec = rep(0,n_yrs)
  Nvec[1]=N0
  for(iyr in 2:n_yrs)Nvec[iyr]=Nvec[iyr-1]*lambda
  sum((Nvec-N)^2)
}

gen.formula.S <- function(df,df_col){
  n_var = length(df_col)
  Var = colnames(df)[df_col]
  mod_name = paste0("Surv_trans~")
  for(i in 1:(n_var-1))mod_name=paste0(mod_name,Var[i],"+")
  mod_name = paste0(mod_name,Var[n_var])
  mod_name
}

gen.formula.lambda <- function(df,df_col){
  n_var = length(df_col)
  Var = colnames(df)[df_col]
  mod_name = paste0("Lambda~")
  for(i in 1:(n_var-1))mod_name=paste0(mod_name,Var[i],"+")
  mod_name = paste0(mod_name,Var[n_var])
  mod_name
}

gen.formula.N <- function(df,df_col){
  n_var = length(df_col)
  Var = colnames(df)[df_col]
  mod_name = paste0("N_log~")
  for(i in 1:(n_var-1))mod_name=paste0(mod_name,Var[i],"+")
  mod_name = paste0(mod_name,Var[n_var])
  mod_name
}



sim.beaver <- function(inputs){
  library(fishSim)
  
  #define some functions - note doing this within sim_beaver so it will work in parallel
  
  #mogogamous mating function
  monoMate <- function (indiv, batchSize, fecundityDist = "poisson", 
                        osr = c(0.5, 0.5), year = "-1", maxClutch = Inf, Partner) 
  {
    mothers <- subset(indiv, indiv[, 2] == "F" & is.na(Partner)==FALSE) #only partnered mothers mate in this version
    mPartner <- subset(Partner,indiv[,2] == "F" & is.na(Partner)==FALSE)
    getClutch <- function(fecundityDist, n, batchSize) {
      if (fecundityDist == "poisson") {
        clutch <- rpois(n = n, lambda = batchSize)
      }
      if (fecundityDist == "truncPoisson") {
        clutch <- rTruncPoisson(n = n, T = batchSize)
      }
      if (fecundityDist == "binomial") {
        clutch <- rbinom(n, 1, prob = batchSize)
      }
      if (fecundityDist == "multinomial") {
        clutch <- sample(0:(length(batchSize) - 1), size = n, 
                         replace = TRUE, prob = batchSize)
      }
      return(clutch)
    }
    clutch <- getClutch(fecundityDist = fecundityDist, n = nrow(mothers), 
                        batchSize = batchSize)
    mothers <- subset(mothers, clutch > 0)
    if(nrow(mothers)>0){  #possible there won't be any mothers with clutch>0
      dads <- mPartner[which(clutch>0)]
      clutch <- clutch[clutch > 0]
      clutch[clutch > maxClutch] <- maxClutch
      n_born = sum(clutch)
      Born <- makeFounders(pop=n_born,stocks=1)
      Born$Mum = rep(mothers$Me,times=clutch)
      Born$Dad = rep(dads,times=clutch)
      Born$BirthY = year
      Born$AgeLast=0
      indiv <- rbind(indiv, Born)
    }
    return(indiv)
  }
  
  sum.sq.lambda <- function(par,N){  #for estimating "realized" lambda for each simulation using least squares
    lambda = exp(par[1])
    N0 = exp(par[2])
    n_yrs = length(N)
    Nvec = rep(0,n_yrs)
    Nvec[1]=N0
    for(iyr in 2:n_yrs)Nvec[iyr]=Nvec[iyr-1]*lambda
    sum((Nvec-N)^2)
  }
  
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
  
  Age_props = Re(eigen(L)$vectors[,1]/sum(eigen(L)$vectors[,1]))  #stable age structure proportions
  
  archive <- make_archive()
  indiv = makeFounders(pop=N_init,osr=c(0.5,0.5),stocks=1,minAge=1,maxAge=15,survCurv=Age_props)
  Partner=rep(NA,N_init)
  flag = 0
  for(iyr in 1:n_yrs){
    if(flag==0){
      indiv <- mort(indiv = indiv, type = "age", ageMort = c(0,curMort), year = iyr,maxAge=n_ages)
      DeathIDs <- indiv$Me[-which(is.na(indiv$DeathY))]
      archive <- archive_dead(indiv = indiv, archive = archive)
      Partner = Partner[is.na(indiv$DeathY)]
      indiv <- remove_dead(indiv = indiv)
      if(iyr>1)Partner[which(Partner%in%DeathIDs)]=NA #for partners that died, replace old partner with missing indicator
      if(sum(indiv$Sex=="F" & indiv$AgeLast>=mat_age)==0 | sum(indiv$Sex=="M" & indiv$AgeLast>=mat_age)==0)flag=1 #ran out of mature breeders
      if(flag==0){
        N_ta[iyr,]=summary(factor(indiv$AgeLast,levels=as.character(1:(n_ages-1))))
        #assign missing monogamous partner for those sexually mature animals missing one
        M_no_partner = which(indiv$Sex=="M" & indiv$AgeLast>mat_age & is.na(Partner))
        F_no_partner = which(indiv$Sex=="F" & indiv$AgeLast>mat_age & is.na(Partner))
        n_gcd = min(length(M_no_partner),length(F_no_partner))
        if(n_gcd>0){
          if(length(M_no_partner)==n_gcd){
            M_partner = sample(indiv$Me[F_no_partner],n_gcd,replace=F)
          }
          else{
            n_diff = length(M_no_partner)-n_gcd
            M_partner = sample(c(indiv$Me[F_no_partner],rep(NA,n_diff)),length(M_no_partner),replace=F)
          }
          Partner[M_no_partner]=M_partner
          for(ipart in 1:length(F_no_partner)){
            which_part = which(Partner==indiv$Me[F_no_partner[ipart]])
            if(length(which_part)>0)Partner[F_no_partner[ipart]]=indiv$Me[which_part]
          }
        }
        N_breed[iyr] = sum(! is.na(Partner))
        if(N_breed[iyr]>0){
          indiv <- monoMate(indiv = indiv,osr = c(0.5,0.5),batchSize=inputs$Const$E_offspring, fecundityDist="poisson",year = iyr,Partner=Partner)
          indiv <- birthdays(indiv = indiv)
          N_ta2[iyr,]=summary(factor(indiv$AgeLast,levels=c(as.character(1:n_ages))))
          diff = nrow(indiv)-length(Partner)
          if(diff>0)Partner=c(Partner,rep(NA,diff))
        }
        else flag=1
      }
    }
  }  
  if(flag==0){
    # sample from archive (dead only)
    for(isampyr in 1:n_yrs_sample){
      cur_yr = N_yrs_sample[isampyr]
      Cur_which = which(archive$DeathY==cur_yr)
      if(length(Cur_which)>=samples_per_year)
        archive[sample(Cur_which,samples_per_year),"SampY"]=N_yrs_sample[isampyr]
      else flag=1
    }
  }
  if(flag==0){
    Sampled = archive[which(!is.na(archive$SampY)),]
    POPs <- quickin(Sampled,max_gen=2)[["POP"]]
    HSPs <- quickin(Sampled,max_gen=2)[["HSP"]]
    FSPs <- quickin(Sampled,max_gen=2)[["FSP"]]
    GGPs <- quickin(Sampled,max_gen=2)[["GGP"]]
    FTPs <- quickin(Sampled,max_gen=2)[["FTP"]]
    
    HSP_like = HSPs
    if(nrow(GGPs)>0)HSP_like = rbind(HSP_like,GGPs)
    if(nrow(FTPs)>0)HSP_like = rbind(HSP_like,FTPs)
    
    ###compile statistics pseudo-likelihood and ABC models
    # here, 'bj' indexes birth year of the younger animal, di gives death year of parent
    n_comp_PO_dibj=n_match_PO_dibj=array(0,dim=c(n_yrs,n_yrs))  #suff stat
    n_comp_PO_bidibj=n_match_PO_bidibj=n_match_PO_bidibj=array(0,dim=c(n_yrs,n_yrs,n_yrs)) #for debugging
    if(length(POPs)>0){
      if(length(POPs)==2)POPs=matrix(POPs,1,2) #make sure this is a matrix
      for(i in 1:nrow(POPs)){
        indiv1 = Sampled[which(Sampled[,1]==POPs[i,1]),]
        indiv2 = Sampled[which(Sampled[,1]==POPs[i,2]),]
        
        FIRST = (indiv1$BirthY < indiv2$BirthY)
        # case 1: parent is animal 1
        if(FIRST){
          p_by = indiv1$BirthY
          p_sy = indiv1$SampY
          o_by = indiv2$BirthY  
        }
        else{
          p_by = indiv2$BirthY
          p_sy = indiv2$SampY
          o_by = indiv1$BirthY  
        }
        if(p_sy>(o_by+2)){  #only tabulate once young are independent from parents
          n_match_PO_dibj[p_sy,o_by]=n_match_PO_dibj[p_sy,o_by]+1 
          n_match_PO_bidibj[p_by,p_sy,o_by]=n_match_PO_bidibj[p_by,p_sy,o_by]+1 
        }
      }
    }
    
    #FSP, HSP+
    n_match_FSP_delta = rep(0,15)
    n_match_FSP_bi = rep(0,30)
    if(length(FSPs)>0){
      if(length(FSPs)==2)FSPs=matrix(FSPs,1,2) #make sure this is a matrix
      for(i in 1:nrow(FSPs)){
        indiv1 = Sampled[which(Sampled[,1]==FSPs[i,1]),]
        indiv2 = Sampled[which(Sampled[,1]==FSPs[i,2]),]
        b_int = abs(indiv1$BirthY-indiv2$BirthY)
        s_int = abs(indiv1$SampY-indiv2$SampY)
        min_yr = min(indiv1$BirthY,indiv2$BirthY)
        if(b_int>2 | s_int>1){
          n_match_FSP_delta[b_int+1]=n_match_FSP_delta[b_int+1]+1
          n_match_FSP_bi[min_yr]=n_match_FSP_bi[min_yr]+1
        }
      }
    }
    n_match_HSPp_delta = rep(0,15)
    n_match_HSPp_bi = rep(0,30)
    if(length(HSP_like)>0){
      if(length(HSP_like)==2)HSP_like=matrix(HSP_like,1,2) #make sure this is a matrix
      for(i in 1:nrow(HSP_like)){
        indiv1 = Sampled[which(Sampled[,1]==HSP_like[i,1]),]
        indiv2 = Sampled[which(Sampled[,1]==HSP_like[i,2]),]
        b_int = abs(indiv1$BirthY-indiv2$BirthY)
        s_int = abs(indiv1$SampY-indiv2$SampY)
        min_yr = min(indiv1$BirthY,indiv2$BirthY)
        if(b_int>2 | s_int>1){
          n_match_HSPp_delta[b_int+1]=n_match_HSPp_delta[b_int+1]+1
          n_match_HSPp_bi[min_yr]=n_match_HSPp_bi[min_yr]+1
        }
      }
    }
    
    lam_mod <- nlminb(c(0,log(100)),sum.sq.lambda,N=N_breed)
    lambda = exp(lam_mod$par[1])
    n_match_PO_bj = colSums(n_match_PO_dibj)
    
  }
  
  if(flag==1)return(NULL)
  else return(list("n_match_PO_dibj"=n_match_PO_dibj,"n_comp_PO_dibj"=n_comp_PO_dibj,"lambda"=lambda,"N_mean"=mean(N_breed),
                   "n_match_FSP_delta"=n_match_FSP_delta,"n_match_FSP_bi"=n_match_FSP_bi,
                   "n_match_HSPp_delta"=n_match_HSPp_delta,"n_match_HSPp_bi"=n_match_HSPp_bi,"N_breed"=N_breed,"N_ta"=N_ta,
                   "n_match_PO_bj"=n_match_PO_bj))
  
}

##############
# Main sim loop
##############
n_sim = 100
cur_seed = 10100
n_ABC_sims = 24000  
n_clusters=12
cl <- makeCluster(n_clusters)
abc_quantile = 0.02

S_ad = 0.77 #adult survival (age 3+)
S_juv = 0.65 #juvenile survival (1-2)
S_0 = 0.38 #kit survival (0)
N0 = 200  #first year abundance
n_yrs = 30
n_ages=15
N_yrs_sample = c(16:30)
n_yrs_sample = length(N_yrs_sample)
samples_per_year = 20
N_ta = matrix(0,n_yrs,n_ages-1)
N_ta2 = matrix(0,n_yrs,n_ages)
Repro = c(0,0,0,rep(4,12)) 
mat_age = 3
Yrs = c(1:n_yrs_sample)
ageMort = 1-c(S_0,S_juv,S_juv,rep(S_ad,11),0)

ABC_data = vector("list",n_sim) #so we can rerun ABC with different algorithms post hoc
Lambda_ABC = N_breed_ABC = rep(0,n_ABC_sims)

Kin_matches = matrix(NA,n_sim,3)

L = matrix(0,15,15)
L[1,] = c(0,0,0,rep(2*S_ad,12))
L[2,1]=S_0
L[3,2]=L[4,3]=S_juv
for(i in 5:15)L[i,i-1]=S_ad

eigen(L)

Const = list("n_yrs"=n_yrs,"n_ages"=n_ages,"N_yrs_sample"=N_yrs_sample,
             "n_yrs_sample"=n_yrs_sample,"N_breed"=rep(0,n_yrs),
             "samples_per_year"=samples_per_year,"N_ta"=N_ta,"N_ta2"=N_ta2,
             "Repro"=Repro,"mat_age"=mat_age,"Yrs"=Yrs,"E_offspring"=4)

inputs = list("Const"=Const)
inputs$N_init=N0
inputs$curMort = ageMort
inputs$L = L
Truth=inputs

Nbreed_true = Lambda_realized = rep(NA,n_sim)
Cover = Bias = MSE = array(NA,dim=c(n_sim,3,2)) 
# 3,2 dimensions are parameters (lambda, N_breed, phiA) and ABC estimator (NN, semi-auto)

st_time <- Sys.time()

inputs = list("Const"=Const)
reduced_lm = vector("list",3)  #holds reduced linear model for previous simulation (used for caculating summary statistics for semi-auto ABC)

for(isim in 1:n_sim){
  set.seed(cur_seed)
  cur_seed = cur_seed+1
  inputs$N_init=Truth$N_init
  inputs$curMort = Truth$curMort
  inputs$L = Truth$L
  
  sim_data <- sim.beaver(inputs)
  
  if(is.null(sim_data)==FALSE){
    Nbreed_true[isim] = mean(sim_data$N_breed)
    lam_mod <- nlminb(c(0,log(100)),sum.sq.lambda,N=sim_data$N_breed)
    Lambda_realized[isim] = exp(lam_mod$par[1])
    
    
    Kin_matches[isim,]=c(sum(sim_data$n_match_PO_bj),sum(sim_data$n_match_FSP_bi),
                         sum(sim_data$n_match_HSPp_bi))
    
    Lims_N0 = c(50,500)
    Lims_phi0 = c(0.2,0.55)
    Lims_phiA = c(0.6,0.94)
    
    ABC_inputs_list = vector("list",n_ABC_sims)
    ABC_pars = matrix(0,n_ABC_sims,3)
    ABC_pars[,1]=runif(n_ABC_sims,Lims_N0[1],Lims_N0[2])
    ABC_pars[,2]=runif(n_ABC_sims,Lims_phi0[1],Lims_phi0[2])
    ABC_pars[,3]=runif(n_ABC_sims,Lims_phiA[1],Lims_phiA[2])
    
    for(i in 1:n_ABC_sims){
      L[2,1]=ABC_pars[i,2]
      L[1,4:15]=2*ABC_pars[i,3]
      for(iage in 4:14)L[iage+1,iage]=ABC_pars[i,3]
      Lambda_ABC[i]=Re(eigen(L)$values[1])
      while(Lambda_ABC[i]>1.1 | Lambda_ABC[i]<0.9){ #try to keep populations from going too crazy
        ABC_pars[i,2]=runif(1,Lims_phi0[1],Lims_phi0[2])
        ABC_pars[i,3]=runif(1,Lims_phiA[1],Lims_phiA[2])
        L[2,1]=ABC_pars[i,2]
        L[1,4:15]=2*ABC_pars[i,3]
        for(iage in 4:14)L[iage+1,iage]=ABC_pars[i,3]
        Lambda_ABC[i]=Re(eigen(L)$values[1])
      }
      ABC_inputs_list[[i]]=list(N_init=ABC_pars[i,1],
                                curMort=c(1-ABC_pars[i,2],ageMort[c(2,3)],rep(1-ABC_pars[i,3],11),1),
                                L=L,
                                Const=Const)
    }
    
    sim_out = parSapply(cl,ABC_inputs_list,sim.beaver)  
    
    Stat_mat = matrix(0,n_ABC_sims,17+11+16+10+16)   
    No_est = rep(0,n_ABC_sims)
    
    for(i in 1:n_ABC_sims){
      if(is.null(sim_out[[i]])==FALSE){
        N_breed_ABC[i]=sim_out[[i]]$N_mean
        Stat_mat[i,]=c(sim_out[[i]]$n_match_PO_bj[11:27],
                       sim_out[[i]]$n_match_FSP_delta[1:11],
                       sim_out[[i]]$n_match_FSP_bi[11:26],
                       sim_out[[i]]$n_match_HSPp_delta[2:11],
                       sim_out[[i]]$n_match_HSPp_bi[11:26]
        )
      }
      else No_est[i]=1
    }
    
    Which_no_est = which(No_est==1)
    Par_dm = data.frame("Lambda"=Re(Lambda_ABC),"Nbreed"=N_breed_ABC,"PhiA"=ABC_pars[,3])
    Par_dm = Par_dm[-Which_no_est,]
    Stat_mat= Stat_mat[-Which_no_est,]
    
    #save some output
    ABC_data[[isim]]$target = c(sim_data$n_match_PO_bj[11:27],
                                sim_data$n_match_FSP_delta[1:11],
                                sim_data$n_match_FSP_bi[11:26],
                                sim_data$n_match_HSPp_delta[2:11],
                                sim_data$n_match_HSPp_bi[11:26])
    ABC_data[[isim]]$param = Par_dm
    ABC_data[[isim]]$sumstat = Stat_mat
    
    logit_bounds = matrix(0,3,2)
    logit_bounds[3,2]=1
    abc_out <- abc(target=c(sim_data$n_match_PO_bj[11:27],
                            sim_data$n_match_FSP_delta[1:11],
                            sim_data$n_match_FSP_bi[11:26],
                            sim_data$n_match_HSPp_delta[2:11],
                            sim_data$n_match_HSPp_bi[11:26]),
                   param=Par_dm,tol=abc_quantile,
                   sumstat=Stat_mat,transf=c("log","log","logit"),
                   logit.bounds=logit_bounds,
                   method="neuralnet",sizenet=5)
    summary_abc = summary(abc_out,print=FALSE,intvl=0.9)
    
    Cover[isim,1,1]= (Lambda_realized[isim]>summary_abc[2,1] & Lambda_realized[isim]<summary_abc[6,1])
    Cover[isim,2,1]= (Nbreed_true[isim]>summary_abc[2,2] & Nbreed_true[isim]<summary_abc[6,2])
    Cover[isim,3,1]= (0.77>summary_abc[2,3] & 0.77<summary_abc[6,3])
    Bias[isim,,1]= (summary_abc[5,]-c(Lambda_realized[isim],Nbreed_true[isim],0.77))/c(Lambda_realized[isim],Nbreed_true[isim],0.77)
    MSE[isim,,1]= (summary_abc[5,]-c(Lambda_realized[isim],Nbreed_true[isim],0.77))^2
    
    #### Semi-automated ABC
    # fit regression models for each parameter
    last_lm = reduced_lm

    # Lambda
    Stats = data.frame(cbind(ABC_data[[isim]]$sumstat[,1:17],ABC_data[[isim]]$sumstat[,29:44],ABC_data[[isim]]$sumstat[,55:70]))
    Stats = cbind(Stats,Stats^2)
    colnames(Stats)=c(paste0("POPlin",c(1:17)),paste0("FSPlin",c(1:16)),paste0("HSPlin",c(1:16)),
                      paste0("POPquad",c(1:17)),paste0("FSPquad",c(1:16)),paste0("HSPquad",c(1:16)))
    Which_zero = which(colSums(Stats)==0)
    Stats$Lambda = Par_dm[,1]
    Which_cols_model = c(1:(ncol(Stats)-1))
    if(length(Which_zero)>0)Which_cols_model=Which_cols_model[-Which_zero]
    global_model = gen.formula.lambda(Stats,Which_cols_model)
    global_lm = lm(global_model,data=Stats,na.action = "na.fail")
    pvals = summary(global_lm)$coefficients[,4]
    cols_keep = names(which(pvals[2:length(pvals)]<0.05))
    which_keep = which(colnames(Stats) %in% cols_keep)
    red_formula = gen.formula.lambda(Stats,which_keep)
    reduced_lm[[1]] = lm(red_formula,data=Stats,na.action="na.fail")
    
    # Mean number of breeders
    Stats = data.frame(cbind(ABC_data[[isim]]$sumstat[,1:17],ABC_data[[isim]]$sumstat[,29:44],ABC_data[[isim]]$sumstat[,55:70]))
    Stats = cbind(Stats,Stats^2)
    colnames(Stats)=c(paste0("POPlin",c(1:17)),paste0("FSPlin",c(1:16)),paste0("HSPlin",c(1:16)),
                      paste0("POPquad",c(1:17)),paste0("FSPquad",c(1:16)),paste0("HSPquad",c(1:16)))
    Which_zero = which(colSums(Stats)==0)
    Stats$Nbreed = Par_dm[,2]
    Stats$N_log = log(Stats$Nbreed)
    Which_cols_model = c(1:(ncol(Stats)-2))
    if(length(Which_zero)>0)Which_cols_model=Which_cols_model[-Which_zero]
    global_model = gen.formula.N(Stats,Which_cols_model)
    global_lm = lm(global_model,data=Stats,na.action = "na.fail")
    pvals = summary(global_lm)$coefficients[,4]
    cols_keep = names(which(pvals[2:length(pvals)]<0.05))
    which_keep = which(colnames(Stats) %in% cols_keep)
    red_formula = gen.formula.N(Stats,which_keep)
    reduced_lm[[2]] = lm(red_formula,data=Stats,na.action="na.fail")
    
    # Adult survival
    Stats_S = data.frame(cbind(ABC_data[[isim]]$sumstat[,18:28],ABC_data[[isim]]$sumstat[,45:54]))
    Stats = cbind(Stats_S,Stats_S^2)
    colnames(Stats)=c(paste0("FSPlin",c(1:11)),paste0("HSPlin",c(1:10)),
                      paste0("FSPquad",c(1:11)),paste0("HSPquad",c(1:10)))
    Stats$Surv = Par_dm[,3]
    Stats$Surv_trans = logit(Stats$Surv)
    global_model = gen.formula.S(Stats,c(1:(ncol(Stats)-2)))
    global_lm = lm(global_model,data=Stats,na.action = "na.fail")
    pvals = summary(global_lm)$coefficients[,4]
    which_keep = which(pvals[2:length(pvals)]<0.05)
    red_formula = gen.formula.S(Stats,which_keep)
    reduced_lm[[3]] = lm(red_formula,data=Stats,na.action="na.fail")
    
    if(isim>1){ #now the ABC part - this will use the models fitted in the previous simulation
                #so as not to "double dip"
      
      #### now, calc ssq for real versus abc parameter predictions, quantile, etc.
      
      #real data
      Real_df = data.frame(matrix(c(ABC_data[[isim]]$target[c(1:17,29:44,55:70)],
                                    (ABC_data[[isim]]$target[c(1:17,29:44,55:70)])^2),nrow=1))
      colnames(Real_df)=c(paste0("POPlin",c(1:17)),paste0("FSPlin",c(1:16)),paste0("HSPlin",c(1:16)),
                          paste0("POPquad",c(1:17)),paste0("FSPquad",c(1:16)),paste0("HSPquad",c(1:16)))
      predict_real_lambda = predict(last_lm[[1]],newdata=Real_df)
      predict_real_N = exp(predict(last_lm[[2]],newdata=Real_df))      
      
      Real_df = data.frame(matrix(c(ABC_data[[isim]]$target[c(18:28,45:54)],
                                    (ABC_data[[isim]]$target[c(18:28,45:54)])^2),nrow=1))
      colnames(Real_df)=c(paste0("FSPlin",c(1:11)),paste0("HSPlin",c(1:10)),
                          paste0("FSPquad",c(1:11)),paste0("HSPquad",c(1:10)))
      predict_real_S = plogis(predict(last_lm[[3]],newdata=Real_df))
      
      #abc replicates
      n_conv_ABC = nrow(ABC_data[[isim]]$sumstat)
      ABC_preds = weighted_SSQ = matrix(NA,n_conv_ABC,3)
      for(iabc in 1:n_conv_ABC){
        ABC_df =  data.frame(matrix(c(ABC_data[[isim]]$sumstat[iabc,c(1:17,29:44,55:70)],
                                    (ABC_data[[isim]]$sumstat[iabc,c(1:17,29:44,55:70)])^2),nrow=1))
        colnames(ABC_df)=c(paste0("POPlin",c(1:17)),paste0("FSPlin",c(1:16)),paste0("HSPlin",c(1:16)),
                            paste0("POPquad",c(1:17)),paste0("FSPquad",c(1:16)),paste0("HSPquad",c(1:16)))
        ABC_preds[iabc,1] = predict(last_lm[[1]],newdata=ABC_df)
        ABC_preds[iabc,2] = exp(predict(last_lm[[2]],newdata=ABC_df))  
        
        ABC_df = data.frame(matrix(c(ABC_data[[isim]]$sumstat[iabc,c(18:28,45:54)],
                                      (ABC_data[[isim]]$sumstat[iabc,c(18:28,45:54)])^2),nrow=1))
        colnames(ABC_df)=c(paste0("FSPlin",c(1:11)),paste0("HSPlin",c(1:10)),
                            paste0("FSPquad",c(1:11)),paste0("HSPquad",c(1:10)))
        ABC_preds[iabc,3] = plogis(predict(last_lm[[3]],newdata=ABC_df))
      }

      weighted_SSQ[,1]=(predict_real_lambda - ABC_preds[,1])^2/var(ABC_preds[,1])
      weighted_SSQ[,2]=(predict_real_N - ABC_preds[,2])^2/var(ABC_preds[,2])
      weighted_SSQ[,3]=(predict_real_S - ABC_preds[,3])^2/var(ABC_preds[,3])
      
      Compound_SSQ = rowSums(weighted_SSQ)
      quant_ssq = quantile(Compound_SSQ,abc_quantile)
      
      Ests_FP = ABC_data[[isim]]$param[which(Compound_SSQ<=quant_ssq),]

      Cover[isim,1,2]= (Lambda_realized[isim]>quantile(Ests_FP[,1],0.05) & Lambda_realized[isim]<quantile(Ests_FP[,1],0.95))
      Cover[isim,2,2]= (Nbreed_true[isim]>quantile(Ests_FP[,2],0.05) & Nbreed_true[isim]<quantile(Ests_FP[,2],0.95))
      Cover[isim,3,2]= (0.77>quantile(Ests_FP[,3],0.05) & 0.77<quantile(Ests_FP[,3],0.95))
      #kde to get mode to use as point estimate
      Post_modes = rep(0,3)
      for(ipar in 1:3){
        kde = density(Ests_FP[,ipar])
        Post_modes[ipar]=kde$x[which(kde$y==max(kde$y))]
      }
      Bias[isim,,2]= (Post_modes-c(Lambda_realized[isim],Nbreed_true[isim],0.77))/c(Lambda_realized[isim],Nbreed_true[isim],0.77)
      MSE[isim,,2]= (Post_modes-c(Lambda_realized[isim],Nbreed_true[isim],0.77))^2
    }
  }
  if((isim %% 1)==0)save.image("Sim_mono_24000_2.RData")
}
parallel::stopCluster(cl)

Sys.time()-st_time

#Performance table
Bias_mean = apply(Bias,c(2,3),'mean',na.rm=TRUE)
Cover_mean = apply(Cover,c(2,3),'mean',na.rm=TRUE)
RMSE_mean = sqrt(apply(MSE,c(2,3),'mean',na.rm=TRUE))
Perf_table = data.frame("Parameter" = c(rep("lambda",2),rep("Nf",2),rep("Sa",2)),
                        "Estimator" = rep(c("nn","FP"),3),
                        "Bias" = as.numeric(t(Bias_mean)),
                        "Cover" = as.numeric(t(Cover_mean)),
                        "RMSE" = as.numeric(t(RMSE_mean)))
library(xtable)
xtable(Perf_table)



