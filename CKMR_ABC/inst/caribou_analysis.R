# caribou analysis

caribou_data = read.csv("./data/caribou.csv",header=T)

#form Leslie matrix to approximate expected fraction of 
# non-calves that are actually yearlings

S_c= .15
S_y= S_s = 0.8
S_a= 0.85
f_y = 0
f_s = .15
f_a = 0.9

n_ages=20
L = matrix(0,n_ages,n_ages)
L[1,2]=0.5*f_y*S_y    #expected # of female calves for yearlings
L[1,3]=0.5*f_s*S_s  # 2 year olds
L[1,4:n_ages]=0.5*f_a*S_a  #expected # of female calves for adults
L[2,1]=S_c
L[3,2]=S_y
L[4,3]=S_s
for(i in 5:n_ages)L[i,i-1]=S_a

eigen(L)$values[1]
stable_stage=eigen(L)$vector[,1]
stable_stage = stable_stage/sum(stable_stage) #standardize
stable_stage

prop_y_repro = as.numeric(stable_stage[2]/sum(stable_stage[2:n_ages]))
prop_s_repro = as.numeric(stable_stage[3]/sum(stable_stage[2:n_ages]))
prop_a_repro = 1-prop_y_repro-prop_s_repro

# proportions of ages shift slightly at time of pellet counts due to differences in 8-month survival
denom = sum(prop_y_repro*S_y^0.666+prop_s_repro*S_s^0.666+prop_a_repro*S_a^0.666)
prop_y_count = prop_y_repro*S_y^0.666 / denom
prop_s_count = prop_s_repro*S_s^0.666 / denom
prop_a_count = 1-prop_y_count-prop_s_count
prop_breed = prop_a_count+prop_s_count*0.15


chapman_est <- function(n_0,n_f,x){
  floor((n_0+1)*(n_f+1)/(x+1)-1)
}

chapman_CI <- function(n_0,n_f,x){
  sigma_hat = sqrt(1/(x+0.5)+1/(n_f-x+0.5)+1/(n_0-x+0.5)+(x+0.5)/((n_0-x+0.5)*(n_f-x+0.5)))
  n_f+n_0-x-0.5+(n_f-x+0.5)*(n_0-x+0.5)/(x+0.5)*exp(c(-1,1)*1.645*sigma_hat)  
}

n_yrs = nrow(caribou_data)
Chapman_ests = matrix(0,n_yrs,3)
for(i in 1:n_yrs){
  Chapman_ests[i,1]=chapman_est(n_0=caribou_data[i,"n_calves"],
                      n_f = caribou_data[i,"n_pot_mo"]*prop_breed,
                      x=caribou_data[i,"n_MOP"])
  Chapman_ests[i,2:3]=chapman_CI(n_0=caribou_data[i,"n_calves"],
                        n_f = caribou_data[i,"n_pot_mo"]*prop_breed,
                        x=caribou_data[i,"n_MOP"])
}

library(ggplot2)

point_ests = data.frame("Year"=as.character(rep(c(2006:2015),3)),
                        "Estimator"=c(rep("CMR",n_yrs),rep("CKMR-pseudo",n_yrs),rep("CKMR-chapman",n_yrs)),
                        "N"=c(caribou_data$CMR_F,caribou_data$CKMR_NF,Chapman_ests[,1]),
                        "Lower"=c(caribou_data$CMR_lower,caribou_data$CKMR_Lower,Chapman_ests[,2]),
                        "Upper"=c(caribou_data$CMR_upper,caribou_data$CKMR_upper,Chapman_ests[,3]))
                        
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot() + geom_errorbar(data=point_ests,position=position_dodge(width=0.5),
                         aes(x=Year,y=N,ymin=Lower,ymax=Upper,group=Estimator,color=Estimator),size=0.8,width=0.3)+
                          coord_cartesian(ylim=c(0, 100))+
           geom_point(data=point_ests,position=position_dodge(width=0.8),aes(x=Year,y=N,color=Estimator,group=Estimator,shape=Estimator),size=1.5)+
           scale_color_manual(values=cbp2)+
           ylab(expression(N[f]))


#now the ABC part

#first, write population simulator and virtual data collector function
sim.caribou <- function(inputs){
  library(fishSim)
  N_init=inputs$N_init
  curMort=inputs$curMort
  phi0=1-inputs$curMort[1]
  L = inputs$L
  n_yrs=inputs$Const$n_yrs
  n_ages=inputs$Const$n_ages
  m_samples_per_year=inputs$Const$m_samples_per_year
  c_samples_per_year=inputs$Const$c_samples_per_year
  N_ta=inputs$Const$N_ta
  N_ta2=inputs$Const$N_ta2
  mat_age=inputs$Const$mat_age
  Yrs=inputs$Const$Yrs
  
  curMort_1 = 1-(1-curMort)^0.666  #8 month mortality probabilities
  curMort_2 = 1-(1-curMort)^0.333  #4 month mortality probabilities
  
  Age_props = Re(eigen(L)$vectors[,1]/sum(eigen(L)$vectors[,1]))  #stable age structure proportions
  
  Matches = rep(0,n_yrs)
  archive <- make_archive()
  indiv = makeFounders(pop=N_init,osr=c(0.5,0.5),stocks=1,minAge=1,maxAge=20,survCurv=Age_props)
  flag = 0
  for(iyr in 1:(n_yrs+1)){  #have to run an additional year so calves all have 'real' moms
    if(flag==0){
      #survival to pellet collection
      which_F = which(indiv$Sex=="F")
      if(iyr>1)N_ta[iyr-1,]=tabulate(indiv$AgeLast[which_F],nbins=20)      
      indiv <- mort(indiv = indiv, type = "age", ageMort = c(0,curMort_1), year = iyr,maxAge=n_ages)
      indiv <- remove_dead(indiv = indiv)
      
      if(sum(indiv$Sex=="F" & indiv$AgeLast>2)==0 | sum(indiv$Sex=="M" & indiv$AgeLast>2)==0)flag=1 #ran out of mature breeders
      if(flag==0){
        which_F = which(indiv$Sex=="F")
        if(iyr>1)N_ta2[iyr-1,]=tabulate(indiv$AgeLast[which_F],nbins=20)
        
        #pellet collection
        if(iyr>1){
          which_calves = which(indiv$AgeLast==1)
          n_calves=length(which_calves)
          if(n_calves>c_samples_per_year[iyr-1]){
            Moms_of_calves = indiv[sample(which_calves,c_samples_per_year[iyr-1],replace=F),"Mum"]
          }
          else flag=1
          
          if(flag==0){
            which_moms = which(indiv$AgeLast>1 & indiv$Sex=="F")
            n_moms = length(which_moms)
            if(n_moms>m_samples_per_year[iyr-1]){
              Moms = indiv[sample(which_moms,m_samples_per_year[iyr-1]),"Me"]  
            }
            else flag=1
            
            if(flag==0)Matches[iyr-1] = sum(Moms_of_calves %in% Moms)
          }
        }
        
        #breeding
        temp = altMate(indiv,batchSize=1,fecundityDist="binomial",maturityCurve=c(0,mat_age),year=iyr,type="age")
        which_new = which(temp$BirthY==iyr)
        if(length(which_new)>0){
          potential_births = temp[which_new,]
        }
        dim(potential_births)
        
        #survival to end of year
        indiv <- mort(indiv = indiv, type = "age", ageMort = c(0,curMort_2), year = iyr,maxAge=n_ages)

        #only retain calves if mother survived preceding interval
        which_dead_moms = which(indiv$DeathY==iyr & indiv$Sex=="F")
        if(length(which_dead_moms)>0){
          deadIDs = indiv[which_dead_moms,"Me"]
          which_fetus_die = which(potential_births$Mum %in% deadIDs)
          if(length(which_fetus_die)>0)potential_births=potential_births[-which_fetus_die,]
          indiv=rbind(indiv,potential_births)
        }
        indiv <- remove_dead(indiv = indiv)
        
        #increment ages
        if(sum(indiv$Sex=="F" & indiv$AgeLast>2)==0 | sum(indiv$Sex=="M" & indiv$AgeLast>2)==0)flag=1 #ran out of mature breeders
        indiv <- birthdays(indiv = indiv)
      }
    }
  }  
  if(flag==1)return(NULL)
  else return(list("N_ta2"=N_ta2,"Matches"=Matches)) 
}

mat_age = c(0,0,0.15,rep(0.9,n_ages-3))
N_ta = matrix(0,n_yrs,n_ages)
N_ta2 = matrix(0,n_yrs,n_ages)
Const = list("n_yrs"=n_yrs,"n_ages"=n_ages,
             "m_samples_per_year"=caribou_data$n_pot_mo,
             "c_samples_per_year"=caribou_data$n_calves,
             "N_ta"=N_ta,"N_ta2"=N_ta2,
             "mat_age"=mat_age,"Yrs"=c(1:n_yrs))

inputs = list("Const"=Const)

N0 = 200
inputs$N_init=N0
ageMort = 1-c(S_c,S_y,S_s,rep(S_a,n_ages-4),0)  #killing off last age group
inputs$curMort = ageMort
inputs$L = L

sim_out <-sim.caribou(inputs)
sim_out

library(parallel)
set.seed(10000)
n_clusters=12
cl <- makeCluster(n_clusters)

n_abc = 1200000
ABC_inputs_list = vector("list",n_abc)
abc_quantile = 0.002
N0_abc = runif(n_abc,50,400)
S0_abc = runif(n_abc,0.05,0.65)
Lambda_abc = rep(NA,n_abc)
N_abc = matrix(NA,n_abc,10)  #includes calves
N_abc2 =  N_abc  #not including calves
Matches_abc = N_abc #holds number of matches
for(iabc in 1:n_abc){
  L[2,1]=S0_abc[iabc] 
  Lambda_abc[iabc]=as.numeric(eigen(L)$values[1])
  ABC_inputs_list[[iabc]]=list(N_init=N0_abc[iabc],
                              curMort=1-c(S0_abc[iabc],S_y,S_s,rep(S_a,n_ages-4),0),
                              L=L,
                              Const=Const)
}

st_time = Sys.time()
sim_out = parSapply(cl,ABC_inputs_list,sim.caribou)  
parallel::stopCluster(cl)
Sys.time()-st_time

# calculate some features of converging abc replicates 
for(iabc in 1:n_abc){
  if(!is.null(sim_out[[iabc]])){
    Matches_abc[iabc,]=sim_out[[iabc]]$Matches
    N_abc[iabc,]=rowSums(sim_out[[iabc]]$N_ta2)
    N_abc2[iabc,]=rowSums(sim_out[[iabc]]$N_ta2[,2:20])
  }
}


which_conv = which(!is.na(N_abc[,1]))
n_conv=length(which_conv)
N_abc_red = N_abc[which_conv,]
N_abc2_red = N_abc2[which_conv,]
Lambda_abc_red = Lambda_abc[which_conv]
Matches_abc_red = Matches_abc[which_conv,]
N0_abc_red = N0_abc[which_conv]
S0_abc_red = S0_abc[which_conv] 

lm_df = data.frame("Lambda"=Lambda_abc_red,"N0_log"=log(N0_abc_red),
                   "S0"=S0_abc_red)
colnames1 = colnames(lm_df)
lm_df = cbind(lm_df,Matches_abc_red)
lm_df = cbind(lm_df,Matches_abc_red^2)
var_names = c(paste0("lin",c(1:10)),paste0("quad",c(1:10)))
colnames(lm_df)=c(colnames1,var_names)

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
  mod_name = paste0("N0_log~")
  for(i in 1:(n_var-1))mod_name=paste0(mod_name,Var[i],"+")
  mod_name = paste0(mod_name,Var[n_var])
  mod_name
}

formula_lambda = gen.formula.lambda(lm_df,c(4:ncol(lm_df)))
formula_logN0 = gen.formula.N(lm_df,c(4:ncol(lm_df)))
lm_lambda = lm(formula_lambda,data=lm_df)
lm_logN0 = lm(formula_logN0,data=lm_df)


#real data
Real_df = data.frame(matrix(c(caribou_data[,"n_MOP"],
                              (caribou_data[,"n_MOP"])^2),nrow=1))
colnames(Real_df)=c(paste0("lin",c(1:10)),paste0("quad",c(1:10)))
predict_real_lambda = predict(lm_lambda,newdata=Real_df)
predict_real_N0 = exp(predict(lm_logN0,newdata=Real_df))      

#abc replicates
pred_abc_lambda=predict(lm_lambda)
pred_abc_N0 = exp(predict(lm_logN0))
weighted_SSQ = matrix(NA,n_conv,2)
weighted_SSQ[,1]=(predict_real_lambda - pred_abc_lambda)^2/var(pred_abc_lambda)
weighted_SSQ[,2]=(predict_real_N0 - pred_abc_N0)^2/var(pred_abc_N0)

Compound_SSQ = rowSums(weighted_SSQ)
quant_ssq = quantile(Compound_SSQ,abc_quantile)

which_select = which(Compound_SSQ<=quant_ssq)
Lambda_ests = Lambda_abc_red[which_select]
N0_ests = N0_abc_red[which_select]
S0_ests = S0_abc_red[which_select] 
N_est = N_abc_red[which_select,]
N_est2 = N_abc2_red[which_select,]

kde = density(Lambda_ests)
post_mode_lambda = kde$x[which(kde$y==max(kde$y))]
kde = density(N0_ests)
post_mode_N0 = kde$x[which(kde$y==max(kde$y))]
kde = density(S0_ests)
post_mode_S0 = kde$x[which(kde$y==max(kde$y))]

N_mode = N_mode2 = rep(0,n_yrs)
Quantile95 = matrix(0,2,n_yrs)
for(i in 1:n_yrs){
  kde=density(N_est[,i])
  N_mode[i]=kde$x[which(kde$y==max(kde$y))]
  kde=density(N_est2[,i])
  N_mode2[i]=kde$x[which(kde$y==max(kde$y))]
  Quantile95[1,i]=quantile(N_est2[,i],0.025)
  Quantile95[2,i]=quantile(N_est2[,i],0.975)
}

save.image('caribou_results.RData')


#plot these results on top of point estimates

line_ests <- data.frame("Year"=c(1:10),"N"=N_mode2,
                      "N_min"=Quantile95[1,],"N_max"=Quantile95[2,])



Big_plot = ggplot() + geom_errorbar(data=point_ests,position=position_dodge(width=0.5),
                         aes(x=Year,y=N,ymin=Lower,ymax=Upper,group=Estimator,color=Estimator),size=0.8,width=0.3)+
  coord_cartesian(ylim=c(0, 100))+
  geom_point(data=point_ests,position=position_dodge(width=0.5),aes(x=Year,y=N,color=Estimator,group=Estimator,shape=Estimator),size=1.5)+
  scale_color_manual(values=cbp2)+
  ylab(expression(N[f]))+
  geom_line(data=line_ests,aes(x=Year,y=N))+
  geom_ribbon(data=line_ests,aes(x=Year,ymin = N_min, ymax = N_max), fill = "grey70",alpha=0.5)

png("Caribou_ests.png",width=6,height=6,units="in",res=800)
  Big_plot
dev.off()


#prior posterior plot
n_abc_red = length(Lambda_abc_red)
n_ests=length(Lambda_ests)
Prior_df = data.frame("value"=c(S0_abc,N0_abc),
                      "Parameter"=c(rep("S[0]",n_abc),
                                    rep("N[0]",n_abc)),
                      "Type" = "Prior")
Prepost_df = data.frame("value"=c(S0_abc_red,N0_abc_red),
                        "Parameter"=c(rep("S[0]",n_abc_red),
                                      rep("N[0]",n_abc_red)),
                        "Type" = "Non-censored")
Posterior_df = data.frame("value"=c(S0_ests,N0_ests),
                        "Parameter"=c(rep("S[0]",n_ests),
                                      rep("N[0]",n_ests)),
                        "Type" = "Posterior")
PP_df = rbind(Prior_df,Prepost_df,Posterior_df)
PP_df$Type = factor(PP_df$Type,levels=c("Prior","Non-censored","Posterior"))



PP_plot = ggplot(PP_df)+geom_density(aes(x=value,group=Type,fill=Type),alpha=0.5)+
  facet_wrap(~Parameter,nrow=1,scales="free",labeller = label_parsed)+
  xlab("Parameter value")+ylab("Density")


png("Caribou_pp.png",width=6,height=6,units="in",res=800)
PP_plot
dev.off()








  
  