library(ggplot2)
library(NMOF)
library(mvtnorm)
library(TwoSampleMR)
library(reshape2)
library(dplyr)
library(stringr)
library(ggpubr)
library(extrafont)
library(SimDesign)
library(GGally)
library(mixtools)
library(data.table)
library(parallel)
library(tidyverse)
library(extRemes)

comp_run=TRUE

grand_dir 
script_dir 
setwd(grand_dir)

input_X <- "RLS_pooled_MR.txt"
input_Y <- "QSM_Right_caudate_MR.txt"
fX = fread(paste0(grand_dir,input_X),header=TRUE,stringsAsFactors=F) ### Exposure file if non-UKBB (overlapping SNPs between EXP-OUT)
fY = fread(paste0(grand_dir,input_Y),header=TRUE,stringsAsFactors=F) ### Outcome file if non-UKBB (overlapping SNPs between EXP-OUT)

LDfile = fread("LHC-MR/RLS_IRON_RUN/data/LDscores.txt", sep="\t", header=TRUE)

attach(LDfile)
mafT = .005
selF = which(info>.99 & abs(mafUK10K-mafUKBB)<(2*sqrt(1/4e3+1/360e3)) & mafUKBB>mafT & mafUK10K>mafT & 
               !(chr==6 & pos>=28.5e6 & pos<=33.5e6))
LDfile = LDfile[selF,]

colnames(fX)[2] <- colnames(fY)[2] <- "rsid"
colnames(fX)[9] <- colnames(fY)[9] <- "n_complete_samples"

snps <- read.table("w_hm3.noMHC.snplist",header=T,stringsAsFactors = F)

snp_tokeep = intersect(intersect(LDfile$rs, fX$rsid),intersect(fY$rsid,snps$SNP))
ld_wd = LDfile[match(snp_tokeep,LDfile$rs),]
fX = fX[match(snp_tokeep,fX$rsid),]	
fY = fY[match(snp_tokeep,fY$rsid),]

all(fX$rsid==ld_wd$rs)	
all(fX$rsid==fY$rsid)

nX = mean(fX$n_complete_samples)  #Get sample size for trait X
nY = mean(fY$n_complete_samples)  #Get sample size for trait Y

bX = fX$beta   #Get standardised beta for trait X
bY = fY$beta   #Get standardised beta for trait Y

downsampleXN <- 100
bX = bX[seq(1, length(bX), downsampleXN)]
bY = bY[seq(1, length(bY), downsampleXN)]
ld = as.numeric(unlist(ld_wd[seq(1, nrow(ld_wd), downsampleXN),'LDSC']))
wd = as.numeric(unlist(ld_wd[seq(1, nrow(ld_wd), downsampleXN),'weight']))

M = length(bX)  #Number of SNPs
M

source_script_dir = "LHC-MR/RLS_IRON_RUN/data/Scripts"
model_script_map <- c(
  "full" = "optim_comp.R",
  "only_latent" = "optim_ab.R",
  "only_causal" = "optim_U.R",
  "only_a" = "optim_bU.R",
  "only_b" = "optim_aU.R",
  "no_a" = "optim_a.R",
  "no_b" = "optim_b.R",
  "NULL"= "optim_abU.R"
)
running_dir = "LHC-MR/RLS_IRON_RUN/"
model_names = c("full","only_latent","only_causal","only_a","only_b","no_a","no_b","NULL")
mLL_model_res_list <- lapply(model_names,function(model_name_i)
                               mLL_model_f(model_name=model_name_i,running_dir=running_dir,source_script_dir=source_script_dir,mc.cores=12,maxit=50,model_script_map=model_script_map,phen_corr=0,nX,nY,bX,bY,ld,wd,M)
)
names (mLL_model_res_list) <- model_names
LRT_to_full <- mLL_model_process_f(mLL_model_res_list,base="full")
LRT_to_NULL <- mLL_model_process_f(mLL_model_res_list,base="NULL")
LRT_res <-cbind(LRT_to_full,LRT_to_NULL[match(LRT_to_full[,1],LRT_to_NULL[,1]),c(17:18)])

mLL_model_f <- function(model_name,running_dir,source_script_dir,mc.cores=12,maxit=50,model_script_map,phen_corr,...){
  
  nXY = (nX/380e3)*(nY/380e3)*380e3 #Estimated fraction of sample overlap
  rho = phen_corr*nXY/(as.numeric(nX)*as.numeric(nY)) #Calculated phenotypic correlation due to sample overlap
  
  parmes_to_remove_map <- list(
    "full" = "",
    "only_latent" = c("a","b"),
    "only_causal" = c("pU","tX","tY"),
    "only_a" = c("pU","tX","tY","b"),
    "only_b" = c("pU","tX","tY","a"),
    "no_a" = c("a"),
    "no_b" = c("b"),
    "NULL"=  c("pU","tX","tY","a","b")
  )
  
  theta.df = matrix( data = NA, nrow=120, ncol=12) #nrows indicates how many different set of starting points to run in rslurm
  theta.df1 = t(apply(theta.df,1,FUN = function(x){c(1,1,0.5*runif(3,0,1),0.5*runif(3,0,1),runif(3,0,1)-0.5,rho)} )) #Random starting points for the 12 parameters
  par.df = data.frame(par=I(apply(theta.df1,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm
  #Set or arguments to be used in optim function, upper and lower bounds, parscale factor. 
  args.df = rbind(lower = c(0,0,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,-1,-1,-1,(-1/max(nX,nY))*10),
                  upper = c(2,2,1,1,1,1,1,1,1,1,1,(1/max(nX,nY))*10),
                  parscale=c(1,1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-7))
  
  param <- c("iX", "iY", "pX","pU","pY","h2X","h2Y","tX","tY","a","b","rho")
  
  if(model_name!="full"){
    args.df1 = args.df[,!(param %in% parmes_to_remove_map[[model_name]])]
    par.df = data.frame(par=I(apply(theta.df1[,!(param %in% parmes_to_remove_map[[model_name]])],1,as.list)))
  }
  
  script_dir_full <- paste0(source_script_dir,"/")
  script <- model_script_map[model_name]
  script_path_full <- paste0(script_dir_full,"/",script)
  output_model <- paste0(running_dir,"/",model_name,"_parscale.csv")
  source(script_path_full,local=environment())
  dir.create(running_dir,showWarnings = F,recursive = F)
  sres_null = mclapply(1:NROW(par.df),function(i,run_optim,...){print(i); run_optim(par.df[i,])},mc.cores=mc.cores,run_optim=run_optim,nX,nY,bX,bY,ld,wd,M,rho,nXY,maxit)
  rm(run_optim)
  res_temp_null = data.frame(Reduce("rbind",sres_null),stringsAsFactors = F)
  write.csv(res_temp_null, output_model, row.names = FALSE) ## Will be used to read new sp_mat (starting points in nested models)
  res_temp_null
}

mLL_model_process_f <- function(mLL_model_res_list,base="NULL"){
  param_full = rep(0,13) ##NA or 0
  names(param_full) = c("mLL","iX", "iY", "pX","pU","pY","h2X","h2Y","tX","tY","a","b","rho")
  
  sp_mat_names <- names(mLL_model_res_list)
  model_base <- mLL_model_res_list[[base]]
  sp_mat <- model_base
  parent_min = sp_mat[which(sp_mat$mLL==min(sp_mat$mLL)),]
  parent_min1 = param_full  ## bigger (always 12 params) is second element in match
  parent_min1[match(names(parent_min),names(parent_min1))] = parent_min
  AIC=2*(length(parent_min)-1)-2*parent_min$mLL
  parent_min2=cbind(model = base, t(parent_min1),AIC,nPar=length(parent_min)-1, pval_m_vs_base = NA,pval_base_vs_m = NA)
  LRT_list <- lapply(setdiff(sp_mat_names,base),function(nm){
    res_temp_null <- data.frame(mLL_model_res_list[[nm]],stringsAsFactors = F)
    res_min = res_temp_null[which(res_temp_null$mLL == min(res_temp_null$mLL)), ]
    res_min1 = param_full
    res_min1[match(names(res_min),names(res_min1))] = res_min
    df1 = length(parent_min) - length(res_min) #get degrees of freedom
    if(df1<0){
      x=res_min$mLL; y = parent_min$mLL
      LRT_m_vs_base = lr.test(x = x, y = y, df = -df1)
      LRT_base_vs_m = lr.test(x = y, y = x, df = 1)
    }else{
      df_i = 1
      x=res_min$mLL; y = parent_min$mLL
      LRT_m_vs_base = lr.test(x = x, y = y, df = 1)
      LRT_base_vs_m = lr.test(x = y, y = x, df = max(df1,1))
    }
    LRT = lr.test(x = x, y = y, df = abs(df1))
    res_pval_LRT = LRT$p.value
    AIC=2*(length(res_min)-1)-2*res_min$mLL
    LRT1_null = cbind(model = nm, t(res_min1),AIC,nPar=length(res_min)-1, pval_m_vs_base = LRT_m_vs_base$p.value,pval_base_vs_m=LRT_base_vs_m$p.value)
    LRT1_null
  })
  lrt_test <- rbind(parent_min2,Reduce("rbind",LRT_list))
  colnames(lrt_test)[grep("_base",colnames(lrt_test))]<- gsub("base",base,colnames(lrt_test)[grep("base",colnames(lrt_test))])
  data.frame(model=unlist(lrt_test[,1]),apply(LRT_to_full[,-1],2,unlist))
}
