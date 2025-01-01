#
#
# R code used in "Penicillin Allergy as an Instrumental Variable for Estimating Antibiotic Effects on Resistance"
# by Yaki Saciuk (2024)
#
#
#
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(rms)
library(tidyr)
library(parallel)
library(doParallel)
library(MatchIt)
library(ivtools)
# 
#
# pairwiase comparson
results.table1     <- function(temp, j, matched=F){
        #
        # sample with replacement with 
        i.cluster.resample     <-   function(dt_){
                #
                #
                temp <-  inner_join(dt_["ID"],
                                    dplyr::select(summarise(group_by(dt_, ID), IPCW=mean(IPCW)), ID, IPCW),
                                    by="ID")
                #
                temp <-  inner_join(dplyr::select(dplyr::sample_n(temp
                                                                  , size=nrow(temp)
                                                                  , replace=T
                                                                  , weight=IPCW
                                    ), ID),
                                    dt_, 
                                    by="ID",
                                    relationship = "many-to-many")
                #
                #
                return(as.data.frame(temp))}
	      #
        worek    <-  list()
        n        <-  3000
        #
        groups   <-  c("PENICILLIN", "CEPHALOSPORINES", "MACROLIDES", "QUINOLONES", "TETRACYCLINES","TRIMETHO-SULPHA", "MONUROL", "MACRODANTIN")
        #
        comb     <-  c(c(NA), groups[groups!=j])
        # 
        #
        #
        for(ant in comb) { 
            #
            #
            # filter subset/all set 
            if(is.na(ant)){
                 temp_   <- temp
             } else {temp_    <- temp %>% filter(., TREATMENT_GROUP_DESC %in% c(j, ant))}
            #
            temp_   <- filter(temp_, DAYS_SINCE_ANTIBIOTICS<=180)
            # 
            #
            # 
            BS_IV      <-   function(x){
                            # 
                            # cluster sampling with replacement
                            sample_temp   <-  i.cluster.resample(temp_)
                            #
                            # full matching for matched analysis (Table S3)
                            if(matched){
                                m.out2     <-  matchit(Z2 ~     rcs(AGE_AT_EPISODE, 5)     +
                                                                rcs(YEAR_INDEX, 5)         +         
                                                                factor(IS_MALE)            +                       
                                                                factor(SES_GROUP)          +
                                                                factor(SOCIAL_SECTOR)      +
                                                                rcs(PREV_CULTURE_TESTS, 5) +
                                                                rcs(PREV_ISSUES, 5)        +
                                                                factor(REG_diabetes)       +
                                                                factor(REG_hypertension)   +
                                                                factor(REG_cardio_gen)     +
                                                                factor(CKD)                +
                                                                factor(REG_immunosup)
                                                          , data = sample_temp
                                                          , method = "quick" #"full"
                                                          , estimand = "ATT"
                                                          , mahvars = ~ AGE_AT_EPISODE + YEAR_INDEX + PREV_CULTURE_TESTS + PREV_ISSUES)
                                #
                                temp.model.data     <- match.data(m.out2)
                                rm(m.out2)
                                gc()}
                            #
                            # first stage
                            fitX.LZ       <-  glm(EXPOSED ~      Z2                           + 
                                                                 rcs(AGE_AT_EPISODE, 5)       +
                                                                 rcs(YEAR_INDEX, 5)           +
                                                                 factor(IS_MALE)              +                       
                                                                 factor(SES_GROUP)            +
                                                                 factor(SOCIAL_SECTOR)        +
                                                                 rcs(PREV_CULTURE_TESTS, 5)   +
                                                                 rcs(PREV_ISSUES, 5)          +
                                                                 factor(REG_diabetes)         +
                                                                 factor(REG_hypertension)     +
                                                                 factor(REG_cardio_gen)       +
                                                                 factor(CKD)                  +
                                                                 factor(REG_immunosup)																 
                                                  ,data=sample_temp
                                                  ,family=binomial(link = "probit") #   either without if linear-linear 2sls is employed
                                                  ,weights=weights
                                                  )
                            #
                            # linear outcome-exposure model
                            fitY.LX   <-  glm(R ~  EXPOSED                      +
                                                   rcs(AGE_AT_EPISODE, 5)       +
                                                   rcs(YEAR_INDEX, 5)           +
                                                   factor(IS_MALE)              +                       
                                                   factor(SES_GROUP)            +
                                                   factor(SOCIAL_SECTOR)        +
                                                   rcs(PREV_CULTURE_TESTS, 5)   +
                                                   rcs(PREV_ISSUES, 5)          +
                                                   factor(REG_diabetes)         +
                                                   factor(REG_hypertension)     +
                                                   factor(REG_cardio_gen)       +
                                                   factor(CKD)                  +
                                                   factor(REG_immunosup)         
                                           , data=sample_temp                                                             
                                           , weights=weights
                                           )
                            # 2sls estimation
                            fitIV_ts <- ivtools::ivglm(estmethod="ts"
                                                       , fitX.LZ=fitX.LZ
                                                       , fitY.LX=fitY.LX
                                                       , data=as.data.frame(model.data)
                            )
                            summary(fitIV_ts)
                            return(as.numeric(coef(fitIV)[2]))}
           #
           n.cores <- 8
           # create the cluster
           my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
           #
           ## socket cluster with 7 nodes on host 'localhost'
           #  register it to be used by %dopar%
           doParallel::registerDoParallel(cl = my.cluster)
           #
           # check if it is registered (optional)
           foreach::getDoParRegistered()
           ## [1] TRUE
           # how many workers are available? (optional)
           foreach::getDoParWorkers()
           #
           # -----------------------------------------------------
           # run in parallel mode
           BS_IVs <- foreach(i = 1:(n+1)
                             ,.combine = "c"
                             ,.packages=c("dplyr", "ivtools", "rms")
           ) %dopar% {BS_IV(i)}
           #
           # close the cluster
           parallel::stopCluster(cl = my.cluster)
           #
           # group name, either pairwise (penicillin vs some single antibiotic group) or penicillin vs all other abtibiotic groups
           G   <- if_else(is.na(ant), "All Groups", ant)
           #
           #  return the mean, 2.5% percentile and 97.5% percentile of the bootstrap sampling
           worek[[G]]     <- data.frame(N=nrow(temp_),
                                        RDiv=round(100*mean(BS_IVs),2),
                                        LLiv=round(100*quantile(BS_IVs, 0.025),2),
                                        ULiv=round(100*quantile(BS_IVs, 0.975), 2))}
           #
           return(do.call(rbind, worek))}



# AMR RD over time points
results.table2     <- function(temp, comb=c(90, 180, 270, 365)){
        #
        i.cluster.resample     <-   function(dt_){
                #
                #
                temp <-  inner_join(dt_["ID"],
                                    dplyr::select(summarise(group_by(dt_, ID), IPCW=mean(IPCW)), ID, IPCW),
                                    by="ID")
                #
                temp <-  inner_join(dplyr::select(dplyr::sample_n(temp
                                                                  , size=nrow(temp)
                                                                  , replace=T
                                                                  , weight=IPCW
                                                                  )
                                                  , ID),
                                    dt_, 
                                    by="ID",
                                    relationship = "many-to-many"
                                    )
                #
                return(as.data.frame(temp))}
        #
        worek    <-  list()
        n        <-  3000
        #
        #                     
        #
        for(ant in comb) { 
            #
            temp_   <- filter(temp, DAYS_SINCE_ANTIBIOTICS<=ant) 
            #
            # 
            #
            #
            #
            BS_IV      <-   function(x){
                            #
                            sample_temp   <-  i.cluster.resample(temp_)
                            #
                            fitX.LZ       <-  glm(EXPOSED ~      Z2                           +         
                                                                 rcs(AGE_AT_EPISODE, 5)       +
                                                                 rcs(YEAR_INDEX, 5)           +         
                                                                 factor(IS_MALE)              +                       
                                                                 factor(SES_GROUP)            +
                                                                 factor(SOCIAL_SECTOR)        +
                                                                 rcs(PREV_CULTURE_TESTS, 5)   +
                                                                 rcs(PREV_ISSUES, 5)          +
                                                                 factor(REG_diabetes)         +
                                                                 factor(REG_hypertension)     +
                                                                 factor(REG_cardio_gen)       +
                                                                 factor(CKD)                  +
                                                                 factor(REG_immunosup)																	 
                                                  ,data=sample_temp
					                                        ,family=binomial(link = "probit")
					                                        ,weights=weights
						                )
                            #
                            #
                            #
                            fiY.LX   <-  lm(R ~    EXPOSED                      +
                                                   rcs(AGE_AT_EPISODE, 5)       +
                                                   rcs(YEAR_INDEX, 5)           +         
                                                   factor(IS_MALE)              +                       
                                                   factor(SES_GROUP)            +
                                                   factor(SOCIAL_SECTOR)        +
                                                   rcs(PREV_CULTURE_TESTS, 5)   +
                                                   rcs(PREV_ISSUES, 5)          +
                                                   factor(REG_diabetes)         +
                                                   factor(REG_hypertension)     +
                                                   factor(REG_cardio_gen)       +
                                                   factor(CKD)                  +
                                                   factor(REG_immunosup)         
                                           , data=sample_temp                                                                
                                           , weights = weights
                                           )
                            #
                            return(as.numeric(coef(fitIV)[2]))}
           #
           n.cores <- 10
           # create the cluster
           my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
           #
           ## socket cluster with 7 nodes on host 'localhost'
           #  register it to be used by %dopar%
           doParallel::registerDoParallel(cl = my.cluster)
           #
           # check if it is registered (optional)
           foreach::getDoParRegistered()
           ## [1] TRUE
           # how many workers are available? (optional)
           foreach::getDoParWorkers()
           #
           # -----------------------------------------------------
           BS_IVs <- foreach(i = 1:(n+1)
                             ,.combine = "c"
                             ,.packages=c("dplyr", "ivtools", "rms")
           ) %dopar% {BS_IV(i)}
           #
           # 
           parallel::stopCluster(cl = my.cluster)
           #
           #
           G   <- paste("Up to:", ant, " days", sep="")
           #
           worek[[G]]     <- data.frame(N=nrow(temp_),
                                        RDiv=round(100*mean(BS_IVs),2),
                                        LLiv=round(100*quantile(BS_IVs, 0.025),2),
                                        ULiv=round(100*quantile(BS_IVs, 0.975), 2))}
           #
           return(do.call(rbind, worek))}








