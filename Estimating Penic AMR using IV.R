#
# R code used in "Penicillin Allergy as an Instrumental Variable for Estimating Antibiotic Effects on Resistance"
# by Yaki Saciuk (2024)
#
rm(list=ls())
gc()
#
library(data.table)
library(dplyr)
library(kableExtra)
library(table1)
library(ggplot2)
library(lubridate)
library(stringr)
library(rms)
library(extrafont)
loadfonts(device = "win")


library(Hmisc)
library(tidyr)
library(knitr)
library(ivtools)
#
#
#
#
# Data extraction
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
#
path_2         <-  "H://phd//Estimating AMR using IV//data"
#
#
replace_na_i   <-  function(x) fifelse(is.na(x), 0, 1)
#
#
#
vitek.results    <-  fread("//yourpath//vitek.results.csv")                                                      %>% 
                     as.data.frame(.)                                                                            %>%                
                     mutate(., Z2=if_else(Penicillin==1, 0, 1))     # Mark IV variable Z2 (1 for allergic, 0 for non-allergic to penicillin)
				   
antibiotics.use    <-  fread("//yourpath//antibiotics.use.csv") %>% as.data.frame(.)
#
#
# 
#
depo   <- list()  # pairwise (exposure: penicillin vs a single other group) analysis (Table S5)
depo.p <- list()  # IPCW probabilities distribution (Table S6)
depo.2 <- list()  # over time (Table 2)
depo.3 <- list()  # penicillins vs quionolones (Table 3)
n       <- 0
#
# Calculate the causal effect of previous penicillin use on the future AMR of five testes antibiotics
set.seed(5555)
#
for(j in c("Amoxicillin_CA", "Ampicillin", "Piperacillin_Tazobactam", "Gentamicin", "Ciprofloxacin")){
    #
    n             <-  n + 1
    j_tests       <-  filter(vitek.results, OMA147==j)  # filter the relevant amr results of the specific antibiotic (j)
	#
	# calculate inverse probability of censoring for this particular antibiotic (j)
    antibiotics.use_a  <-  mutate(antibiotics.use, CENSORED=if_else(CENSORED==1, 1, if_else(REQUEST_NUMBER %in% j_tests$REQUEST_NUMBER, 0, 1))) 
    #
	# choose only the first amr test after the antibiotic use
    antibiotics.use_a  <-  as.data.frame(mutate(group_by(antibiotics.use_a, ID, CODE_ID, ISSUE_DATE), VITEK_NUM=row_number()))
    antibiotics.use_a  <-  filter(antibiotics.use_a, VITEK_NUM==1 | CENSORED==0)
    #
    # calculate IPCW
    iwpc.logit <-    glm(CENSORED           ~ rcs(ISSUE_YEAR, 5)                                  +
	                                            factor(IS_MALE)                                     +
                                              factor(TREATMENT_GROUP_DESC)                        + 
                                              factor(SES_GROUP)                                   + 
                                              rcs(AGE_AT_EPISODE, 5)                              +
                                              factor(REG_cardio_gen)                              +
                                              factor(CKD)                                         +
                                              factor(REG_diabetes)                                +
                                              factor(REG_hypertension)                            +
                                              factor(REG_immunosup)                               +
                                              factor(REG_obesity)                                 +
                                              factor(SOCIAL_SECTOR)
                         , data=antibiotics.use_a 
                         , family="binomial")
    #
	  # calculate the marginal prob. for stabilizing the IPCW
    marginal_treat      <- 1- mean(antibiotics.use_a$CENSORED)
    antibiotics.use_a   <-    mutate(antibiotics.use_a, P=predict(iwpc.logit, newdata=antibiotics.use_a, type="response"))           %>%
                              mutate(., IPCW=(1-P)^-1, SIPCW=marginal_treat*((1-P)^-1))
    
    # keep the results of IPCW probabilities distribution of j (table s6)
    depo.p[[n]]  <-   do.call(rbind, tapply(100*(1-antibiotics.use_a$P), antibiotics.use_a$CENSORED, summary))                       %>%
                      as.data.frame(.)                                                                                               %>%
                      mutate(., AMR_anti=j)
    
    # 
    temp.model.data  <-  filter(vitek.results, OMA147==j)                                                                            %>%
	                       mutate(., P=predict(iwpc.logit, newdata=., type="response"))                                                %>%
                         mutate(., IPCW=(1-P)^-1, SIPCW=marginal_treat*((1-P)^-1))                                                   %>%
                         mutate(., EXPOSED=if_else(TREATMENT_GROUP_DESC=="PENICILLIN", 1, 0))

    # Calculate RD with 95% bootstrap C.I. using IV analysis. result.table(i) can be found in aux.R in this depository
	  # results.table(i) function can be found below
    # (i) pairwise (table s5)
    res3       <-  results.table1(temp.model.data, "PENICILLIN") 
    depo[[n]]  <-  mutate(res3, amr_anti=j)
    #
    # (ii) overtime (RD in 90, 180, 270, 360 days after the antibiotic use) (table 2)
    res4  <-  results.table2(temp.model.data)
    depo.2[[n]]  <- mutate(res4, amr_anti=j)
    #
    # penic + quino. (table 3)
    temp.model.data   <- filter(temp.model.data, TREATMENT_GROUP_DESC %in% c("PENICILLIN", "QUINOLONES"))
    res5  <-  results.table2(temp.model.data)
    depo.3[[n]]  <- mutate(res5, amr_anti=j)
    #
    #
}
#
#
#
#
# ===============================================================================================================================================================
# read penicillin sub-classification (amozicillin, amoxicillin ca and other penicillin) of exposure for sub-analysis
antibiotics          <-  fread("//yourpath//antibiotics.use.csv") %>% as.data.frame(.)
#
#
vitek.results_b      <-  left_join(mutate(vitek.results_b_a, ISSUE_DATE=as.Date(ISSUE_DATE))
                                  ,mutate(antibiotics, ISSUE_DATE=as.Date(ISSUE_DATE))
                                  ,by=c("ID", "CODE_ID", "ISSUE_DATE"))
#
worek  <- list()
n      <- 0
#
for(i in c(1:4)){
    for(j in c("Amoxicillin_CA", "Ampicillin", "Piperacillin_Tazobactam", "Gentamicin")){
    #
    n             <-  n + 1
    j_tests       <-  filter(vitek.results_b, OMA147==j)
    antibiotics.use_ba  <-  mutate(antibiotics.use_b, CENSORED=if_else(CENSORED==1, 1,
                                                            if_else(REQUEST_NUMBER %in% j_tests$REQUEST_NUMBER, 0, 1)))
    #
    antibiotics.use_ba  <-  as.data.frame(mutate(group_by(antibiotics.use_ba, ID, CODE_ID, ISSUE_DATE), VITEK_NUM=row_number()))
    antibiotics.use_ba  <-  filter(antibiotics.use_ba, VITEK_NUM==1 | CENSORED==0)
    #
    #
    iwpc.logit <-    glm(CENSORED           ~ rcs(ISSUE_YEAR, 3)                                  +
                                              factor(TREATMENT_GROUP_DESC)                        + 
                                              factor(IS_MALE)                                     +
                                              factor(SES_GROUP)                                   + 
                                              rcs(AGE_AT_EPISODE, 3)                              +
                                              factor(REG_cardio_gen)                              +
                                              factor(CKD)                                         +
                                              factor(REG_diabetes)                                +
                                              factor(REG_hypertension)                            +
                                              factor(REG_immunosup)                               +
                                              factor(REG_obesity)                                 +
                                              factor(SOCIAL_SECTOR)
                         ,data=antibiotics.use_ba 
                         ,family="binomial")
    #
    # calculate the particular P and IPCW
    marginal_treat    <- 1- mean(antibiotics.use_ba$CENSORED)
    #
    #
    vitek.results_ba <-   mutate(vitek.results_ba, P=predict(iwpc.logit, newdata=vitek.results_ba, type="response"))                        %>%
                          mutate(., IPCW=(1-P)^-1, SIPCW=marginal_treat*((1-P)^-1))
    
    #
    temp.model.data   <-  vitek.results_ba
    temp.model.data   <-  mutate(temp.model.data, EXPOSED=if_else(TREATMENT_GROUP_DESC=="PENICILLIN", 1, 0))
    temp.model.data   <-  filter(temp.model.data, OMA147==j)
    #
    if(i==1){
         temp.model.data   <- filter(temp.model.data, is.na(AMOX_CA) | AMOX_CA==1)
         x <-  "amox ca"                            
    } else if(i==2){
         temp.model.data   <- filter(temp.model.data, is.na(AMOX) | AMOX==1)
         x <-  "amox"
    } else{
         temp.model.data   <- filter(temp.model.data, is.na(N) | AMOX==1 | PENIC==1)
         x <-  "amox_penic"
    }
    #
    
    # (ii) overtime
    res4  <-  results.table2(temp.model.data)
    # 
    worek[[n]]  <- mutate(res4, amr_anti=j, exp_penic=x, exp_penic_cd=i)
    #
    #
    #
    }
}





