library(JMbayes)
library(timeROC)
library(dplyr)
library(purrr)

#Varibales and Dataset
train = 0.7
t_r=t_index_r=3
scal_r=1
data("pbc2")
pbc2$event <- as.numeric(pbc2$status != "alive")
pbc2$id    <- as.numeric(as.character(pbc2$id))

#Remove missing values for dependent and independent covariates#
pbc2 <- pbc2 %>% filter(!is.na(serBilir), !is.na(spiders), !is.na(age), !is.na(drug) )

pbc2 <- pbc2 %>% group_by(id) %>% 
                 mutate(Flag=if_else((years - year) < t_index_r & event==0 , NA_real_, abs(years - year))) %>%  
                 filter(!is.na(Flag)) %>% ungroup() 

set.seed(123)
train_ids      <- sample(unique(pbc2$id), size = ceiling(n_distinct(pbc2$id)*train), replace = FALSE)
pbc2_longtrain <- filter(pbc2, id %in% train_ids)
pbc2_ep1train  <- pbc2_longtrain %>% filter(year==0) 
pbc2_longval   <- filter(pbc2, !(id %in% train_ids))

###############
#A. Training (deriving the risk tools)#
##############
#A1. based on earliest observations
CoxFit_earl <- coxph(Surv(years, event) ~ drug + age + log(serBilir) + spiders, data = pbc2_ep1train, model = TRUE)

#A2. Joint Model
MixedModelFit <- mvglmer(list(log(serBilir) ~ year + (year | id),
                              spiders ~ year + (1 | id)), 
                         data = pbc2_longtrain,
                         families = list(gaussian, binomial))
CoxFit        <- coxph(Surv(years, event) ~ drug + age, data = pbc2_ep1train, model = TRUE)
JMFit         <- mvJointModelBayes(MixedModelFit, CoxFit, timeVar = "year")

#################
#B. Validation ##
#################
#Validation Data set to be used#
df_val <- pbc2_longval %>% group_by(id) %>% 
                           mutate(age=age+year, years_org=years, years=t_r/scal_r, Y.s=years_org - year) %>% 
                           slice(n()) %>% ungroup()

#B1. Prediction using current values (of validation data) but using the equations derived from model "Earliest"#
df_val$S_earl <- exp(-predict(CoxFit_earl, newdata = df_val, type="expected")) 
earl          <- timeROC(T=df_val$Y.s, delta=df_val$event, marker=I(1 - df_val$S_earl), cause=1, times=I(t_r/scal_r))

#B2. Joint Model
#B2.1. Using age at baseline
df_age_earl <- pbc2_longval %>% group_by(id) %>% 
                                mutate(s = last(year), t_plus_s = s + (t_r/scal_r), Y.s = years-s) %>%
                                ungroup()

d.sJointE   <- df_age_earl %>% group_split(id) %>% map_dfr(function(sgr) {
  #browser()
  fit <- JMbayes::survfitJM(object = JMFit, newdata = sgr, last.time = unique(sgr$s) , survTimes = unique(sgr$t_plus_s), idVar = "id")
  res <- data.frame(id=sgr$id[1], years=sgr$years[1], s=sgr$s[1], t_plus_s=sgr$t_plus_s[1], event=sgr$event[1], fit[['summaries']][[1]])
  res
})

joint <- timeROC(T=d.sJointE$Y.s, delta=d.sJointE$event, marker=I(1 - d.sJointE$Mean), cause=1, times=I(t_r/scal_r))

#################
#C. Results    ##
#################
earl
joint