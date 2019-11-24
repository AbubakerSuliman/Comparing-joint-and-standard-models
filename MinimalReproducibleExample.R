library(JMbayes)
library(timeROC)
library(dplyr)
library(purrr)
#Varibales and Dataset
train = 0.7
t_r=t_index_r=3
scal_r=1

pbc2$event <- as.numeric(pbc2$status != "alive")
pbc2$id <- as.numeric(as.character(pbc2$id))
#Remove missing values for dependent and independent covariates#
pbc2 = pbc2 %>% filter(!is.na(serBilir), !is.na(spiders), !is.na(age), !is.na(drug) )

pbc2 <- pbc2 %>% group_by(id) %>% 
  mutate(Flag=if_else((years - year) < t_index_r & event==0 ,NA_real_,abs(years - year))) %>%  
  filter(!is.na(Flag)) %>% 
  ungroup() 

set.seed(123)
train_ids <- sample(unique(pbc2$id), size = ceiling(n_distinct(pbc2$id)*train), replace = FALSE)
pbc2_longtrain <- filter(pbc2, id %in% train_ids)
pbc2_ep1train <- pbc2_longtrain %>% filter(year==0) %>% ungroup()
pbc2_longval <- filter(pbc2, !(id %in% train_ids))


###############
#A. Training (deriving the risk tools)#
##############
#A1. based on earliest observations
CoxFit_earl <- coxph(Surv(years, event) ~ drug + age + log(serBilir) + spiders, data = pbc2_ep1train, model = TRUE)

#A2. based on current observations
pbc2_longtrain_curr <- pbc2_longtrain %>% group_by(id) %>% 
  mutate(age=age+year, years = years - year) %>% 
  slice(n())

CoxFit_curr <- coxph(Surv(years, event) ~ drug + age + log(serBilir) + spiders, data = pbc2_longtrain_curr, model = TRUE)

#A3. Joint Model
MixedModelFit <- mvglmer(list(log(serBilir) ~ year + (year | id),
                              spiders ~ year + (1 | id)), data = pbc2_longtrain,
                         families = list(gaussian, binomial))

CoxFit <- coxph(Surv(years, event) ~ drug + age, data = pbc2_ep1train, model = TRUE)
JMFit <- mvJointModelBayes(MixedModelFit, CoxFit, timeVar = "year")



#################
#B. Validation ##
#################

#Validation Data set to be used#
pbc2_longval.temp = pbc2_longval %>% group_by(id) %>% 
  mutate(age=age+year, years_org=years, years=t_r/scal_r) %>% slice(n())

#B1. Prediction using current values (of validation data) but using the equations derived from model "Earliest"#
pbc2_longval.temp$S_t = exp(-predict(CoxFit_earl,newdata = pbc2_longval.temp , type="expected")) 
index = which(names(pbc2_longval.temp) %in% c("id", "years_org", "year", "event", "S_t"))
d.s = pbc2_longval.temp[, index]
d.s$marker = d.s$S_t
d.s$Y.s <- d.s$years_org-d.s$year
a.earl =  timeROC(T=d.s$Y.s, delta=d.s$event, marker=I(1 - d.s$marker), cause=1, times=I(t_r/scal_r))

#B2.On current values but using the equations derived from model "Current"

pbc2_longval.temp$S_t = exp(-predict(CoxFit_curr,newdata = pbc2_longval.temp , type="expected"))  #see 
index = which(names(pbc2_longval.temp) %in% c("id", "years_org", "year", "event", "S_t"))
d.s = pbc2_longval.temp[, index]
d.s$marker = d.s$S_t
d.s$Y.s <- d.s$years_org-d.s$year
a.curr =  timeROC(T=d.s$Y.s, delta=d.s$event, marker=I(1 - d.s$marker), cause=1, times=I(t_r/scal_r))


#B3. Joint Model
df.temp1 <- pbc2_longval %>% group_by(id) %>% arrange(desc(year), .by_group = TRUE) %>% slice(1) %>% ungroup() %>% mutate(s = year) %>% mutate(t_plus_s= s + (t_r/scal_r)) %>% dplyr::select(id, s, t_plus_s)

# #B3.1 New
# d.sJoint <- df.temp2 %>% group_split(id) %>% map_dfr(function(sg) { 
#   #browser()
#   fit <- JMbayes::survfitJM(object = JMFit, newdata = sg, last.time = unique(sg$s) , survTimes = unique(sg$t_plus_s), idVar = "id")
#   res <- data.frame(id=sg$id[1], years=sg$years[1], s=sg$s[1], t_plus_s=sg$t_plus_s[1], event=sg$event[1], fit[['summaries']][[1]])
#   res
#   })
# 
# d.sJoint$Y.s <- d.sJoint$years-d.sJoint$s
# indexJ = which(names(d.sJoint) %in% c("Mean"))
# d.sJoint$marker = d.sJoint[ ,indexJ, drop=TRUE]
# a.Joint = timeROC(T=d.sJoint$Y.s, delta=d.sJoint$event, marker=I(1 - d.sJoint$marker), cause=1, times=I(t_r/scal_r))

#B3.2 Using age at baseline
df.temp0_age_earl = pbc2_longval %>% arrange(id)
df.temp2_age_earl = merge(df.temp0_age_earl, df.temp1, by.x = "id", by.y = "id", all = TRUE)
sub = unique(df.temp2_age_earl$id)
df_result1 = list()
for(i in 1:length(sub))
{
  a = filter(df.temp2_age_earl, id %in% sub[i])
  sprobs <- JMbayes::survfitJM(object = JMFit, newdata = a, last.time = unique(a$s) , survTimes = unique(a$t_plus_s), idVar = "id")
  stopifnot("summaries" %in% names(sprobs))
  df_result1[[i]]  <- bind_rows(sprobs[['summaries']]) %>% t() %>% as.data.frame() %>% rownames_to_column('id') %>% dplyr::select(id, V2, V3) %>% rename(marker_mean=V2, marker_median=V3)
  rm(a); rm(sprobs)
}
df_result =  bind_rows(df_result1)
sprob_st = df_result
others = df.temp1
D_ep1_s = pbc2_longval[pbc2_longval$year==0,]

D_new = merge(D_ep1_s, sprob_st, by.x = "id", by.y = "id", all = TRUE)
D_new = merge(D_new, others,by.x = "id", by.y = "id", all = TRUE )
index = which(names(D_new) %in% c("id", "years", "event", "marker_mean", "marker_median","s","t_plus_s"))
d.s = D_new[, index]
d.s$Y.s <- d.s$years-d.s$s
d.s = dplyr::filter(d.s, Y.s > 0)

index = which(names(d.s) %in% c("marker_mean"))
d.s$marker = d.s[,index, drop=TRUE]
a.joint.earl =  timeROC(T=d.s$Y.s,delta=d.s$event,marker=I(1 - d.s$marker),cause=1,times=I(t_r/scal_r))


#B3.3 Using age at last visit
#93.39
#%>% group_by(id) %>% mutate(spiders=last(spiders), serBilir=last(serBilir)) %>% ungroup()
df.temp0_age_curr = pbc2_longval %>% arrange(id) %>% mutate(age=age+year) 
df.temp2_age_curr = merge(df.temp0_age_curr, df.temp1, by.x = "id", by.y = "id", all = TRUE)

sub = unique(df.temp2_age_curr$id)
df_result1 = list()
for(i in 1:length(sub))
{
  a = filter(df.temp2_age_curr, id %in% sub[i])
  sprobs <- JMbayes::survfitJM(object = JMFit, newdata = a, last.time = unique(a$s) , survTimes = unique(a$t_plus_s), idVar = "id")
  stopifnot("summaries" %in% names(sprobs))
  df_result1[[i]]  <- bind_rows(sprobs[['summaries']]) %>% t() %>% as.data.frame() %>% rownames_to_column('id') %>% dplyr::select(id, V2, V3) %>% rename(marker_mean=V2, marker_median=V3)
  rm(a); rm(sprobs)
}
df_result =  bind_rows(df_result1)
sprob_st = df_result
others = df.temp1
D_ep1_s = pbc2_longval[pbc2_longval$year==0,]

D_new = merge(D_ep1_s, sprob_st, by.x = "id", by.y = "id", all = TRUE)
D_new = merge(D_new, others,by.x = "id", by.y = "id", all = TRUE )
index = which(names(D_new) %in% c("id", "years", "event", "marker_mean", "marker_median","s","t_plus_s"))
d.s = D_new[, index]
d.s$Y.s <- d.s$years-d.s$s
d.s = dplyr::filter(d.s, Y.s > 0)

index = which(names(d.s) %in% c("marker_mean"))
d.s$marker = d.s[,index, drop=TRUE]
a.joint.curr =  timeROC(T=d.s$Y.s,delta=d.s$event,marker=I(1 - d.s$marker),cause=1,times=I(t_r/scal_r))

a.earl
a.curr
a.joint.earl
a.joint.curr