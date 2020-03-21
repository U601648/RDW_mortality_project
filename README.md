
##########################################################################################
##                                                                                      ##
##                        Appendix 2 - Statistical analysis R code                      ##
##                                                                                      ##
##########################################################################################

## RDW_mortality_project
## Describing the added predicted value for RDW in emergency laparotomy patients
## Load R packages ##
library(car)
library(DescTools)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(readr)
library(rms)
library(survival)
library(survminer)

## Clear global environment ##
rm(list = ls())

## Import data set ##
rdw.vars<-read_csv("~/Desktop/RDW_project/RDW_csv_data/rdw_project_data.csv")

## Change names of data column ##
rdw.vars<-rdw.vars %>% 
  rename(time30=time_1,
         time1year=time_2,
         time2year=time3)

## Change categorical variables to factors ##
rdw.vars<-mutate_at(rdw.vars,vars(sex,nela_yr, asa,rdw_quartiles,gcs,ecg,cardiac_signs,
                                  resp_signs,op_sev,procedures_no, pred_tbl,pred_perit_soil,
                                  malignancy, ncepod_cat, indc_class),as.factor)


##########################################################################################
## (1) Survivors vs non-survivors by RDW quartile                                       ##
##########################################################################################

## age ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(age_at_adm, probs=0.5),
            `25%`=quantile(age_at_adm, probs=0.25),
            `75%`=quantile(age_at_adm, probs=0.75),
            avg=mean(age_at_adm),
            n=n())
kruskal.test(age_at_adm~rdw_quartiles, data=rdw.vars)

## sex ##
rdw.sex<-table(rdw.vars$sex,rdw.vars$rdw_quartiles)
addmargins(rdw.sex)
round(100*prop.table(rdw.sex,2),digits = 1)
fisher.test(rdw.sex)

## surgical indication ##
ind_count<-rdw.vars %>% 
  select(rdw_quartiles, indc_class) %>% 
  group_by(rdw_quartiles,indc_class) %>% 
  summarise(count=n())

ind.table<-matrix(ind_count$count, ncol=4, nrow = 6)
colnames(ind.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(ind.table)<-c("bleeding", "colitis", "ischaemia", "obstruction", "other", "sepsis")
fisher.test(ind.table, simulate.p.value = TRUE, B=1e7)

## nela risk score ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(nela_risk, probs=0.5),
            `25%`=quantile(nela_risk, probs=0.25),
            `75%`=quantile(nela_risk, probs=0.75),
            avg=mean(nela_risk),
            n=n())
kruskal.test(nela_risk~rdw_quartiles, data=rdw.vars)

## asa ##
rdw.asa<-table(rdw.vars$asa, rdw.vars$rdw_quartiles)
addmargins(rdw.asa)
round(100*prop.table(rdw.asa,2),digits = 1)

asa_count<-rdw.vars %>% 
  select(rdw_quartiles, asa) %>% 
  group_by(rdw_quartiles,asa) %>% 
  summarise(count=n())

asa.table<-matrix(c(62,51,33,29,31,39,53,58), ncol = 4, byrow = TRUE)
colnames(asa.table)<-c("rdw1", "rdw2","rdw3","rdw4")
rownames(asa.table)<-c("ASA<3", "ASA≥3")
asa.table
chisq.test(asa.table)

## urgency of surgery ##
rdw.ncepod<-table(rdw.vars$ncepod_cat,rdw.vars$rdw_quartiles)
addmargins(rdw.ncepod)
round(100*prop.table(rdw.ncepod,2),digits = 1)
round(100*(24/86), digits = 1)

urg_count<-rdw.vars %>% 
  select(rdw_quartiles, ncepod_cat) %>% 
  group_by(rdw_quartiles,ncepod_cat) %>% 
  summarise(count=n())

urg.table<-matrix(c(12,12,23,19,39,31,31,31,36,42,24,34,6,5,8,3),ncol = 4, byrow = TRUE)
colnames(urg.table)<-c("rdw1", "rdw2","rdw3","rdw4")
row.names(urg.table)<-c("expedited", "urgent6_18", "urgen2_6","immediate")
urg.table<-as.table(urg.table)
urg.table
chisq.test(urg.table)#chi square test

## ecg ##
rdw.ecg<-table(rdw.vars$ecg,rdw.vars$rdw_quartiles)
addmargins(rdw.ecg)
round(100*prop.table(rdw.ecg,2),digits = 1)

ecg_count<-rdw.vars %>% 
  select(rdw_quartiles, ecg) %>% 
  group_by(rdw_quartiles,ecg) %>% 
  summarise(count=n())

ecg.table<-matrix(ecg_count$count, ncol=4, nrow = 3)
colnames(ecg.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(ecg.table)<-c("none", "af1", "af2")
fisher.test(ecg.table)

## cardiac signs ##
rdw.cardiac<-table(rdw.vars$cardiac_signs,rdw.vars$rdw_quartiles)
addmargins(rdw.cardiac)
round(100*prop.table(rdw.cardiac,2),digits = 1)

cardiac_count<-rdw.vars %>% 
  select(rdw_quartiles, cardiac_signs) %>% 
  group_by(rdw_quartiles,cardiac_signs) %>% 
  summarise(count=n())

cardiac.table<-matrix(cardiac_count$count, ncol=4, nrow = 4)
colnames(cardiac.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(cardiac.table)<-c("none", "failure1", "failure2", "failure3")
chisq.test(cardiac.table)
fisher.test(cardiac.table,simulate.p.value=TRUE, B=1e7)

## respiratory signs ##
rdw.resp<-table(rdw.vars$resp_signs,rdw.vars$rdw_quartiles)
addmargins(rdw.resp)
round(100*prop.table(rdw.resp,2),digits = 1)

resp_count<-rdw.vars %>% 
  select(rdw_quartiles, resp_signs) %>% 
  group_by(rdw_quartiles, resp_signs) %>% 
  summarise(count=n())

resp.table<-matrix(resp_count$count, ncol=4, nrow = 4)
colnames(resp.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(resp.table)<-c("no dysp", "dysp_ex", "limit_dysp", "dysp_rest")
fisher.test(resp.table,simulate.p.value=TRUE, B=1e7)

## hb ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(hb, probs=0.5),
            `25%`=quantile(hb, probs=0.25),
            `75%`=quantile(hb, probs=0.75),
            avg=mean(hb),
            n=n())
kruskal.test(hb~rdw_quartiles, data=rdw.vars)

## mcv ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(mcv, probs=0.5),
            `25%`=quantile(mcv, probs=0.25),
            `75%`=quantile(mcv, probs=0.75),
            avg=mean(mcv),
            n=n())
kruskal.test(mcv~rdw_quartiles, data=rdw.vars)

## cr ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(cr, probs=0.5),
            `25%`=quantile(cr, probs=0.25),
            `75%`=quantile(cr, probs=0.75),
            avg=mean(cr),
            n=n())
kruskal.test(cr~rdw_quartiles, data=rdw.vars)

## urea ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(urea, probs=0.5),
            `25%`=quantile(urea, probs=0.25),
            `75%`=quantile(urea, probs=0.75),
            avg=mean(urea),
            n=n())
kruskal.test(urea~rdw_quartiles, data=rdw.vars)

## sodium ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(sodium, probs=0.5),
            `25%`=quantile(sodium, probs=0.25),
            `75%`=quantile(sodium, probs=0.75),
            avg=mean(sodium),
            n=n())
kruskal.test(sodium~rdw_quartiles, data=rdw.vars)

## wcc ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(wcc, probs=0.5),
            `25%`=quantile(wcc, probs=0.25),
            `75%`=quantile(wcc, probs=0.75),
            avg=mean(wcc),
            n=n())
kruskal.test(wcc~rdw_quartiles, data=rdw.vars)

## sbp ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(sbp, probs=0.5, na.rm = TRUE),
            `25%`=quantile(sbp, probs=0.25,na.rm = TRUE),
            `75%`=quantile(sbp, probs=0.75, na.rm = TRUE),
            avg=mean(sbp),
            n=n())
kruskal.test(sbp~rdw_quartiles, data=rdw.vars)

## pulse ##
rdw.vars %>% group_by(rdw_quartiles) %>%
  summarise(`50`=quantile(pulse, probs=0.5,na.rm = TRUE),
            `25%`=quantile(pulse, probs=0.25,na.rm = TRUE),
            `75%`=quantile(pulse, probs=0.75,na.rm = TRUE),
            avg=mean(pulse),
            n=n())
kruskal.test(pulse~rdw_quartiles, data=rdw.vars)

## op severity ##
rdw.sev<-table(rdw.vars$op_sev, rdw.vars$rdw_quartiles)
addmargins(rdw.sev)
round(100*prop.table(rdw.sev,2),digits = 1)

sev_count<-rdw.vars %>% 
  select(rdw_quartiles, op_sev) %>% 
  group_by(rdw_quartiles,op_sev) %>% 
  summarise(count=n())

sev.table<-matrix(cardiac_count$count, ncol=4, nrow = 2)
colnames(cardiac.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(cardiac.table)<-c("major", "major+")
fisher.test(sev.table)

## soiling ##
rdw.soiling<-table(rdw.vars$pred_perit_soil, rdw.vars$rdw_quartiles)
addmargins(rdw.soiling)
round(100*prop.table(rdw.soiling,2),digits = 1)

soiling_count<-rdw.vars %>% 
  select(rdw_quartiles, pred_perit_soil) %>% 
  group_by(rdw_quartiles,pred_perit_soil) %>% 
  summarise(count=n())

soiling.table<-matrix(cardiac_count$count, ncol=4, nrow = 4)
colnames(sev.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(sev.table)<-c("none", "serous", "local", "bowel_cont")
fisher.test(cardiac.table, workspace = 2e8)

## tbl ##
rdw.tbl<-table(rdw.vars$pred_tbl, rdw.vars$rdw_quartiles)
addmargins(rdw.tbl)
round(100*prop.table(rdw.tbl,2),digits = 1)

tbl_count<-rdw.vars %>% 
  select(rdw_quartiles, pred_tbl) %>% 
  group_by(rdw_quartiles,pred_tbl) %>% 
  summarise(count=n())

tbl.table<-matrix(cardiac_count$count, ncol=4, nrow = 4)
colnames(tbl.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(tbl.table)<-c("100", "100_500", "500_1000", "1000")
fisher.test(tbl.table, workspace = 2e8)

## malignancy ##
rdw.mal<-table(rdw.vars$malignancy, rdw.vars$rdw_quartiles)
addmargins(rdw.mal)
round(100*prop.table(rdw.mal,2),digits = 1)

mal_count<-rdw.vars %>% 
  select(rdw_quartiles, malignancy) %>% 
  group_by(rdw_quartiles,malignancy) %>% 
  summarise(count=n())

mal.table<-matrix(cardiac_count$count, ncol=4, nrow = 4)
colnames(mal.table)<-c("rdw1", "rdw2", "rdw3", "rdw4")
rownames(mal.table)<-c("none", "local", "distant", "mets")
fisher.test(mal.table, simulate.p.value = TRUE, B=1e7)

## 30-day mortality ##
rdw.30<-table(rdw.vars$mort30days, rdw.vars$rdw_quartiles)
addmargins(rdw.30)
round(100*prop.table(rdw.30,2),digits = 1)
fisher.test(rdw.30)

## overall mortality ##
rdw.tot<-table(rdw.vars$mort_stat1, rdw.vars$rdw_quartiles)
addmargins(rdw.tot)
round(100*prop.table(rdw.tot,2),digits = 1)
fisher.test(rdw.tot)

##########################################################################################
## (2) Outcome plots                                                                    ##
##########################################################################################

dd<-datadist(rdw.vars); options(datadist = "dd")

## 30-day survival plot following emergency laparotomy ##
## 30 day cumulative events plot ##
units(rdw.vars$time30)<-"Day" # important for calibration plot

S<-Surv(rdw.vars$time30, rdw.vars$mort30days)

fit<-survfit(S~1, data=rdw.vars)
ggsurv<-ggsurvplot(fit, fun = "event", conf.int = T, risk.table = T, 
                   title=   "30-day cumulative mortality rate following emergency laparotomy",
                   xlim=c(0,30), xlab="Follow-up time in days",
                   ylim=c(0,0.13),
                   break.y.by=0.01,
                   surv.plot.height=0.8,
                   tables.height=0.10,
                   cumevents = T,
                   cumevents.height=0.10,
                   font.main=c(12,"bold","black"),
                   fontsize=4,
                   palette="black", legend="none", 
                   tables.theme = theme_cleantable(),
                   tables.y.text=F)

ggsurv$table<-ggpar(ggsurv$table,
                    font.title = c(10,"bold", "black"), 
                    tickslab = F)

ggsurv$cumevents<-ggsurv$cumevents+labs(
  title="Cumulative mortality rate")

ggsurv$cumevents<-ggpar(ggsurv$cumevents,
                        font.title = c(10,"bold","black"),
                        tickslab = F)

print(ggsurv)

## Overall cumulative events plot - to account for low mortality rate ##
units(rdw.vars$time)<-"Year"
S1<-Surv(rdw.vars$time/365.25, rdw.vars$mort_stat3)

fit1<-survfit(S1~1, data=rdw.vars)

ggsurv1<-ggsurvplot(fit1, fun = "event", conf.int = T, risk.table = T,
                    title= "Overall cumulative mortality rate following emergency laparotomy",
                    xlim=c(0,1580/365.25),xlab="Follow-up time in years",
                    surv.plot.height=0.8,
                    tables.height = 0.10,
                    cumevents = T,
                    cumevents.height=0.10,
                    font.main=c(12,"bold","black"),
                    fontsize=4,
                    palette = "black", legend="none",
                    tables.theme = theme_cleantable(),
                    tables.y.text=FALSE)


ggsurv1$table<-ggpar(ggsurv1$table,
                     font.title = c(10,"bold", "black"), 
                     tickslab = F)

ggsurv1$cumevents<-ggsurv1$cumevents+labs(
  title="Cumulative mortality rate")

ggsurv1$cumevents<-ggpar(ggsurv1$cumevents,
                         font.title = c(10,"bold","black"),
                         tickslab = F)

print(ggsurv1)

##########################################################################################
## (3.1) RDW and mortality - boxplots                                                   ##
##########################################################################################

## RDW as a continous variable - comparing median values survivors vs non-survivors ##
dbd<-Dual_boxplot_data <- read_csv("~/Desktop/RDW_project/RDW_csv_data/Dual_boxplot_data.csv")

## Rename row names ##
dbd<-dplyr::mutate(dbd,Mortality1=ifelse(Mortality=="mort30days", "30 day mortality", "Overall mortality"))

dbd<-dbd %>% 
  mutate(RDW=Winsorize(RDW,probs = c(0.00,0.95)))

## Basic boxplot ##
rdw_box3<-ggplot(dbd, aes(x= Mortality1, y=RDW, fill=YesNo))+
  stat_boxplot(geom = "errorbar", width=0.7)+
  geom_boxplot(width=0.7, outlier.shape = NA)
theme_bw()


## Change the colour to grey ##
rdw_box3<-rdw_box3+scale_fill_brewer(palette = "Greys")+theme_pubr(legend = "right")

## Add the Wilcoxon test onto the graph ##
rdw_box3<-rdw_box3+stat_compare_means(method = "wilcox.test", label.y = 22.5)

## Change the legend label ##
rdw_box3<-rdw_box3+labs(fill="Survival status")

## Change x-axis and y-axis labels ##
rdw_box3<-rdw_box3+scale_x_discrete(name= "Follow-up period")+
  scale_y_continuous(name="RDW (%)")

## Add title ##
rdw_box3<-rdw_box3+ggtitle("RDW distribution by survival status")

## Centre title ##
rdw_box3<-rdw_box3+theme(plot.title = element_text(hjust = 0.5))

## Final version of publishable boxplot ##
rdw_box3

rm(dbd) ## remove dual boxplot data frame 

##########################################################################################
## (3.2) RDW and mortality - Kaplan-Meier curves                                        ##
##########################################################################################

## Kaplan-Meier survival curves by RDW quartiles ##
## Comparing RDW as a categorical variable for survivors vs non-survivors ##

km<-KM_30_day <- read_csv("~/Desktop/RDW_project/RDW_csv_data/KM_30_day.csv")

fit.km<-survfit(Surv(time_1,mort30days)~rdw_quartiles, data = km)

## 30-day event rate plot ##
ggsurvplot(fit.km, 
           pval=F,
           fun = "event",
           pval.coord=c(2,0.125),
           legend.labs=c("RDW1", "RDW2","RDW3","RDW4"),
           title="30-day event rate post emergency laparotomy by RDW quartiles",
           font.title=c(12,"bold","black"),
           xlab="Survival time in days",
           ylim=c(0, 0.15),
           xlim=c(0,30),
           linetype = c(1,2,3,4),
           palette = c("grey"),
           risk.table = TRUE,
           risk.table.height=0.25,
           risk.table.y.text.col=F,
           risk.table.y.text=T,
           data = km,
           ggtheme = theme_pubr()+theme(plot.title = element_text(hjust = 0.5)))

## Fit mortality for total study period ##

fit.km.total<-survfit(Surv(time,mort_stat3)~rdw_quartiles, data = rdw.vars)

## Cumulative event rate for study period by RDW quartiles ##

fit.km.total<-survfit(Surv(time,mort_stat3)~rdw_quartiles, data = rdw.vars)
ggsurvplot(fit.km.total, 
           pval=F,
           fun = "event",
           pval.coord=c(1.8,0.4),
           legend.labs=c("RDW1", "RDW2","RDW3","RDW4"),
           title="Overall event rate post emergency laparotomy by RDW quartiles",
           font.title=c(12,"bold","black"),
           xlab="Survival time in years",
           ylim=c(0,0.4),
           xscale=365.25,
           break.x.by=365.25,
           break.y.by=0.1,
           linetype = c("strata"),
           palette = c("grey"),
           risk.table = TRUE,
           risk.table.height=0.25,
           risk.table.y.text.col=F,
           risk.table.y.text=T,
           data = rdw.vars,
           ggtheme = theme_pubr()+theme(plot.title = element_text(hjust = 0.5)))

##########################################################################################
## (4) Single variable analysis                                                         ##
##########################################################################################

## Single variable Cox PH analysis for 30-day mortality for candidate variables ##
## Analysis of variables in their original form and original scale ##

summary(coxph(S~rdw.vars$rdw)) 
summary(coxph(S~rdw.vars$nela_risk))
summary(coxph(S~rdw.vars$hb))
summary(coxph(S~rdw.vars$cr))
summary(coxph(S~rdw.vars$indc_class))


## Single variable Cox PH analysis for overall study mortality for candidate variables ##
## Analysis of variables in their original form and original scale ##

summary(coxph(S1~rdw.vars$rdw)) 
summary(coxph(S1~rdw.vars$nela_risk))
summary(coxph(S1~rdw.vars$hb))
summary(coxph(S1~rdw.vars$cr))s
summary(coxph(S1~rdw.vars$indc_class))

##########################################################################################
## (5) Checking for outliers in candidate variables                                     ##
##########################################################################################

## Screening for outliers by examining boxplots for continous variables ##

par(mfrow=c(2,2))
boxplot(rdw.vars$rdw)
title("RDW %")
boxplot(rdw.vars$hb)
title("Hb g/l")
boxplot(rdw.vars$nela_risk)
title("NELA risk score")
boxplot(rdw.vars$cr)
title("Cr µmol/l")

## Winsorising at the 1st & 99th percentile ##

rdw.vars<-rdw.vars %>% 
  mutate(rdw=Winsorize(rdw,probs = c(0.00,0.95)),
         hb=Winsorize(hb,probs = c(0.05,1.00)),
         nela_risk=Winsorize(nela_risk, probs = c(0.00,0.95)),
         cr=Winsorize(cr,probs = c(0.00,0.95)))

##########################################################################################
## (6) Checking for linearity & transformation of variables                             ##
##########################################################################################

## shorten rdw.vars & attach as new dataframe ##

data<-rdw.vars %>% 
  dplyr::select(time30, mort30days, time, mort_stat3, rdw,nela_risk,hb, cr, sex, indc_class)
attach(data)

######################### Transformation for laparotomy indication #######################

levels(data$indc_class)[levels(data$indc_class)%in% # combine colitis, ischaemia & other
                          c("colitis","ischaemia","other")]<-"other"

############################# Transformations for rdw ####################################

dd<-datadist(data)
describe(rdw) # range is obtained

########### Variable fit for 30-day mortality ###########

# linear model 
rdwfit.linear<-cph(S~rdw, data=data)
plot(Predict(rdwfit.linear, rdw=seq(11.7,19.5, by=1)),
         xlab="RDW %, linear fit", anova=anova(rdwfit.linear),
         pval=T, data= llist(rdw))
anova(rdwfit.linear)

# dichotomous model at the median
data$rdw.dichotomised<-ifelse(rdw<=13.9,0,1)
dd<-datadist(data)
options(datadist = "dd")
rdwfit.dichotomous<-cph(S~rdw.dichotomised, data=data)
plot(Predict(rdwfit.dichotomous), xlab="Dichotomised RDW at 13.9% [median]", anova=anova(rdwfit.dichotomous), pval=T)
anova(rdwfit.dichotomous)

# categorical model 
data$rdw.quartiles<-ifelse(rdw<13.1,1,ifelse(rdw<13.9,2,ifelse(rdw<15.1,3,4)))
dd<-datadist(data)
options(datadist = "dd")
data$rdw.quartiles<-as.factor(data$rdw.quartiles)
rdwfit.quartiles<-cph(S~rdw.quartiles, data = data)
plot(Predict(rdwfit.quartiles), anova=anova(rdwfit.quartiles), xlab = "RDW quartiles", pval=T)
anova(rdwfit.quartiles)

# Log fit
rdwfit.log<-cph(S~log(rdw+1), data = data)
plot(Predict(rdwfit.log,rdw=seq(11,20, by=1)), anova=
           anova(rdwfit.log), pval=T, xlab="RDW %", main="Logarithmic fit for 30-day mortality",data=llist(rdw))
anova(rdwfit.log)

# Two-degree polynomial 
rdwfit.poly<-cph(S~pol(rdw,2), data = data)
plot(Predict(rdwfit.poly,rdw=seq(11,20, by=1)), 
         anova=anova(rdwfit.poly),pval=T, xlab="RDW %, quadratic fit",data=llist(rdw))
anova(rdwfit.poly)

# Restricted cubic spline 
rdwfit.spline<-cph(S~rcs(rdw,3), data=data)
plot(Predict(rdwfit.spline, rdw=seq(11,20, by=1)),
         anova=anova(rdwfit.spline), pval=T, xlab = "RDW %, spline fit", data=llist(rdw))
anova(rdwfit.spline)


########### Variable fit for overall study period ###########

# linear model 
rdwfit.linear<-cph(S1~rdw, data=data)
plot(Predict(rdwfit.linear, rdw=seq(11,20, by=1)),
         xlab="RDW %, linear fit", anova=anova(rdwfit.linear),
         pval=T, data= llist(rdw))
anova(rdwfit.linear)

# dichotomous model at the median
data$rdw.dichotomised<-ifelse(rdw<=13.9,0,1)
dd<-datadist(data)
options(datadist = "dd")
rdwfit.dichotomous<-cph(S1~rdw.dichotomised, data=data)
plot(Predict(rdwfit.dichotomous), xlab="Dichotomised RDW at 13.9% [median]", anova=anova(rdwfit.dichotomous), pval=T)
anova(rdwfit.dichotomous)

# categorical model 
data$rdw.quartiles<-ifelse(rdw<13.1,1,ifelse(rdw<13.9,2,ifelse(rdw<15.1,3,4)))
dd<-datadist(data)
options(datadist = "dd")
data$rdw.quartiles<-as.factor(data$rdw.quartiles)
rdwfit.quartiles<-cph(S1~rdw.quartiles, data = data)
plot(Predict(rdwfit.quartiles), anova=anova(rdwfit.quartiles), xlab = "RDW quartiles", pval=T)
anova(rdwfit.quartiles)

# Log fit
rdwfit.log<-cph(S1~log(rdw+1), data = data)
plot(Predict(rdwfit.log,rdw=seq(11,20, by=1)), anova=
            anova(rdwfit.log), pval=T, xlab="RDW %",main="Logarithmic fit for overall mortality",data=llist(rdw))
anova(rdwfit.log)

# Two-degree polynomial 
rdwfit.poly<-cph(S1~pol(rdw,2), data = data)
plot(Predict(rdwfit.poly,rdw=seq(11,20, by=1)), 
          anova=anova(rdwfit.poly),pval=T, xlab="RDW %, quadratic fit",data=llist(rdw))
anova(rdwfit.poly)

# Restricted cubic spline 
rdwfit.spline<-cph(S1~rcs(rdw,3), data=data)
plot(Predict(rdwfit.spline, rdw=seq(11,20, by=1)),
          anova=anova(rdwfit.spline), pval=T, xlab = "RDW %, spline fit", data=llist(rdw))
anova(rdwfit.spline)

############################# Transformations for nela-risk ##################################

dd<-datadist(data)
describe(nela_risk) # range is obtained

########### Variable fit for 30-day mortality ###########

# linear model 
nelafit.linear<-cph(S~nela_risk, data=data)
plot(Predict(nelafit.linear, nela_risk=seq(0.1,42, by=1)),
          xlab="NELA risk score, linear fit", anova=anova(nelafit.linear),
          pval=T, data= llist(nela_risk))
anova(nelafit.linear)

# dichotomous model at the median
data$nela.dichotomised<-ifelse(data$nela_risk<=2.5,0,1)
dd<-datadist(data)
options(datadist = "dd")
nelafit.dichotomous<-cph(S~nela.dichotomised, data=data)
plot(Predict(nelafit.dichotomous),
          anova=anova(nelafit.dichotomous), pval=T)
anova(nelafit.dichotomous)

# Categorical model 
data$nela.quartiles<-ifelse(nela_risk<0.7,1,ifelse(nela_risk<2.5,2,ifelse(nela_risk<8.8,3,4)))
dd<-datadist(data)
options(datadist = "dd")
data$nela.quartiles<-as.factor(data$nela.quartiles)
nelafit.categorical<-cph(S~nela.quartiles, data = data)
plot(Predict(nelafit.categorical), anova=anova(nelafit.categorical), xlab = "NELA risk score quartiles", pval=T)
anova(nelafit.categorical)

# Log fit
nelafit.log<-cph(S~log(nela_risk+1), data = data)
plot(Predict(nelafit.log,nela_risk=seq(0.1,42, by=1)), anova=
            anova(nelafit.log), pval=T, xlab="NELA risk score, log fit",data=llist(nela_risk))
anova(nelafit.log)

# Two-degree polynomial 
nelafit.poly<-cph(S~pol(nela_risk,2), data = data)
plot(Predict(nelafit.poly,nela_risk=seq(0.1,42, by=1)), 
          anova=anova(nelafit.poly),pval=T, xlab="NELA risk score, quadratic fit",data=llist(nela_risk))
anova(nelafit.poly)

# Restricted cubic spline 
nelafit.spline<-cph(S~rcs(nela_risk,3), data=data)
plot(Predict(nelafit.spline,nela_risk=seq(0.1,42, by=1)),
          anova=anova(nelafit.spline), pval=T, xlab = "NELA risk score, spline fit", data=llist(nela_risk))
anova(nelafit.spline)


########### Variable fit for overall study period ###########

dd<-datadist(data)
describe(nela_risk) # range is obtained

# linear model 
nelafit.linear<-cph(S1~nela_risk, data=data)
plot(Predict(nelafit.linear, nela_risk=seq(0.1,42, by=1)),
          xlab="NELA risk score, linear fit", anova=anova(nelafit.linear),
          pval=T, data= llist(nela_risk))
anova(nelafit.linear)

# dichotomous model at the median
data$nela.dichotomised<-ifelse(data$nela_risk<=2.5,0,1)
dd<-datadist(data)
options(datadist = "dd")
nelafit.dichotomous<-cph(S1~nela.dichotomised, data=data)
plot(Predict(nelafit.dichotomous), xlab="Dichotomised nela 30 days mortality risk at 13.9%", 
          anova=anova(nelafit.dichotomous), pval=T)
anova(nelafit.dichotomous)

# Categorical model 
data$nela.quartiles<-ifelse(nela_risk<0.7,1,ifelse(nela_risk<2.5,2,ifelse(nela_risk<8.8,3,4)))
dd<-datadist(data)
options(datadist = "dd")
data$nela.quartiles<-as.factor(data$nela.quartiles)
nelafit.categorical<-cph(S1~nela.quartiles, data = data)
plot(Predict(nelafit.categorical), anova=anova(nelafit.categorical), xlab = "NELA quartiles", pval=T)
anova(nelafit.categorical)

# Log fit
nelafit.log<-cph(S1~log(nela_risk+1), data = data)
plot(Predict(nelafit.log,nela_risk=seq(0.1,42, by=1)), anova=
            anova(nelafit.log), pval=T, xlab="NELA risk score, log fit",data=llist(nela_risk))
anova(nelafit.log)

# Two-degree polynomial 
nelafit.poly<-cph(S1~pol(nela_risk,2), data = data)
plot(Predict(nelafit.poly,nela_risk=seq(0.1,42, by=1)), 
          anova=anova(nelafit.poly),pval=T, xlab="NELA risk score, quadratic fit",data=llist(nela_risk))
anova(nelafit.poly)

# Restricted cubic spline 
nelafit.spline<-cph(S1~rcs(nela_risk,3), data=data)
plot(Predict(nelafit.spline,nela_risk=seq(0.1,42, by=1)),
          anova=anova(nelafit.spline), pval=T, xlab = "NELA risk score, spline fit ", data=llist(nela_risk))
anova(nelafit.spline)

############################# Transformations for hb ##################################

dd<-datadist(data)
describe(hb) # range is obtained

## Variable fit for 30-day mortality ##

# linear model 
hbfit.linear<-cph(S~hb, data=data)
plot(Predict(hbfit.linear, hb=seq(86,182, by=1)),
          xlab="Haemoglobin g/l, linear fit", anova=anova(hbfit.linear),
          pval=T, data= llist(hb))
anova(hbfit.linear)

# dichotomous model 
data$anaemia<-ifelse(sex=="M", ifelse(hb<130, "yes", "no"), ifelse(hb<120, "yes", "no"))
dd<-datadist(data)
options(datadist = "dd")
hbfit.dichotomous<-cph(S~anaemia, data=data)
plot(Predict(hbfit.dichotomous), xlab="Anaemia", anova=anova(hbfit.dichotomous), pval=T)
anova(hbfit.dichotomous)

# categorical model 
data$hbcat<-ifelse(hb<117,1,ifelse(hb<133,2,ifelse(hb<146.2,3,4)))                                                                                                                             
dd<-datadist(data)
options(datadist = "dd")
data$hbcat<-as.factor(data$hbcat)
hbfit.categorical<-cph(S~hbcat, data = data)
plot(Predict(hbfit.categorical), anova=anova(hbfit.categorical), xlab = "Haemoglobin g/l", pval=T)
anova(hbfit.categorical)

# Log fit
hbfit.log<-cph(S~log(hb+1), data = data)
plot(Predict(hbfit.log, hb=seq(86,182, by=1)), anova=
            anova(hbfit.log), pval=T, xlab="Haemoglobin g/l, log fit",data=llist(hb))
anova(hbfit.log)

# Two-degree polynomial 
hbfit.poly<-cph(S~pol(hb,2), data = data)
plot(Predict(hbfit.poly,hb=seq(86,182, by=1)), 
          anova=anova(hbfit.poly),pval=T, xlab="Haemoglobin g/l, quadratic fit",data=llist(hb))
anova(hbfit.poly)

# Restricted cubic spline 
hbfit.spline<-cph(S~rcs(hb,3), data=data)
plot(Predict(hbfit.spline, hb=seq(86,182, by=1)),
          anova=anova(hbfit.spline), pval=T, xlab = "Haemoglobin g/l, spline fit", data=llist(hb))
anova(hbfit.spline)

########### Variable fit for overall study period ###########

dd<-datadist(data)
describe(hb) # range is obtained

# linear model 
hbfit.linear<-cph(S1~hb, data=data)
plot(Predict(hbfit.linear, hb=seq(86,182, by=1)),
          xlab="Haemoglobin g/l, linear fit", anova=anova(hbfit.linear),
          pval=T, data= llist(hb))
anova(hbfit.linear)

# dichotomous model 
data$anaemia<-ifelse(sex=="M", ifelse(hb<130, "yes", "no"), ifelse(hb<120, "yes", "no"))
dd<-datadist(data)
options(datadist = "dd")
hbfit.dichotomous<-cph(S1~anaemia, data=data)
plot(Predict(hbfit.dichotomous), xlab="Anaemia", anova=anova(hbfit.dichotomous), pval=T)
anova(hbfit.dichotomous)

# categorical model 
quantile(data$hb) 
data$hbcat<-ifelse(hb<117,1,ifelse(hb<133,2,ifelse(hb<146.25,3,ifelse(hb<161,4,5))))                                                                                                                             
dd<-datadist(data)
options(datadist = "dd")
data$hbcat<-as.factor(data$hbcat)
hbfit.categorical<-cph(S1~hbcat, data = data)
plot(Predict(hbfit.categorical), anova=anova(hbfit.categorical), xlab = "Haemoglobin g/l", pval=T)
anova(hbfit.categorical)

# Log fit
hbfit.log<-cph(S1~log(hb+1), data = data)
plot(Predict(hbfit.log, hb=seq(86,182, by=1)), anova=
            anova(hbfit.log), pval=T, xlab="Haemoglobin g/l, log fit",data=llist(hb))
anova(hbfit.log)

# Two-degree polynomial 
hbfit.poly<-cph(S1~pol(hb,2), data = data)
plot(Predict(hbfit.poly,hb=seq(86,182, by=1)), 
          anova=anova(hbfit.poly),pval=T, xlab="Haemoglobin g/l, quadratic fit",data=llist(hb))
anova(hbfit.poly)

# Restricted cubic spline 
hbfit.spline<-cph(S1~rcs(hb,3), data=data)
plot(Predict(hbfit.spline, hb=seq(86,182, by=1)),
          anova=anova(hbfit.spline), pval=T, xlab = "Haemoglobin g/l, spline fit", data=llist(hb))
anova(hbfit.spline)


############################# Transformations for cr  ##################################

dd<-datadist(data)
describe(cr) # range is obtained

########### Variable fit for 30-day mortality ###########

# linear model 
crfit.linear<-cph(S~cr, data=data)
plot(Predict(crfit.linear, cr=seq(41,203, by=10)),
          xlab="Cr µmol/l, linear fit", anova=anova(crfit.linear),
          pval=T, data= llist(cr))
anova(crfit.linear)

# dichotomous model at the median
data$cr.dichotomised<-ifelse(cr<=75.5,0,1)
dd<-datadist(data)
options(datadist = "dd")
crfit.dichotomous<-cph(S~cr.dichotomised, data=data)
plot(Predict(crfit.dichotomous), xlab="Dichotomised Cr at 75.5", anova=anova(crfit.dichotomous), pval=T)
anova(crfit.dichotomous)

# categorical model 
data$cr.cat<-ifelse(cr<65,1,ifelse(cr<75.5,2,ifelse(cr<102,3,4)))
dd<-datadist(data)
  options(datadist = "dd")
  data$cr.cat<-as.factor(data$cr.cat)
  crfit.categorical<-cph(S~cr.cat, data = data)
  plot(Predict(crfit.categorical), anova=anova(crfit.categorical), xlab = "Cr", pval=T)
  anova(crfit.categorical)
  
  # Log fit
  crfit.log<-cph(S~log(cr+1), data = data)
  plot(Predict(crfit.log,cr=seq(41,203, by=10)), anova=
              anova(crfit.log), pval=T, xlab="Cr µmol/l, log fit",data=llist(cr))
  anova(crfit.log)
  
  # Two-degree polynomial 
  crfit.poly<-cph(S~pol(cr,2), data = data)
  plot(Predict(crfit.poly,cr=seq(52,203, by=10)), 
            anova=anova(crfit.poly),pval=T, xlab="Cr µmol/l, quadratic fit",data=llist(cr))
  anova(crfit.poly)
  
  # Restricted cubic spline 
  crfit.splines<-cph(S~rcs(cr,3),data=data)
  plot(Predict(crfit.splines, cr=seq(52,203, by=10)),
            anova=anova(crfit.splines), pval=T, xlab = "Cr µmol/l, spline fit", data=llist(cr))
  anova(crfit.splines)
  
  ########### Variable fit for overall study period ###########
  
  dd<-datadist(data)
  describe(cr) # range is obtained
  
  # linear model 
  crfit.linear<-cph(S1~cr, data=data)
  plot(Predict(crfit.linear, cr=seq(41,203, by=10)),
            xlab="Cr µmol/l, linear fit", anova=anova(crfit.linear),
            pval=T, data= llist(cr))
  anova(crfit.linear)
  
  # dichotomous model at the median
  data$cr.dichotomised<-ifelse(cr<=75.5,0,1)
  dd<-datadist(data)
  options(datadist = "dd")
  crfit.dichotomous<-cph(S1~cr.dichotomised, data=data)
  plot(Predict(crfit.dichotomous), xlab="Dichotomised Cr at 75.5", anova=anova(crfit.dichotomous), pval=T)
  anova(crfit.dichotomous)
  
  # categorical model 
  data$cr.cat<-ifelse(cr<65,1,ifelse(cr<75.5,2,ifelse(cr<102,3,4)))
  dd<-datadist(data)
  options(datadist = "dd")
  data$cr.cat<-as.factor(data$cr.cat)
  crfit.categorical<-cph(S1~cr.cat, data = data)
  plot(Predict(crfit.categorical), anova=anova(crfit.categorical), xlab = "Cr", pval=T)
  anova(crfit.categorical)
  
  # Log fit
  crfit.log<-cph(S1~log(cr+1), data = data)
  plot(Predict(crfit.log,cr=seq(52,203, by=10)), anova=
              anova(crfit.log), pval=T, xlab="Cr µmol/l,log fit",data=llist(cr))
  anova(crfit.log)
  
  # Two-degree polynomial 
  crfit.poly<-cph(S1~pol(cr,2), data = data)
  plot(Predict(crfit.poly,cr=seq(41,203, by=10)), 
            anova=anova(crfit.poly),pval=T, xlab="Cr µmol/l, quadratic fit",data=llist(cr))
  anova(crfit.poly)
  
  # Restricted cubic spline 
  crfit.splines<-cph(S1~rcs(cr,3),data=data)
  plot(Predict(crfit.splines, cr=seq(41,203, by=10)),
            anova=anova(crfit.splines), pval=T, xlab = "Cr µmol/l, spline fit", data=llist(cr))
  anova(crfit.splines)


##########################################################################################
## (7.1) Fitting the model for 30-day mortality                                         ##
##########################################################################################

## Fitting model 30-day mortality - with RDW and without

model1<-cph(S~nela_risk+rcs(hb,3)+cr+indc_class,x=T, y=T, surv=T, time.inc = 30, data=data)
model2<-cph(S~log(rdw+1)+nela_risk+rcs(hb,3)+cr+indc_class,x=T, y=T, surv=T, time.inc = 30, data=data)

##########################################################################################
## (7.2) Fitting the model for overall mortality                                        ##
##########################################################################################

## Fitting model overall mortality - with RDW and without 

model3<-cph(S1~nela_risk+rcs(hb,3)+cr+indc_class,x=T, y=T, surv=T, time.inc =1580/365.25 , data=data)
model4<-cph(S1~log(rdw+1)+nela_risk+rcs(hb,3)+cr+indc_class,x=T, y=T, surv=T, time.inc =1580/365.25, 
            data=data)

##########################################################################################
## (8) Checking for multicollinearity                                                   ##
##########################################################################################

names(data)
data.cor<-data[,-c(1:4,9:18)] # removing non-numeric variables from dataframe
cor(data.cor)
vif(model1)
vif(model2)
vif(model3)
vif(model4)

##########################################################################################
## (9.1) Checking for interactions - 30-mortality model                                ##
##########################################################################################

## global tests for the model without rdw ##

z1<-predict(model1, type = "terms")

nela.ia<-z1[,"nela_risk"]
all.others<-z1[,-1]
anova(cph(S~nela.ia*all.others))
cph(S~nela.ia*all.others)

hb.ia<-z1[,"hb"]
all.others<-z1[,-2]
anova(cph(S~hb.ia*all.others))
cph(S~hb.ia*all.others)

cr.ia<-z1[,"cr"]
all.others<-z1[,-3]
anova(cph(S~cr.ia*all.others))
cph(S~cr.ia*all.others)

indc.ia<-z1[,"indc_class"]
all.others<-z1[,-4]
anova(cph(S~indc.ia*all.others))
cph(S~indc.ia*all.others)


## checking the full model for interactions ##

z2<-predict(model2, type="terms")

rdw.ia<-z2[,"rdw"]
all.others.2<-z2[,-1]
anova(cph(S~rdw.ia*all.others.2))
cph(S~rdw.ia*all.others.2)

nela.ia2<-z2[,"nela_risk"]
all.others.2<-z2[,-2]
anova(cph(S~nela.ia2*all.others.2))
cph(S~nela.ia2*all.others.2)

hb.ia2<-z2[,"hb"]
all.others.2<-z2[,-3]
anova(cph(S~hb.ia2*all.others.2))
cph(S~hb.ia2*all.others.2)

cr.ia2<-z2[,"cr"]
all.others.2<-z2[,-4]
anova(cph(S~cr.ia2*all.others.2))
cph(S~cr.ia2*all.others.2)

indc.ia2<-z2[,"indc_class"]
all.others.2<-z2[,-5]
anova(cph(S~indc.ia2*all.others.2))
cph(S~indc.ia2*all.others.2)

##########################################################################################
## (9.2) Checking for interactions - overall mortality model                           ##
##########################################################################################

## global tests for the model without rdw ##

z3<-predict(model3, type = "terms")

nela.ia3<-z3[,"nela_risk"]
all.others.3<-z3[,-1]
anova(cph(S1~nela.ia3*all.others.3))
cph(S1~nela.ia3*all.others.3)

hb.ia3<-z3[,"hb"]
all.others.3<-z3[,-2]
anova(cph(S1~hb.ia3*all.others.3))
cph(S1~hb.ia3*all.others.3)

cr.ia3<-z3[,"cr"]
all.others.3<-z3[,-3]
anova(cph(S1~cr.ia3*all.others.3))
cph(S1~cr.ia3*all.others.3)

indc.ia3<-z3[,"indc_class"]
all.others.3<-z3[,-4]
anova(cph(S1~indc.ia3*all.others.3))
cph(S1~indc.ia3*all.others.3)

## checking the full model for interactions ##

z4<-predict(model4, type="terms")

rdw.ia4<-z4[,"rdw"]
all.others.4<-z4[,-1]
anova(cph(S1~rdw.ia4*all.others.4))
cph(S1~rdw.ia4*all.others.4)

nela.ia4<-z4[,"nela_risk"]
all.others.4<-z4[,-2]
anova(cph(S1~nela.ia4*all.others.4))
cph(S1~nela.ia4*all.others.4)

hb.ia4<-z4[,"hb"]
all.others.4<-z4[,-3]
anova(cph(S1~hb.ia4*all.others.4))
cph(S1~hb.ia4*all.others.4)

cr.ia4<-z4[,"cr"]
all.others.4<-z4[,-4]
anova(cph(S1~cr.ia4*all.others.4))
cph(S1~cr.ia4*all.others.4)

indc.ia4<-z4[,"indc_class"]
all.others.4<-z4[,-5]
anova(cph(S1~indc.ia4*all.others.4))
cph(S1~indc.ia4*all.others.4)

##########################################################################################
## (10) Updating models with interactions                                               ##
##########################################################################################

## Updating the 30-day mortality model with interactions ##

model1<-cph(S~nela_risk+rcs(hb,3)*cr+indc_class,x=T, y=T, surv=T, time.inc = 30, data=data)
model2<-cph(S~log(rdw+1)+nela_risk+rcs(hb,3)*cr+indc_class,x=T, y=T, surv=T, time.inc = 30, data=data)

## Updating the overall mortality model with interactions ##

model3<-cph(S1~nela_risk*cr+indc_class+rcs(hb,3)*cr,x=T, y=T, surv=T, time.inc =1580/365.25 , data=data)
model4<-cph(S1~log(rdw+1)+nela_risk*cr+rcs(hb,3)*cr+indc_class,x=T, y=T, surv=T, time.inc =1580/365.25, 
            data=data)
summary(model4)

##########################################################################################
## (11) Checking proportional hazards assumption                                        ##
##########################################################################################

## Checking PH assumption in 30-day mortality model  ##
z5<-predict(model1, type="terms")
model1.short<-cph(S~z5, x=T, y=T)
ph1<-cox.zph(model1.short, transform = "identity")
ph1
plot(ph1)

z6<-predict(model2, type="terms")
model2.short<-cph(S1~z6, x=T, y=T)
ph2<-cox.zph(model2.short, transform = "identity")
ph2
plot(ph2)

## Checking PH assumption in overall mortality model ##

z7<-predict(model3, type="terms")
model3.short<-cph(S~z7, x=T, y=T)
ph3<-cox.zph(model3.short, transform = "identity")
ph3
plot(ph3)

z8<-predict(model4, type="terms")
model4.short<-cph(S1~z8, x=T, y=T)
ph4<-cox.zph(model4.short, transform = "identity")
ph4
plot(ph4)

##########################################################################################
## (12) Influential observations                                                        ##
##########################################################################################

inf1<-which.influence(model1)
inf2<-which.influence(model2)
inf3<-which.influence(model3)
inf4<-which.influence(model4)

## Sensitivity analysis without influential values for RDW model2 ##

subset<-data[-c(36,98,157,202,296),]
attach(subset)

S.sens<-Surv(time30, mort30days)
sensitivity.model<-cph(S.sens~log(rdw+1)+nela_risk+rcs(hb,3)*cr+indc_class,x=T, y=T, surv=T, time.inc = 30, data=subset)
sensitivity.model
anova(sensitivity.model)

anova(model2) #compare the models with and without influential outliers 
detach(subset)

## Sensitivity analysis without influential values for RDW model4 ##

subset1<-data[-c(51,52,98,163),]
attach(subset1)

S.sens1<-Surv(time/365.25, mort_stat3)
sensitivity.model1<-cph(S.sens1~log(rdw+1)+nela_risk*cr+rcs(hb,3)*cr+indc_class,x=T, y=T, surv=T, time.inc =1580/365.25, 
                        data=subset1)
sensitivity.model1
anova(sensitivity.model1)

anova(model4)  #compare the models with and without influential outliers 
detach(subset1)

##########################################################################################
## (13) Relative contribution of RDW                                                    ##
##########################################################################################

## 30-day mortality model ##
plot(anova(model1), margin = "P", rm.ia = TRUE)
title(main="30-day mortality Cox regression model 
      without RDW", cex.main=0.8)

plot(Predict(model1), anova=anova(model1), pval=T,
     main="30-day mortality Cox regression model without RDW")

plot(anova(model2), margin="P", rm.ia = TRUE)
title(main="30-day mortality Cox regression model 
      with RDW", cex.main=0.8)

plot(Predict(model2), anova = anova(model2), pval = TRUE,
     main="30-day mortality Cox regression model with RDW")


## Overall mortality model ##
plot(anova(model3), margin="P", rm.ia = TRUE)
title(main="Overall mortality Cox regression model 
      without RDW",cex.main=0.8)

plot(Predict(model3), anova = anova(model3), pval = TRUE,
     main="Overall mortality Cox regression model without RDW")

plot(anova(model4), margin="P", rm.ia = TRUE)
title(main="Overall mortality Cox regression model 
      with RDW", cex.main=0.8)

plot(Predict(model4), anova = anova(model4), pval = TRUE,
     main="Overall mortality Cox regression model with RDW")


##########################################################################################
## (14) Added value of RDW                                                              ##
##########################################################################################

## Likelihood ratio chi square test for 30-day mortality model ##
lrtest(model1,model2)

## Likelihood ratio chi square test for overall mortality model ##
lrtest(model3,model4)

##########################################################################################
## (15) Internal validation                                                             ##
##########################################################################################

set.seed(123)

validate(model3, B=1000)
validate(model4, B=1000)

##########################################################################################
## (16) Calibration plots - not adjusted for overfitting                                ##
##########################################################################################

## calibration of model 4 - overall study period

units(data$time)<-"Year"
model4<-cph(S1~log(rdw+1)+nela_risk*cr+rcs(hb,3)*cr+indc_class,x=T, y=T, surv=T, time.inc =1580/365.25, 
            data=data)

set.seed(717) # to ensure reproducibility
cal2<-calibrate(model4, u=4,B=300, cmethod = "hare")
plot(cal2,subtitles=FALSE)
title("Fig. 5 Risk prediction calibration plot for overall survival model", cex.main=1.1)

##########################################################################################
## (17) Determining heuristic shrinkage                                                 ##
##########################################################################################

## Calculating heuristic shrinkage based on LR chi2 AIC ##

model4 # LR chi2 79.89, d.f. 11
(79.89-11)/79.89 # 0.8623107

##########################################################################################
## (18) LASSO shrinkage to adjust for overfitting                                       ##
##########################################################################################

library(glmpath)

########## Formal determination of shrinkage factor - for overall model #################

## determining shrinkage factor by AIC ## 

mydata1<-list(x=predict(model4, type="ccterms"), time=data$time, status=data$mort_stat3)
path1<-coxpath(data = mydata1)
plot(path1, mar = c(5,4,4,10)-0.5)
plot(path1, type="aic")

## Shrinkage factor for overall model ## 
lasso.factors1<-path1$b.predictor[path1$aic==min(path1$aic),] # shrinkage factor

## Shrinking the model coefficients ## 

lasso.coefs1<-model4$coef

lasso.coefs1["rdw"]<-lasso.coefs1["rdw"]*lasso.factors1[1]
lasso.coefs1["nela_risk"]<-lasso.coefs1["nela_risk"]*lasso.factors1[2]
lasso.coefs1["cr"]<-lasso.coefs1["cr"]*lasso.factors1[2]
lasso.coefs1["hb"]<-lasso.coefs1["hb"]*lasso.factors1[2]
lasso.coefs1["hb'"]<-lasso.coefs1["hb'"]*lasso.factors1[2]
lasso.coefs1["indc_class=other"]<-lasso.coefs1["indc_class=other"]*lasso.factors1[3]
lasso.coefs1["indc_class=obstruction"]<-lasso.coefs1["indc_class=obstruction"]*lasso.factors1[3]
lasso.coefs1["indc_class=sepsis"]<-lasso.coefs1["indc_class=sepsis"]*lasso.factors1[3]
lasso.coefs1[9]<-lasso.coefs1[9]*lasso.factors1[2]*lasso.factors1[2]
lasso.coefs1[10]<-lasso.coefs1[10]*lasso.factors1[2]*lasso.factors1[2]
lasso.coefs1[11]<-lasso.coefs1[11]*lasso.factors1[2]*lasso.factors1[2]

## Updating the model coefficients ## 

lassomodel1<-model4
lassomodel1$coefficients<-lasso.coefs1

##########################################################################################
## (19) Forest plot & nomogram                                                        ##
##########################################################################################

plot(summary(lassomodel1,rdw=c(13.9,17.4),hb=c(85,133), cr=c(75,150)), at=c(seq(-1,10,by=1)), q=c(0.95),
     col = "darkgrey", col.points = "black",pch = 3, lwd = 1.5,
     main = "Fig. 6 Hazard ratio - Overall mortality model")

## nomogram for overall mortality model ## 

surv1<-Survival(lassomodel1)
surv.total<-function(x) surv1(4.325804*0.9246836, lp=x)
ss<-c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95)

nomo1<-nomogram(lassomodel1, rdw=c(seq(10,20, by=0.5)), nela_risk=c(seq(0,100,by=100)),
                cr=c(80,100,120,140,160,180,200), hb=c(80,130), interact = list(nela_risk=c(0.7,2.5,8.8),
                                                                                hb=c(80,130),cr),
                fun=surv.total, funlabel = "Predicted 4 year survival", fun.at = ss,
                lp=T, lp.at = c(seq(-2,4, by=0.5)), nint = 10, maxscale = 100)

plot(nomo1, xfrac = .65, lmgp = .35, cex.axis = 0.65)
title(main="Fig. 7 Nomogram predicting 4-year survival post emergency laparotomy", cex.main=1)






