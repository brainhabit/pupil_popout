#####################################################################################################################################################################################################
#
# This script is used to analyse data for the infant pupillometry/eye-tracking study based on a pop-out visual attention task.
# 
# Giorgia Bussu, 2023
#
######################################################################################################################################################################################################

rm(list=ls())

setwd('C:/Users/your_name/Documents/pupil_popout')
list.files()

library(readxl)
library(glmnet)

use_imputed <- 1

if(use_imputed==1){load("data_imputation_final.RData")}

if(use_imputed==0){

data<-read_xlsx('data/pupil_sound_dataset_long_original.xlsx')

# missing gaze (-1) set to NaN

data$inAOILatency_1[which(data$inAOILatency_1==-1)]<-NA
data$inAOILatency_2[which(data$inAOILatency_2==-1)]<-NA
data$inAOILatency_3[which(data$inAOILatency_3==-1)]<-NA
data$inAOILatency_4[which(data$inAOILatency_4==-1)]<-NA

data$inAOIRawLatency_1[which(data$inAOIRawLatency_1==-1)]<-NA
data$inAOIRawLatency_2[which(data$inAOIRawLatency_2==-1)]<-NA
data$inAOIRawLatency_3[which(data$inAOIRawLatency_3==-1)]<-NA
data$inAOIRawLatency_4[which(data$inAOIRawLatency_4==-1)]<-NA

data_selected<-data

validgaze<-which(data_selected$ValidGaze==0)

data_selected$latencycheck<-data_selected$MinLatency
data_selected$MinLatency[validgaze]<-NA

data_selected$firstlookcheck<-data_selected$FirstLookFace
data_selected$FirstLookFace[validgaze]<-NA

data_selected$pupilcheck<-data_selected$PupilDilation
data_selected$PupilDilation[validgaze]<-NA
data_selected$MeanBaseline[validgaze]<-NA
data_selected$MeanResponse[validgaze]<-NA

data_selected$PupilDilation<-as.numeric(data_selected$PupilDilation)
data_selected$FirstLookFace<-as.numeric(data_selected$FirstLookFace)
data_selected$MinLatency<-as.numeric(data_selected$MinLatency)

data_selected$firstlookcheck<-as.numeric(data_selected$firstlookcheck)
data_selected$pupilcheck<-as.numeric(data_selected$pupilcheck)
data_selected$latencycheck<-as.numeric(data_selected$latencycheck)

id_name<-unique(data_selected$id_name)

prop_invalid_minlat<-rep(0,length(id_name))
for(ii in 1:length(id_name)){
  prop_invalid_minlat[ii]<-length(which(is.na(data_selected$MinLatency[which(data_selected$id_name==id_name[ii])])))/length(data_selected$MinLatency[which(data_selected$id_name==id_name[ii])])
}
excessive_invalid_minlat<-which(prop_invalid_minlat>0.3)

prop_invalid_face<-rep(0,length(id_name))
for(ii in 1:length(id_name)){
  prop_invalid_face[ii]<-length(which(is.na(data_selected$FirstLookFace[which(data_selected$id_name==id_name[ii])])))/length(data_selected$FirstLookFace[which(data_selected$id_name==id_name[ii])])
}
excessive_invalid_face<-which(prop_invalid_face>0.3)

prop_invalid_pupil<-rep(0,length(id_name))
for(ii in 1:length(id_name)){
  prop_invalid_pupil[ii]<-length(which(is.na(data_selected$PupilDilation[which(data_selected$id_name==id_name[ii])])))/length(data_selected$PupilDilation[which(data_selected$id_name==id_name[ii])])
}
excessive_invalid_pupil<-which(prop_invalid_pupil>0.3)


excluded_id<-id_name[excessive_invalid_pupil]
excluded_id<-substr(excluded_id,1,3)
substr(excluded_id[which(substr(excluded_id,1,1)=='p')],1,1)<-'P'


for(ii in 1:length(excessive_invalid_pupil)){
  data_selected<-data_selected[-which(data_selected$id_name==id_name[excessive_invalid_pupil[ii]]),]
}



# Remove participant based on technical error reported by students
data_selected<-data_selected[-which(data_selected$id_name==id_name[41]),]



################## DATA IMPUTATION 

data_to_impute<-data_selected[-c(2,25:33)]

data_to_impute$MeanBaseline<-as.numeric(data_to_impute$MeanBaseline)
data_to_impute$MeanBaseline[which(data_to_impute$MeanBaseline=='NaN')]<-NA
data_to_impute$MeanResponse<-as.numeric(data_to_impute$MeanResponse)
data_to_impute$MeanResponse[which(data_to_impute$MeanResponse=='NaN')]<-NA
data_to_impute$PupilDilation[which(data_to_impute$PupilDilation=='NaN')]<-NA

data_to_impute$pupilcheck[which(data_to_impute$pupilcheck=='NaN')]<-NA
data_to_impute$latencycheck[which(data_to_impute$latencycheck=='NaN')]<-NA
data_to_impute$firstlookcheck[which(data_to_impute$firstlookcheck=='NaN')]<-NA
which(data_to_impute=='NaN')

# clean columns not necessary for imputation anyway (sampling rate & validity flags) or embedded in NAs of mean variables
data_to_impute<-data_to_impute[-c(8,10:17,26,27,30:32,34,36)]

# age & sex
demo<-read_xlsx('data/demo.xlsx')

## check missing data excluded participants
substr(demo$id[which(substr(demo$id,1,1)=='p')],1,1)<-'P'
demo$excluded<-rep(0,length(demo$id))
demo$excluded[match(excluded_id,demo$id)]<-1
t.test(age~excluded,demo)
demo$sex_num<-as.numeric(as.factor(demo$sex))
t.test(sex_num~excluded,demo)

# match demo and idname
id_name<-unique(data_to_impute$id_name)
id_name_str<-substr(id_name,1,3)
substr(id_name_str[which(substr(id_name_str,1,1)=='p')],1,1)<-'P'

demo_ordered<-demo[match(id_name_str,demo$id),]

data_to_impute$age<-rep(0,dim(data_to_impute)[1])
for(ii in 1:length(id_name)){
data_to_impute$age[which(data_to_impute$id_name==id_name[ii])]<-demo_ordered$age[ii]
}

data_to_impute$sex<-rep(0,dim(data_to_impute)[1])
for(ii in 1:length(id_name)){
data_to_impute$sex[which(data_to_impute$id_name==id_name[ii])]<-demo_ordered$sex[ii]
}

data_to_impute$sex[which(data_to_impute$sex=='Boy')]<-1
data_to_impute$sex[which(data_to_impute$sex=='Girl')]<-0
data_to_impute$sex<-as.numeric(data_to_impute$sex)



# imputation NB do not impute last 3 columns as those are for missing data checks (RM)
library(mice)

data_original<-data
data<-data.frame(data_to_impute)

tempData <- mice(data,m=100,maxit=50,meth='cart',seed=2022)

# pooling 
complete_data<-complete(tempData)

}

#save.image("~/pupil_sound/data_imputation_final.RData")

####################################################### ANALYSIS ###########################################################################################

id_name<-unique(complete_data$id)

ntrials<-rep(0,length(id_name));for(i in 1:length(id_name)){
ntrials[i]<-length(which(complete_data$ValidGaze[which(complete_data$id==id_name[i])]==1))
}

complete_data$validtrials<-rep(0,dim(complete_data)[1])
for(ii in 1:length(id_name)){
  complete_data$validtrials[which(complete_data$id==id_name[ii])]<-ntrials[ii]
}


library(lmerTest)

## prep variables

#robust ranges
complete_data$PupilDilation[which(complete_data$PupilDilation>(mean(complete_data$PupilDilation,na.rm=T)+3*sd(complete_data$PupilDilation,na.rm=T)))]<-NA
complete_data$PupilDilation[which(complete_data$PupilDilation<(mean(complete_data$PupilDilation,na.rm=T)-3*sd(complete_data$PupilDilation,na.rm=T)))]<-NA
complete_data$MinLatency[which(complete_data$MinLatency>(mean(complete_data$MinLatency,na.rm=T)+3*sd(complete_data$MinLatency,na.rm=T)))]<-NA
complete_data$MinLatency[which(complete_data$MinLatency<(mean(complete_data$MinLatency,na.rm=T)-3*sd(complete_data$MinLatency,na.rm=T)))]<-NA


id_name<-unique(data_selected$id_name)

missingpupil_idmean<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  missingpupil_idmean[ii]<-mean(complete_data$MissingRawPupil[which(complete_data$id==id_name[ii])])
}
complete_data$missingpupil_idmean<-rep(0,dim(complete_data)[1])
for(ii in 1:length(id_name)){
  complete_data$missingpupil_idmean[which(complete_data$id==id_name[ii])]<-missingpupil_idmean[ii]
}

complete_data$idwc_missingpupil<-complete_data$MissingRawPupil-complete_data$missingpupil_idmean

complete_data$missingpupil_idmean<-log(complete_data$missingpupil_idmean)
complete_data$idwc_missingpupil<-log(complete_data$idwc_missingpupil+0.3)

isi_idmean<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  isi_idmean[ii]<-mean(complete_data$isi[which(complete_data$id==id_name[ii])])
}
complete_data$isi_idmean<-rep(0,dim(complete_data)[1])
for(ii in 1:length(id_name)){
  complete_data$isi_idmean[which(complete_data$id==id_name[ii])]<-isi_idmean[ii]
}

complete_data$idwc_isi<-complete_data$isi-complete_data$isi_idmean


# set to same scale
complete_data$validtrials<-scale(complete_data$validtrials)
complete_data$missingpupil_idmean<-scale(complete_data$missingpupil_idmean)
complete_data$isi_idmean<-scale(complete_data$isi_idmean)


# binary
complete_data$Order<-complete_data$Order-1
complete_data$Order<-as.factor(complete_data$Order)

complete_data$age_std<-complete_data$age
complete_data$age_std[which(complete_data$age_std>5)]<-1
complete_data$age_std[which(complete_data$age_std>1)]<-0
complete_data$age_std<-as.factor(complete_data$age_std)

# 0=girl; 1=boy
complete_data$sex<-as.factor(complete_data$sex)

# categorical contrasts set to reference level (silence for type and volume)
complete_data$id<-as.factor(complete_data$id)

# 0=silent; 1=social; 2=non-social
complete_data$qualitative<-as.factor(complete_data$AudioType)
contrasts(complete_data$qualitative)<-contr.treatment(3,base=1)

complete_data$quantitative<-as.factor(complete_data$AudioVolume)
contrasts(complete_data$quantitative)<-contr.treatment(3,base=1)

complete_data$isi_std<-scale(complete_data$isi)

###################################################################################################################################################################################################
#### MODELLING PUPIL DILATION ############################################################################################################################################################################
###################################################################################################################################################################################################

# step-wise model definition
anova(lmer(PupilDilation ~ (1|id),data = complete_data),lm(PupilDilation ~ 1,data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + missingpupil_idmean + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + idwc_isi + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + isi_std + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + validtrials + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + age_std + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + sex + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + qualitative + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + qualitative + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + age_std*qualitative + (1|id),data = complete_data))
anova(lmer(PupilDilation ~ idwc_missingpupil + age_std*qualitative + (1|id),data = complete_data),lmer(PupilDilation ~ idwc_missingpupil + age_std*qualitative + quantitative + (1|id),data = complete_data))

# sound type
anova(lmer(PupilDilation ~ idwc_missingpupil + age_std*qualitative + (1|id),data = complete_data,REML=F))
summary(lmer(PupilDilation ~ idwc_missingpupil + age_std*qualitative + (1|id),data = complete_data,REML=F))

# sound volume
anova(lmer(PupilDilation ~ idwc_missingpupil + quantitative + (1|id),data = complete_data,REML=F))
summary(lmer(PupilDilation ~ idwc_missingpupil + quantitative + (1|id),data = complete_data,REML=F))

# sound repetition
anova(lmer(PupilDilation ~ idwc_missingpupil + Order + (1|id),data = complete_data,REML=F))
summary(lmer(PupilDilation ~ idwc_missingpupil + Order + (1|id),data = complete_data,REML=F))

# total model; separate model tested for main effects
pupil_model<-lmer(PupilDilation ~ idwc_missingpupil + age_std*qualitative + quantitative+ (1|id),data = complete_data,REML=F) # try with REML?

anova(pupil_model)
summary(pupil_model) 

# model diagnostics
model_pupil<-pupil_model

qqnorm(residuals(model_pupil))
qqline(residuals(model_pupil))
shapiro.test(residuals(model_pupil))

library(MuMIn)
r.squaredGLMM(model_pupil)
print(anova(model_pupil))

library(multcomp)
glht(model_pupil, linfct = mcp(qualitative="Tukey"))

performance::icc(model_pupil)

ranova(model_pupil)

plot(model_pupil)

###################################################################################################################################################################################################
#### MODELLING LATENCY ############################################################################################################################################################################
###################################################################################################################################################################################################

id_name<-unique(data_selected$id_name)

missingraw_idmean<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  missingraw_idmean[ii]<-mean(complete_data$MissingRawGazeTrial[which(complete_data$id==id_name[ii])])
}
complete_data$missinggaze_idmean<-rep(0,dim(complete_data)[1])
for(ii in 1:length(id_name)){
  complete_data$missinggaze_idmean[which(complete_data$id==id_name[ii])]<-missingraw_idmean[ii]
}

complete_data$idwc_missinggaze<-complete_data$MissingRawGazeTrial-complete_data$missinggaze_idmean

complete_data$missinggaze_idmean<-log(complete_data$missinggaze_idmean)
complete_data$idwc_missinggaze<-log(complete_data$idwc_missinggaze+0.28)

onscreen_idmean<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  onscreen_idmean[ii]<-mean(complete_data$onScreenArray[which(complete_data$id==id_name[ii])])
}
complete_data$onscreen_idmean<-rep(0,dim(complete_data)[1])
for(ii in 1:length(id_name)){
  complete_data$onscreen_idmean[which(complete_data$id==id_name[ii])]<-onscreen_idmean[ii]
}

complete_data$idwc_onscreen<-complete_data$onScreenArray-complete_data$onscreen_idmean

complete_data$onscreen_idmean<-I(complete_data$onscreen_idmean)^3

# same scale
complete_data$latency_std<-log(complete_data$MinLatency)
complete_data$missinggaze_idmean<-scale(complete_data$missinggaze_idmean)
complete_data$onscreen_idmean<-scale(complete_data$onscreen_idmean)


# step-wise model definition
anova(lmer(latency_std ~ (1|id),data = complete_data),lm(latency_std ~ 1,data = complete_data))
anova(lmer(latency_std ~ (1|id),data = complete_data),lmer(latency_std ~ age_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + sex + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + idwc_missinggaze + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + missinggaze_idmean + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + (1|id),data = complete_data),lmer(latency_std ~ age_std +idwc_onscreen + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+  (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean+ validtrials + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+ (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean+ isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+  (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean+ idwc_isi + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+ isi_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean+isi_std + qualitative + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+ isi_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean+isi_std + quantitative + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ age_std + onscreen_idmean+isi_std + Order + (1|id),data = complete_data))

anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ qualitative *quantitative + age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ Order*quantitative + age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ Order*qualitative + age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ Order*age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ qualitative*age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ quantitative*age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ quantitative*sex+ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ qualitative*sex+ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ Order*sex+ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ quantitative*isi_std+ age_std + onscreen_idmean + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ qualitative*isi_std+ age_std + onscreen_idmean + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ Order*isi_std+ age_std + onscreen_idmean + (1|id),data = complete_data))
anova(lmer(latency_std ~ age_std + onscreen_idmean+isi_std + (1|id),data = complete_data),lmer(latency_std ~ isi_std*age_std + onscreen_idmean + (1|id),data = complete_data))


# main effects
anova(lmer(latency_std ~ age_std*isi_std + onscreen_idmean + qualitative+(1|id),data = complete_data, REML=F))
anova(lmer(latency_std ~ age_std*isi_std + onscreen_idmean + quantitative+(1|id),data = complete_data, REML=F))
anova(lmer(latency_std ~ age_std*isi_std + onscreen_idmean + Order+(1|id),data = complete_data, REML=F))

# model diagnostics
model_latency<-lmer(latency_std ~ age_std*isi_std + onscreen_idmean + (1|id),data = complete_data, REML=F)

complete_data$pupil_scaled<-scale(complete_data$PupilDilation)

model_latency_pupil <- lmer(latency_std ~ age_std*isi_std + onscreen_idmean + pupil_scaled + (1|id),data = complete_data, REML=F)
 
qqnorm(residuals(model_latency))
qqline(residuals(model_latency))
shapiro.test(residuals(model_latency))

library(MuMIn)
r.squaredGLMM(model_latency)
print(anova(model_latency))

performance::icc(model_latency)

ranova(model_latency)

plot(model_latency)

###################################################################################################################################################################################################
#### MODELLING SOCIAL ORIENTING ############################################################################################################################################################################
###################################################################################################################################################################################################

#library(glmer)
pupil_idmean<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  pupil_idmean[ii]<-mean(complete_data$PupilDilation[which(complete_data$id==id_name[ii])],na.rm=T)
}

face_id_volume0a<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  face_id_volume0a[ii]<-mean(complete_data$FirstLookFace[which(complete_data$id==id_name[ii] & complete_data$AudioVolume==0 & complete_data$Order==0)],na.rm=T)
}
face_id_volume0b<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  face_id_volume0b[ii]<-mean(complete_data$FirstLookFace[which(complete_data$id==id_name[ii] & complete_data$AudioVolume==0 & complete_data$Order==1)],na.rm=T)
}

face_id_volume1a<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  face_id_volume1a[ii]<-mean(complete_data$FirstLookFace[which(complete_data$id==id_name[ii] & complete_data$AudioVolume==0.7 & complete_data$Order==0)],na.rm=T)
}
face_id_volume1b<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  face_id_volume1b[ii]<-mean(complete_data$FirstLookFace[which(complete_data$id==id_name[ii] & complete_data$AudioVolume==0.7 & complete_data$Order==1)],na.rm=T)
}


face_id_volume2a<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  face_id_volume2a[ii]<-mean(complete_data$FirstLookFace[which(complete_data$id==id_name[ii] & complete_data$AudioVolume==1 & complete_data$Order==0)],na.rm=T)
}
face_id_volume2b<-rep(0,length(id_name))

for(ii in 1:length(id_name)){
  face_id_volume2b[ii]<-mean(complete_data$FirstLookFace[which(complete_data$id==id_name[ii] & complete_data$AudioVolume==1 & complete_data$Order==1)],na.rm=T)
}

faceplot<-data.frame(cbind(rep(id_name,6),c(face_id_volume0a,face_id_volume0b,face_id_volume1a,face_id_volume1b,face_id_volume2a,face_id_volume2b),c(rep(0,46*2),rep(0.7,46*2),rep(1,46*2)),rep(c(rep(0,46),rep(1,46)),3)))


complete_data$pupil_idmean<-rep(0,dim(complete_data)[1])

for(ii in 1:length(id_name)){
  complete_data$pupil_idmean[which(complete_data$id==id_name[ii])]<-pupil_idmean[ii]
}

complete_data$idwc_pupil<-complete_data$PupilDilation-complete_data$pupil_idmean
complete_data$pupil_mean_scaled<-scale(complete_data$pupil_idmean)

complete_data$pupil2<-I(complete_data$PupilDilation)^2
complete_data$pupil2<-scale(sqrt(complete_data$pupil2))

anova(glmer(FirstLookFace ~ isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ 1 + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ sex+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ onscreen_idmean+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ missinggaze_idmean+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_missinggaze+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ qualitative+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ quantitative+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ Order*quantitative+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ Order+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ Order*qualitative+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ quantitative+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ Order*quantitative+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ Order+idwc_onscreen+age_std*quantitative+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std*quantitative+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std*qualitative+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std*Order+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std+Ordr*isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std+Order*isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std+qualitative*isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ idwc_onscreen+age_std+quantitative*isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ quantitative*Order+idwc_onscreen+age_std+isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))

model_face<-glmer(FirstLookFace ~ quantitative * Order + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit"))
model_face_pupil<-glmer(FirstLookFace ~ pupil_scaled + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit"))

anova(glmer(FirstLookFace ~ pupil_scaled + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ pupil_scaled*quantitative + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ pupil_scaled + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ pupil_scaled*qualitative + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))
anova(glmer(FirstLookFace ~ pupil_scaled + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")),glmer(FirstLookFace ~ pupil_scaled*Order + idwc_onscreen + age_std + isi_std + (1 |id),data = complete_data, family = binomial(link = "logit")))


save.image("~/pupil_sound/data_analysis.RData")


###################################################################################################################################################################################################
#### PLOTS ############################################################################################################################################################################
###################################################################################################################################################################################################

library(sdamr)
complete_data$soundtype<-ordered(complete_data$qualitative,levels=c('0','2','1'),labels=('silent','non-social','social'))

# FIG 2: sound type ON pupil dilation BY age
png('fig2_new.png',res=300,width=3000,height=2500)
rain_height <- .1
ggplot(complete_data, aes(x = "", y = PupilDilation, fill = soundtype)) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,position = position_nudge(x = rain_height+.05)) +
  # rain
  geom_point(aes(colour = soundtype), size = 2, alpha = .5, show.legend = FALSE,position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE,outlier.shape = NA,position = position_nudge(x = -rain_height*2)) +
  # mean and SE point in the cloud
  stat_summary(fun.data = mean_cl_normal, mapping = aes(color = soundtype), show.legend = FALSE,position = position_jitternudge(nudge.x = rain_height * 3)) +
  # adjust layout
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Pupil dilation") +
  coord_flip() +
  facet_wrap(~factor(age_std,levels = c(0,1),labels = c("5 months", "10 months")),nrow = 2) +
  # custom colours and theme
  scale_fill_brewer(palette = "Dark2", name = "Sound Type",labels = c("silent","non-social", "social")) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(text=element_text(size=26),panel.grid.major.y = element_blank(),legend.position = c(0.9,0.5),legend.background = element_rect(fill = "white", color = "white"))
dev.off()

# FIG 3: volume level ON pupil dilation 
png('fig3_new.png',res=300,width=3500,height=2500)
rain_height <- .1
ggplot(complete_data, aes(x = "", y = PupilDilation, fill = quantitative)) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,position = position_nudge(x = rain_height+.05)) +
  # rain
  geom_point(aes(colour = quantitative), size = 2, alpha = .5, show.legend = FALSE,position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE,outlier.shape = NA,position = position_nudge(x = -rain_height*2)) +
  # mean and SE point in the cloud
  stat_summary(fun.data = mean_cl_normal, mapping = aes(color = quantitative), show.legend = FALSE,position = position_nudge(x = rain_height * 3)) +
  # adjust layout
  scale_x_discrete(name = "Volume level") +
  scale_y_continuous(name = "Pupil dilation") +
  coord_flip() +
  facet_wrap(~factor(quantitative,levels = c(0,0.7,1),labels = c("silent","low volume", "high volume")),nrow = 3) +
  # custom colours and theme
  scale_fill_brewer(palette = "Dark2", name = "Sound volume",labels = c("silent","low volume", "high volume")) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(text=element_text(size=26),panel.grid.major.y = element_blank(),legend.position = "none")
dev.off()

# FIG 4: age ON first look latency
png('fig4_bis_new.png',res=300,width=3500,height=2500)
rain_height <- .1
ggplot(complete_data, aes(x = "", y = MinLatency, fill = age_std)) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,position = position_nudge(x = rain_height+.05)) +
  # rain
  geom_point(aes(colour = age_std), size = 2, alpha = .5, show.legend = FALSE,position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE,outlier.shape = NA,position = position_nudge(x = -rain_height*2)) +
  # mean and SE point in the cloud
  stat_summary(fun.data = mean_cl_normal, mapping = aes(color = age_std), show.legend = FALSE,position = position_nudge(x = rain_height * 3)) +
  # adjust layout
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Latency of first look to any target (ms)") +
  coord_flip() +
  facet_wrap(~factor(age_std,levels = c(0,1),labels = c("5 months", "10 months")),nrow = 2) +
  # custom colours and theme
  scale_fill_brewer(palette = "Dark2", name = "Age",labels = c("5 months", "10 months")) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(text=element_text(size=26),panel.grid.major.y = element_blank(),legend.position = "none")
dev.off()

# FIG 5: volume level ON face selection BY repetition
png('fig5.png',res=300,width=3500,height=3500)
ggplot(faceplot,aes(x=volume,y=face,fill=repetition))+geom_boxplot()+facet_wrap(~repetition)+theme_minimal()+scale_fill_brewer(palette = "Set2")+labs(x='Sound volume',y='Proportion First Look to Face')+theme(text=element_text(size=26),panel.grid.major.y = element_blank(),legend.position = "none")
dev.off()

# FIG 6: pupil & face selection
png('fig6.png',res=300,width=3500,height=3500)
ggplot(complete_data,aes(x=pupil_dilation,y=face))+geom_density()+theme_minimal()+scale_fill_brewer(palette = "Set2")+labs(x='Pupil dilation',y='Probability')+theme(text=element_text(size=26),panel.grid.major.y = element_blank())
dev.off()

# FIG 7: pupil dilation ON latency: (A) all; (B) by volume; (C) by sound type; (D) by repetition
png('fig7.png',res=400,width=4000,height=4000)
ggplot(complete_data, aes(x=pupil_scaled, y=latency_std)) +geom_point(aes(color=id),shape = 16, size = 3, show.legend = FALSE, alpha = .4) +
  stat_smooth(color='black',size=2,method = "lm", formula = y~x+I(x^2))+
  scale_x_continuous(name='Pupil dilation (scaled)')+
  scale_y_continuous(name='First look latency (scaled)')+
  theme_minimal()+theme(text=element_text(size=26),panel.grid.major.y = element_blank(),legend.position = "none")
dev.off()
