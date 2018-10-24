#load libraries
library(lmerTest)
library(mosaic)
library(lme4)
library(MASS)
library(ggplot)
library(rgeos)
library(ggplot2)
library(devtools)
library(rcompanion)
library(coefplot)
library(VCA)
library(MASS)
library(lme4)
library(mosaic)

#load data
load("cov.RData")

#initialize spatiotemporal model
clim<-c("range","avg","max","min","rain","hum", "avg.sq", "max.sq", "min.sq", "pa", "b")
st.vars<-lapply(cov[clim], function(x) unlist(x))
st.vars$pa<-unlist(cov$pa)
st.vars$inc<-unlist(cov$b)
st.vars$country<-unlist(lapply(1:length(cov$names), function(i) rep(cov$country[[i]],length(cov$b[[i]]))))
st.vars$province<-unlist(lapply(1:length(cov$names), function(i) rep(cov$names[[i]],length(cov$b[[i]]))))

#function to perform backward selection on an initial formula

#inputs: 
#temp-"max", "min, or "avg" for temperature variable; 
#resp-"pa" or "b" for response variable (presence/absence or incidence respectively)
#ran-TRUE or FALSE, whether to include range in model selection

#outputs:
#aic-the Akaike information criterion for the final model
#steps-traces back explanatory variables in the order they were dropped from the model
#data-dataframe used to generate model (normalized, only includes explanatory variables selected for final model)
#org-original dataframe (not normalized, includes all explanatory variables)
#call-vector of explanatory variables used in the model selected via backward selection
#coeffs-the coefficients corresponding to the selected model
#fe-the nested province in country fixed effects

backwardselect<-function(temp, resp, ran){
  sq<-paste(temp,".sq",sep="")
  vars<-c("rain","hum",temp,sq)
  if(ran==TRUE){
    vars<-c(vars, "range")
  }
  var.remove<-c()
  aic<-c()
  
  
  #initialize dataframe, normalize variables
  if(resp=="pa"){
    fam="binomial"
    data<-as.data.frame(st.vars[c("pa","b","range", "rain", "hum","max","max.sq","min","min.sq", "avg", "avg.sq","country","province")])
    data<-data[-which(apply(data[,1:11], 1, function(x) sum(!is.finite(x)))>0),]
    org<-data
    data[,2:11]<-lapply(data[,2:11], zscore)
  } else{
    data<-as.data.frame(st.vars[c("pa","b","range", "rain", "hum","max","max.sq","min","min.sq", "avg", "avg.sq","country","province")])
    data<-data[-which(apply(data[,1:10], 1, function(x) sum(!is.finite(x)))>0),]
    data<-data[which(data$pa==1),]
    org<-data
    data[,2:11]<-lapply(data[,2:11], zscore)
    fam="gaussian"
  }
  
  #fit a mixed effect model with nested effect of province in country with remaining response variables
  for(i in 1:length(vars)){
    if(i==1){
      use.vars<-vars
    } else {
      use.vars<-vars[-which(vars %in% var.remove)]
    }
    use.vars<-c(use.vars, "country","province")
    df<-data[use.vars]
    df$resp<-data[[resp]]
    
    if(resp=="pa"){
      mod<-glmer(resp~.+(1|country/province)-1, data=df, fam=fam, nAGQ=0)
      k=4
    } else {
      mod<-lmer(resp~.+(1|country/province)-1,data=df)
      k=5
    }
    
    #if any explanatory variables have p>.05, drop the explanatory variable with the greatest p-value (backward selection)
    if(sum(summary(mod)$coefficients[use.vars[-(-1:0+length(use.vars))],k]>.05)>0){
      remove.me<-rownames(summary(mod)$coefficients)[which.max(summary(mod)$coefficients[use.vars[-(-1:0+length(use.vars))],k])]
      var.remove<-c(var.remove, remove.me)
      aic<-c(aic,extractAIC(mod)[2])
      
    } else {
      call<-use.vars
      coeff<-summary(mod)$coefficients[use.vars[-(-1:0+length(use.vars))]]
      return(list("aic"=aic, "data"=data, "org"=org, "steps"=var.remove, "call"=call, "coeffs"=summary(mod)$coefficients[use.vars[-(-1:0+length(use.vars))],], "fe"=tail(summary(mod)$coefficients,-length(use.vars)+2), "residuals"=residuals(mod)))
    }
    
  }
  call<-use.vars
  return(list("aic"=aic, "steps"=var.remove, "data"=data, "org"=org, "call"=call, "coeffs"=summary(mod)$coefficients[use.vars[-(-1:0+length(use.vars))],], "fe"=tail(summary(mod)$coefficients,-length(use.vars)+2)), "residuals"=residuals(mod))
}


#use backward selection to generate the best model for each combination of temperature factors (i.e.: min, max, avg, with or without range)
max.r.pa<-backwardselect("max","pa",TRUE)
avg.r.pa<-backwardselect("avg","pa",TRUE)
min.r.pa<-backwardselect("min","pa",TRUE)
max.pa<-backwardselect("max","pa",FALSE)
avg.pa<-backwardselect("avg","pa",FALSE)
min.pa<-backwardselect("min","pa",FALSE)

max.r.inc<-backwardselect("max","b",TRUE)
avg.r.inc<-backwardselect("avg","b",TRUE)
min.r.inc<-backwardselect("min","b",TRUE)
max.inc<-backwardselect("max","b",FALSE)
avg.inc<-backwardselect("avg","b",FALSE)
min.inc<-backwardselect("min","b",FALSE)

the.aic<-matrix(ncol=6, nrow=2)
rownames(the.aic)<-c("pa","b")
colnames(the.aic)<-c("min.r","max.r","avg.r","min", "max","avg")

#matrix of AIC values corresponding to each model. Presence/absence models on top row,
#incidence models on bottom. Temperature variables used to initialize models are 
#indicated across the top (.r suffix indicates range was included)
the.aic[1,1]<-tail(min.r.pa$aic,1)
the.aic[1,2]<-tail(max.r.pa$aic,1)
the.aic[1,3]<-tail(avg.r.pa$aic,1)
the.aic[1,4]<-tail(min.pa$aic,1)
the.aic[1,5]<-tail(max.pa$aic,1)
the.aic[1,6]<-tail(avg.pa$aic,1)
the.aic[2,1]<-tail(min.r.inc$aic,1)
the.aic[2,2]<-tail(max.r.inc$aic,1)
the.aic[2,3]<-tail(avg.r.inc$aic,1)
the.aic[2,4]<-tail(min.inc$aic,1)
the.aic[2,5]<-tail(max.inc$aic,1)
the.aic[2,6]<-tail(avg.inc$aic,1)

#variance component analysis for both steps of spatiotemporal model
anovaVCA(pa~country/province, max.pa$data)
anovaVCA(b~country/province, data)

#SPATIAL MODEL
read.csv("covprov.csv")

#response variables
rv<-c("local.trans","max.b","avg.b","cum.cases","avg.cases")

#provinces that "pass the hurdle" of having any local transmission
ph<-which(cov.prov$local.trans==1)

#list of spatial models using average temperature wth a country fixed effect
#the tag for each model is its response variable, following the vector "rv"
sp.avg.fe<-list()

#initialize matrix of coefficients, with explanatory variables across the 
#columns and response variables across the rows.
cf<-matrix(nrow=5,ncol=6)
colnames(cf)<-c("city.pop","mean.rain","mean.hum","mean.range","mean.mean", "mean.mean.sq")
rownames(cf)<-rv[1:5]

for(i in 1:5){
  if(i ==1){
    fam<-"binomial"
  }
  if(i %in% 2:5){
    fam<-"gaussian"
  }
  if(i%in%c(2,3)){
    data<-lapply(cov.prov, function(x) x[ph])
    data$resp<-as.numeric(unlist(data[rv[i]]))
    data[2:13]<-lapply(data[2:13], zscore)
  }
  else{
    data<-cov.prov
    data$resp<-as.numeric(unlist(data[rv[i]]))
    if(i==1){
      data[c(2:10,13)]<-lapply(data[c(2:10,13)], zscore)
    }
    else{
      data[2:13]<-lapply(data[2:13], zscore)
    }
  }
  
  #use bidirection selection based on AIC to select the best model for a given response variable
  sp.avg.fe[[rv[i]]]<-stepAIC(glm(data[[rv[i]]]~city.pop+mean.rain+mean.hum+mean.mean+mean.range+mean.mean.sq+as.factor(country)-1,family=fam,na.action="na.omit",data=data),trace=F)

  for(c in names(sp.avg.fe[[i]]$coefficients)){
    if(c %in% colnames(cf)){
      ind<-which(colnames(cf)==c)
      cf[i,ind]<-sp.avg.fe[[i]]$coefficients[which(names(sp.avg.fe[[i]]$coefficients)==c)]
    }
  }
}

##PseudoR2 calculations
#matrix of marginal & condition R2 for spatiotemporal hurdle model (calculated using code outline by Nakagawa & Schielzeth at doi: 10.1111/j.2041-210x.2012.00261.x)
rsq<-matrix(nrow=2,ncol=2)
colnames(rsq)<-c("marg","cond")
rownames(rsq)<-c("pa","b")

mF <- lmer(inc ~ rain + avg.sq + range +(1|country/province), data = avg.r.inc$data)
VarF <- var(as.vector(fixef(mF) %*% t(model.matrix(mF))))

rsq[1,1]<-VarF/(VarF + VarCorr(mF)$'province:country'[1] + VarCorr(mF)$'country'[1] + attr(VarCorr(mF), "sc")^2)
rsq[1,2]<-(VarF + VarCorr(mF)$'province:country'[1] + VarCorr(mF)$'country'[1])/(VarF + VarCorr(mF)$'province:country'[1] + VarCorr(mF)$'country'[1] + (attr(VarCorr(mF), "sc")^2))

mF <- glmer(pa ~ hum + max + (1|country/province), family = "binomial", data = max.pa$data)
VarF <- var(as.vector(fixef(mF) %*% t(model.matrix(mF))))

rsq[2,1]<-VarF/(VarF + VarCorr(mF)$'province:country'[1] + VarCorr(mF)$'country'[1] + pi^2/3)
rsq[2,2]<-(VarF + VarCorr(mF)$'province:country'[1] + VarCorr(mF)$'country'[1])/(VarF + VarCorr(mF)$'province:country'[1] + VarCorr(mF)$'country'[1] + pi^2/3)

#Pseudo R2 for spatial models
sapply(sp.avg.fe, function(x) nagelkerke(x)$Pseudo.R.squared.for.model.vs.null[3])
