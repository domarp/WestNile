##=========================================================================================
## Script for Rdige regression
##=========================================================================================

rm(list = ls())

library(glmnet)
library(MASS)
library(caTools)

##=========================================================================================
## Set working directory
##=========================================================================================
setwd("C:\\Users\\prbalakrishnan\\Desktop\\New folder\\")
source("eda_utilities.r")

##=========================================================================================
## read Files
##=========================================================================================
score.data = readRDS("score.data.rds")
score.data$Date <- NULL
train.data = readRDS("model.data.rds")
train.data.sub = train.data[,!(names(train.data) %in% c("Id","WnvPresent","Species"))]
train.data = cbind(train.data[,c("Id","WnvPresent","Species")], train.data.sub)
train.data[is.na(train.data)] = 0

##=========================================================================================
## Function to compute risk/aggregations for categorical values
##=========================================================================================
createRisk <- function(dt,dts,list.vars){
  n_cols = ncol(dts)
  cnt = n_cols
  for (var in list.vars){
    if(!(var %in% names(dt))){
      next
    }
    cnt = cnt + 1
    dt[[var]] <- trimws(dt[[var]])
    dts[[var]] <- trimws(dts[[var]])
    cf = paste("rsk_mean_",var,sep=".")
    lookup <- as.data.frame(dt %>% group_by_(var) %>% 
                              summarize_(rsk=paste0("log(mean(as.numeric(WnvPresent),na.rm=TRUE))")))
    names(lookup)[2] <- paste0("rsk_mean_",var)
    dt <- merge(dt,lookup,by=var)
    dts <- merge(dts, lookup,by=var,all.x=TRUE)
    dts[which(is.na(dts[cnt])==T),cnt] = 0
  }
  dt <- dt[,-which(names(dt) %in% list.vars)]
  dts <- dts[,-which(names(dts) %in% list.vars)]
  
  rm(lookup)
  
  return(list(dt,dts))
}

trn.r2 <- c()
val.r2 <- c()
k = 2
imp.vars <- c()

##=========================================================================================
## Run Model : Ridge
##=========================================================================================

for(i in 1:1){
  set.seed(6)
  
  var.subsample = sample(4:ncol(train.data),ncol(train.data)-3)
  train.data.1  = train.data[,c(1,2,3,var.subsample)]
  
  train.data.1$cross.rand <- runif(nrow(train.data.1))
  for(cross in 1:k){
    if(cross == 2){
      next
    }
    train.val.indx = which(train.data.1$cross.rand < (cross * (1/k)) & 
                             train.data.1$cross.rand > ((cross-1) * (1/k)))
    train.trn.indx = which(!(1:nrow(train.data.1) %in% train.val.indx))
    train.trn = train.data.1[train.trn.indx,2:ncol(train.data.1)]
    train.val = train.data.1[train.val.indx,2:ncol(train.data.1)]
    train.trn$cross.rand <- NULL
    train.val$cross.rand <- NULL
    train.trn$Date <- NULL
    train.val$Date <- NULL
    
    rsk = createRisk(train.trn, train.val,c("orig.Trap", "Trap","Block","Species","Zipcode"))
    train.trn = rsk[[1]]
    train.val = rsk[[2]]
    
    train.trn.act = train.trn[train.trn$rsk_mean_Species != 0.000,]
    train.val.act = train.val[train.val$rsk_mean_Species != 0.000,]
    
    cv = cv.glmnet(as.matrix(train.trn.act[,2:ncol(train.trn)]), as.matrix(train.trn.act$WnvPresent),
                family = "binomial",alpha = 0)
    fit = glmnet(as.matrix(train.trn.act[,2:ncol(train.trn)]), as.matrix(train.trn.act$WnvPresent),
                 family = "binomial",alpha = 0, lambda = cv$lambda.1se)
    y.pred.trn = predict(fit,as.matrix(train.trn[,2:ncol(train.trn)]),type="response",
                         lambda = cv$lambda.1se)
    y.pred.trn[train.trn$rsk_mean_Species == 0.00,] = 0.0
    t.r2 <- colAUC(as.matrix(y.pred.trn), train.trn$WnvPresent, alg = "ROC")
    cat("\n\nTRAIN",i,": ")
    cat(t.r2, " Rate: ", mean(as.numeric(train.trn$WnvPresent)-1, na.rm = TRUE ))
    cat("\nVAL",i,": ")
    y.pred.val = predict(fit,as.matrix(train.val[,2:ncol(train.val)]),type="response",
                         lambda = cv$lambda.1se)
    y.pred.val[train.val$rsk_mean_Species == 0.00,] = 0.0
    r2 <- colAUC(as.matrix(y.pred.val), train.val$WnvPresent, alg = "ROC")
    cat(r2, " Rate: ", mean(as.numeric(train.val$WnvPresent)-1, na.rm = TRUE))
    
    mi <- coef(fit)[2:length(coef(fit))]
    ms <- sapply(train.trn[,2:ncol(train.trn)], sd)
    
    imp.vars <- cbind(imp.vars, mi * ms)
  }
}

##=========================================================================================
## Write submission output
##=========================================================================================
rsk = createRisk(train.data, score.data,c("orig.Trap", "Trap","Block","Species","Zipcode"))
score.data = rsk[[2]]

options(scipen=999)
submit.predict = predict(fit, as.matrix(score.data[,names(train.val[,2:ncol(train.val)])]), 
                         type = "response")
submission = data.frame(cbind(score.data$Id, submit.predict))
names(submission) = c("Id","WnvPresent")
# submission$WnvPresent[score.data$rsk_mean_Species == 0.00000] = 0.0
submission[is.na(submission)] = 0
write.csv(submission,"ridge_baseline.csv",row.names = FALSE)