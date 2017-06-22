##=========================================================================================
## Script for Random Forest Model
##=========================================================================================

rm(list = ls())

library(randomForest)
library(MASS)
library(caTools)

##=========================================================================================
## Set working directory
##=========================================================================================
setwd("C:\\Users\\prbalakrishnan\\Desktop\\New folder\\")
source("eda_utilities.r")

##=========================================================================================
## read Files & arrange data
##=========================================================================================
score.data = readRDS("score.data.rds")
train.data.1 = readRDS("model.data.rds")
train.data.sub = train.data.1[,!(names(train.data.1) %in% c("Id","WnvPresent"))]
train.data.1 = cbind(train.data.1[,c("Id","WnvPresent")], train.data.sub)

##=========================================================================================
## Imputation
##=========================================================================================
# Zipcode : Impute with TrapID
# Issue to fix : impute with actual zipcode values
train.data.1$Zipcode[which(is.na(train.data.1$Zipcode))] <- train.data.1$Trap[which(is.na(train.data.1$Zipcode))]
score.data$Zipcode[which(is.na(score.data$Zipcode))] <- score.data$Trap[which(is.na(score.data$Zipcode))]
score.data$Zipcode <- as.factor(score.data$Zipcode)
train.data.1$Zipcode <- as.factor(train.data.1$Zipcode)

count.missing <- sapply(rbind(train.data.1,score.data), is.na)
count.missing <- colSums(count.missing)

# For some weather variables ~1% data is missing, fix them by mean imputation
for(cols in which(count.missing > 0 & count.missing < nrow(train.data.1))){
  train.data.1[which(is.na(train.data.1[[cols]])),cols] <- 
    mean(train.data.1[which(!is.na(train.data.1[[cols]])),cols])
}

##=========================================================================================
## Function to compute risk/aggregations for categorical values
##=========================================================================================

createRisk <- function(dt,dts,list.vars, target = "WnvPresent"){
  n.cols = ncol(dts)
  cnt = n.cols
  for (var in list.vars){
    if(!(var %in% names(dt))){
      next
    }
    cnt = cnt + 1
    dt[[var]] <- trimws(dt[[var]])
    dts[[var]] <- trimws(dts[[var]])
    lookup <- as.data.frame(dt %>% group_by_(var) %>% 
                              summarize_(rsk=paste0("mean(as.numeric(", target,"),na.rm=TRUE)")))
    names(lookup)[2] <- paste0("avg_",var)
    dt <- merge(dt,lookup,by=var, sort=TRUE)
    dts <- merge(dts, lookup,by=var,all.x=TRUE, sort=TRUE)
    dts[which(is.na(dts[cnt])==T),cnt] = 0
  }
  dt <- dt[,-which(names(dt) %in% list.vars)]
  dts <- dts[,-which(names(dts) %in% list.vars)]
  
  rm(lookup)
  
  return(list(dt,dts))
}

##=========================================================================================
## Quick check on predictability of numeric variables
##=========================================================================================
factor.vars <- which(sapply(train.data.1, class) %in% c("factor"))
all.var.auc <- colAUC(as.matrix(train.data.1[,-factor.vars]), train.data.1$WnvPresent, alg = "ROC")

##=========================================================================================
## Quick check on correlation amongst numeric variables
##=========================================================================================
corr <- cor(train.data.1[,-c(1,factor.vars)])
write.csv(corr,"corr2.csv",row.names = FALSE)

##=========================================================================================
## Run Model : RF CV
##=========================================================================================

set.seed(1000)

k = 4 #3 foldCV
train.data.1$cross.rand <- runif(nrow(train.data.1))
drop.cols <- c("Id","WetBulb","cross.rand","Date","Trap","Tmin","WoY","isSat",
               "orig.Trap","Month")

for(cross in 1:k){
  train.val.indx = which(train.data.1$cross.rand < (cross * (1/k)) & 
                           train.data.1$cross.rand > ((cross-1) * (1/k)))
  train.trn.indx = which(!(1:nrow(train.data.1) %in% train.val.indx))
  train.trn = train.data.1[train.trn.indx,2:ncol(train.data.1)]
  train.val = train.data.1[train.val.indx,2:ncol(train.data.1)]
  train.trn <- train.trn[,-which(names(train.trn) %in% drop.cols)]
  train.val <- train.val[,-which(names(train.val) %in% drop.cols)]
  
  rsk = createRisk(train.trn, train.val,c("orig.Trap","Block","Zipcode"))
  train.trn = rsk[[1]]
  train.val = rsk[[2]]
  
  train.trn.act = train.trn[train.trn$Species != "MISC",]
  train.val.act = train.val[train.val$Species != "MISC",]
  
  fit = randomForest(WnvPresent ~., data=train.trn, ntree = 1000, 
                     mtry = 20, maxnodes = 8, classwt = c(0.05,0.95))
  y.pred.trn = predict(fit,train.trn, type = "prob")
  y.pred.trn[train.trn$Species == "MISC",] = 0.0
  t.auc <- colAUC(as.matrix(y.pred.trn), train.trn$WnvPresent, alg = "ROC")[1]
  
  y.pred.val = predict(fit,train.val, type = "prob")
  y.pred.val[train.val$Species == "MISC",] = 0.0
  v.auc <- colAUC(as.matrix(y.pred.val), train.val$WnvPresent, alg = "ROC")[1]
  
  cat("\n\nTRAIN",cross,": ")
  cat(t.auc, " Rate: ", mean(as.numeric(train.trn$WnvPresent)-1, na.rm = TRUE ))
  cat("\nVAL",cross,": ")
  cat(v.auc, " Rate: ", mean(as.numeric(train.val$WnvPresent)-1, na.rm = TRUE))
  
  varImpPlot(fit)
}

##=========================================================================================
## Write submission output
##=========================================================================================
train.data.1 = train.data.1[,-which(names(train.data.1) %in% drop.cols)]
rsk = createRisk(train.data.1, score.data,c("Block","Zipcode","orig.Trap"))
train.data.1 = rsk[[1]]
score.data = rsk[[2]]

fit = randomForest(WnvPresent ~., data=train.data.1, ntree = 1000, 
                   mtry = 20, maxnodes = 8, classwt = c(0.05,0.95))

train.score = predict(fit, train.data.1, type = "prob")
train.score[train.data.1$Species == "MISC",] = 0.0
colAUC(as.matrix(train.score), train.data.1$WnvPresent, alg = "ROC")[1]

options(scipen=999)
submit.predict = predict(fit, score.data, type = "prob")
submission = data.frame(cbind(score.data$Id, submit.predict[,2]))
submission$WnvPresent[score.data$Species == "MISC"] = 0.0
names(submission) = c("Id","WnvPresent")
submission[is.na(submission)] = 0
write.csv(submission,"randomForest_baseline.csv",row.names = FALSE)





