##=========================================================================================
## Utility functions for exploratory data analysis
## The script performs univariate and bivariate analysis between different columns
## statsNum : For profiling numeric/ integer variables
## statsNonNum : For profiling non numeric variables
## plotNumHist : Histogram for numeric variables
## plotNonNumFreq : Frequency box-plot for non-numeric variables
## bivariate : Bivariate analysis based on class type
##=========================================================================================

library(dplyr)
library(ggplot2)
library(MASS)
library(lazyeval)

##=========================================================================================
## Function to compute stats for numeric columns
##=========================================================================================

statsNum <- function(iData){
  col.names <- colnames(iData)
  stats.num <- NULL
  for(col in col.names){
    dat <- iData[[col]]
    if(class(iData[[col]])[1] %in% c("integer","numeric")){
      stats.num <- rbind(stats.num,
                     data.frame(Variable = col,
                                Num.Records = length(dat),
                                dType = class(iData[[col]])[1],
                                Num.Unique = length(unique(dat)),
                                Pct.Missing=100 * sum(is.na(dat))/length(dat),
                                Mean = mean(dat, na.rm = TRUE),
                                Median = median(dat, na.rm = TRUE),
                                Min = min(dat, na.rm = TRUE),
                                perc.25 = quantile(dat,c(0.25), na.rm=TRUE),
                                perc.75 = quantile(dat,c(0.75),na.rm=TRUE),
                                Max = max(dat, na.rm = TRUE)))
    }
  }
  row.names(stats.num) <- NULL
  return(stats.num)
}

##=========================================================================================
## Function to compute stats for non-numeric columns 
##=========================================================================================

statsNonNum <- function(iData){
  col.names <- colnames(iData)
  stats.non.num <- NULL
  for(col in col.names){
    dat <- iData[[col]]
    if(class(iData[[col]])[1] %in% c("factor","logical","character")){
      freq <- table(dat)
      freq <- sort(freq, decreasing=TRUE)
      if(length(freq) == 1){
        most.freq.1 <- paste0(names(freq)[1]," : ",freq[[1]])
        most.freq.2 <- ""
        most.freq.3 <- ""
      }
      if(length(freq) == 2){
        most.freq.1 <- paste0(names(freq)[1]," : ",freq[[1]])
        most.freq.2 <- paste0(names(freq)[2]," : ",freq[[2]])
        most.freq.3 <- ""
      }
      if(length(freq) >= 3){
        most.freq.1 <- paste0(names(freq)[1]," : ",freq[[1]])
        most.freq.2 <- paste0(names(freq)[2]," : ",freq[[2]])
        most.freq.3 <- paste0(names(freq)[3]," : ",freq[[3]])
      }
      stats.non.num <- rbind(stats.non.num,
                     data.frame(Variable = col,
                                Num.Records = length(dat),
                                dType = class(iData[[col]])[1],
                                Num.Levels = length(unique(dat)),
                                Pct.Missing = 100 * sum(is.na(dat))/length(dat),
                                Top.1 = most.freq.1,
                                Top.2 = most.freq.2,
                                Top.3 = most.freq.3))
    }
  }
  row.names(stats.non.num) <- NULL
  return(stats.non.num)
}

##=========================================================================================
## Function for univariate plots for numeric variables
##=========================================================================================

plotNumHist <- function(iData, var.list){
  col.names <- colnames(iData)
  for(col in col.names){
    if(!(col %in% var.list)){
      next
    }
    dat <- iData[[col]]
    if(class(iData[[col]])[1] %in% c("integer","numeric")){
      truehist(dat, main = paste0("Histogram : ", col), xlab = col)
    }
  }
}

##=========================================================================================
## Function for univariate plots for categorical variables
##=========================================================================================

plotNonNumFreq <- function(iData, var.list, n.factors = 25){
  col.names <- colnames(iData)
  for(col in col.names){
    if(!(col %in% var.list)){
      next
    }
    dat <- iData[[col]]
    if(class(iData[[col]])[1] %in% c("factor","character","logical")){
      if(length(unique(dat)) >= n.factors){
        freq = table(dat)
        freq = sort(freq, decreasing = TRUE)
        sel.values = names(freq)[1:n.factors]
        dat = dat[which(dat %in% sel.values)]
        print("more than 10 factors; Reducing factors")
      }
      dat <- data.frame(dat)
      names(dat) <- c(col)
      plt = ggplot(dat, aes_string(x=paste0("reorder(",colnames(dat)[1],",",colnames(dat)[1],
                                            ",function(x) length(x))"))) + xlab(col) + 
        geom_bar(position = position_stack(reverse=TRUE)) + coord_flip() 
      print(plt)
    }
  }
}

##=========================================================================================
## Function for bivariate analysis : Numerical x Numerical
##=========================================================================================

# call from main bivariate function
bivariateNumNum <- function(iData, var.list, do.plot = 0){
  print("Correlation:")
  print(cor(as.numeric(iData[[var.list[1]]]), as.numeric(iData[[var.list[2]]])))
  if(do.plot){
    plt <- ggplot(iData, aes_string(var.list[1], var.list[2])) + geom_point()
    print(plt)
  }
}

##=========================================================================================
## Function for bivariate analysis : Numerical x Factor
##=========================================================================================

# call from main bivariate function
bivariateNumNon <- function(iData, var.list, do.plot = 0){
  agg.data = iData %>% group_by_(var.list[2]) %>% 
    summarize_(Mean = interp(~mean(var, na.rm = TRUE), var = as.name(var.list[1])))
  print(agg.data)
  if(do.plot){
    plt <- ggplot(agg.data, aes_string(var.list[2])) + geom_bar(aes(weight = Mean))
    plt <- plt + ylab(var.list[1])
    print(plt)
  }
}

##=========================================================================================
## Function for bivariate analysis : Factor x Factor
##=========================================================================================

# call from main bivariate function
bivariateNonNon <- function(iData, var.list, do.plot = 0){
  agg.data = iData %>% group_by_(var.list[1], var.list[2]) %>% 
    summarize(Count= n())
  print(agg.data)
  if(do.plot){
    plt <- ggplot(agg.data, aes_string(var.list[2],var.list[1])) 
    plt <- plt + geom_tile(aes(fill = Count))
    print(plt)
  }
}

##=========================================================================================
## Function for all bivariate analysis
##=========================================================================================

bivariate <- function(iData, var.list, do.plot = 0){
  class.1 <- class(iData[[var.list[1]]])
  class.2 <- class(iData[[var.list[2]]])
  print(paste0("Calling bivariate function for ", class.1, " x ", class.2))
  if(class.1 %in% c("numeric", "integer")){
    if(class.2 %in% c("numeric","integer")){
      bivariateNumNum(iData, var.list, do.plot)
    } else {
      bivariateNumNon(iData, var.list, do.plot)
    }
  } else {
    bivariateNonNon(iData, var.list, do.plot)
  }
}

# setwd("C:\\Users\\prbalakrishnan\\Desktop\\New folder\\")
# id <- read.csv("train.csv")
# compute.stats.num(id)
# compute.bivariate(id, c("Trap", "Species"),1)
