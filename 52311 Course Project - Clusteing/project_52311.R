library(moments)
library(numbers)

source('compute_features_functions.R')

#scanning the dataset
stripped <- read.csv("C:/Users/a/Desktop/modernstats/data/stripped", header=FALSE, comment.char="#")


r = nrow(stripped)
c = ncol(stripped)
cnames <- c("Index", "Mean", "Variance", "Median", "Series' Length",
            "skewness", "Kurtosis", "Log|Mean|", "Log(Variance)", 
            "Sign Changes", "Differences' Mean", "Differences' Median",
            "Even", "Primes' Percentage", "Sqares' Percentage", "Cubes' Percentage",
            "Mean/Max", "Mean (Mod 2)", "Mean (Mod 3)", "Mean )Mod 5)")



Sys.time()
colnames(features)<-cnames
features <- matrix(data = rep(0, r*20), nrow=r, ncol = 20)
k=267633
l=266123
#k runs over stripped
#l runs over features' rows
while(k <= r){
  sequence <- getSequence(k)
  if (is.null(sequence)!=TRUE){
    features[l,] <- adjust_features(sequence, l)
    l = l+1
  }
  k = k+1
}
features <- features[1:(l-1),]
Sys.time()
results <- kmeans(features, 100, 100)
results2 <- kmeans(features[,-c(1,2,3,4,5,11,12), 100 , 100])