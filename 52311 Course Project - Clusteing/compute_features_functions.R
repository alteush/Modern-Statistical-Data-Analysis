#features

lastinrow<-function(s){
  #computes the last place in the row that isn't NA
  for (i in 1:c){
    if(is.na(s[i])==TRUE){
      return(i-1)
    }
  }
  return(c) 
}

getSignChanges <- function(s){
  #sign change defined to be when one element is positive and the next is negative
  # or vice versa
  #function gets a sequence and counts it's sign changes
  len <- length(s)
  sgnchange <- 0
  if (len>1){
    for (i in 1:(len-1)) {
      if(!(sign(s[i]) == sign(s[i+1])))
      {sgnchange <- sgnchange+1}
    }
  }
  return(sgnchange)
}

getModulo_2 <- function(s){
  return (s %% 2)
}

getModulo_3 <- function(s){
  return (s %% 3)
}

getModulo_5 <- function(s){
  return (s %% 5)
}

getPercentageSquared <- function(s){
  #returns the percentage of squared numbers of sequence's elements
  sqrt_qntt = 0
  len <- length(s)
  if (len == 0) {return (0)}
  for (i in 1:len){
    if (s[i]>=0){
      root <- sqrt(s[i])
      if (root == round(sqrt(root))){
        sqrt_qntt = sqrt_qntt+1
      }
    }
  }
  percent <- sqrt_qntt*100/len
  return (percent)
}

getPercentageCubed <- function(s){
  #returns the percentage of cubed numbers of sequence's elements
  cb_qntt = 0
  len <- length(s)
  if (len == 0) {return (0)}
  for (i in 1:len){
    root <- (abs(s[i]))^(1/3)
    if (root == round((abs(root))^(1/3))){
      cb_qntt = cb_qntt+1}
  }
  percent <- cb_qntt*100/len
  return (percent)
}

getPercentagePrimes <- function(s){
  #computes the percentage od prime numbers of sequence's elements
  #among the elements smaller than 2^53-1
  computable <- which(s < 2^23-1 & s > 0)
  if (length(computable) > 0){
    primes_qntt <- sum(as.double(isPrime(s[computable])))
    len_computable <- length(computable)
    return(primes_qntt/len_computable*100)}
  return (0)
}

getDifferences <- function(s){
  #given a sequence, calculates the diferences and returning the mean of them
  len <- length(s)
  differences <- s[1]
  if (len > 1){
    differences <- rep(0, len-1)
    for (i in 1:(len-1)){
      differences[i] <- s[i+1]-s[i]
    }
  }
  return(differences)
}

even_odd <- function(s){
  #computes number of even elements in a sequence
  return(sum(as.double(s %% 2 == 0)))
}

getSequence <- function(i){
  #some sequences contains more than 110 elements, therefore spread on two rows
  #that function identifies them, since all sequences start with "A".
  if (substr(as.character(stripped[i,1]),1,1) == "A" & is.na(stripped[i,2]) != TRUE){
    #read.csv reads sequences with more than 110 elements in two rows.
    #the 2-nd row doesn't start with A, therfore, shouldn't be counted.
    #there are empty sequences that shouldn't be counted
    seq <- as.double(stripped[i,])
    seq_len <- lastinrow(seq)
    seq <- seq[2:(seq_len)]
    if (!(substr(as.character(stripped[i+1,1]),1,1) == "A")){
      #in case the series expands over two rows
      extra_length <- lastinrow(stripped[i+1,])
      for (j in 1:extra_length){
        extra_int <- type.convert(as.character(stripped[i+1,j]))
        seq <- c(seq,extra_int)}
    }
    return (seq)
  }
}



adjust_features <- function(s, index){
  #function returning a 20-entry vector of sequence's features.
  #we chose to return .Machine$double.xmax value insted of Inf, -Inf,
  #in order we could preform the kmeans algorithm.
  #that's the maximal value a double can recieve.
  feat <- rep(0, 20)
  feat[1] <- index
  if (is.finite(mean(s))==TRUE) {feat[2] <- mean(s)}
  else {feat[3] <- .Machine$double.xmax}
  if (is.finite(var(s))==TRUE) {feat[3] <- var(s)}
  else {feat[3] <- .Machine$double.xmax}
  if (is.finite(median(s))==TRUE) {feat[4] <- median(s)}
  else {feat[4] <- .Machine$double.xmax}
  feat[5] <- length(s)
  if (is.finite(skewness(s))==TRUE) {feat[6] <- skewness(s)}
  if (is.finite(kurtosis(s))==TRUE) {feat[7] <- kurtosis(s)}
  #The default value for skewness and kurtosis is 0
  if (is.finite(log(abs((mean(s)))))==TRUE) {feat[8] <- log(abs(mean(s)))}
  #There are sequences in which the mean is zero or close to zero, whereas
  #the largest mean is about e^+283
  if (is.finite(log(var(s)))==TRUE) {feat[9] <- log(var(s))}
  else {feat[9] <- log(.Machine$double.xmax)}
  feat[10] <- getSignChanges(s)
  if (is.finite(mean(getDifferences(s)))==TRUE) {feat[11] <- mean(getDifferences(s))}
  else {feat[11] <- .Machine$double.xmax*sign(mean(getDifferences(s)))}
  if (is.finite(median(getDifferences(s)))==TRUE) {feat[12] <- median(getDifferences(s))}
  else {feat[12] <- .Machine$double.xmax*sign(median(getDifferences(s)))}
  feat[13] <- even_odd(s) #Even numbers
  feat[14] <- getPercentagePrimes(s) #primes percentage
  feat[15] <- getPercentageSquared(s) #squares percentage
  feat[16] <- getPercentageCubed(s) #cubes percentage
  if (max(s)!=0) {feat[17] <- mean(abs(s))/max(abs(s))}
  feat[18] <- mean(getModulo_2(s))
  feat[19] <- mean(getModulo_3(s))
  feat[20] <- mean(getModulo_5(s))
  return (feat)
}