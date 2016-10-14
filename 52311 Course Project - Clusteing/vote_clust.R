library(proxy)
library(fastcluster)
library(stringr)

dist_means <- function(v1,v2){
  return(abs(mean(v1)-mean(v2)))
}

euc.dist <- function(v1, v2) sqrt(sum((v1 - v2) ^ 2))

readAsList <- function(path,skip=4,n=-5,delim=",",name=T){
  
  ## example usage: 
  ## path <- "OEIS/data/stripped"
  ## x <- readAsList(path,n=100) # will provide first 100 elements
  ## x <- readAsList(path,n=-5) # will provide all elements
  
  #open connection to file
  pathToFile <- path.expand(path)
  f <- file(pathToFile, "rb")  
  
  # read in raw text
  if (n<0){n <- (-10)}
  rawText <- readLines(f,n=(n+skip))
  close(f)
  rawText <- rawText[(skip+1):length(rawText)]
  
  # convert to list
  data <- sapply(rawText, str_split,delim)
  if (name){
    names(data) <- sapply(data,function(x) x[1])
    data <- sapply(data,function(x) as.numeric(x[2:(length(x)-1)]))
  }
  return(data)
}

init_clusters <- function(x_sample,fun = dist_means,k=10,method='complete'){
  print('calculating distance matrix')
  if (length(unique(sapply(x_sample,length)))==1){
    x_sample <- matrix(unlist(x_sample), ncol = length(x_sample[[1]]), byrow = TRUE)
  }
  d <- proxy::dist(x_sample, method = fun)
  print('cutting tree to receive desired number of clusters')
  clusters <- cutree(hclust(d,method),k)
  return(clusters)
}

get_score <- function(dist_vector,method='complete'){
  methods <- c('complete','single','average')
  if (!(method %in% methods)){ stop(paste('method',method,'not valid'))
  } else if (method=='complete'){ return(max(dist_vector))
  } else if (method=='single') { return(min(dist_vector))
  } else if (method=='average') { return(mean(dist_vector))
  }                                       
}
voting <- function(el,init.clusters,x_sample,v,fun=dist_means,method='complete'){
  k <- length(unique(init.clusters))
  vote <- numeric(k)
  for (i in 1:k){
    k.idx <- which(init.clusters==i)
    k.idx.n <- length(k.idx)
    if (k.idx.n > v){
      k.idx <- k.idx[sample(1:k.idx.n,v)]
    }
    dist.k <- sapply(x_sample[k.idx],function(x) fun(el,x))
    vote[i] <- get_score(dist.k,method)
  }
  cluster <- which(vote==min(vote))
  if (length(cluster)>1) { cluster <- 0}
  return(cluster)
}

vote_clust <- function(x,k=10,s=100,v=10,fun = dist_means,method='complete'){
  ## This is a heuristic approach to hierarchical clustering.
  ## In the first step, s random elements of x are clustered based on the 
  ## distance function fun, and the clustering method method, into k clusters.
  ## In the second step, the remaining n-s elements are attributed to the
  ## original clusters using a voting method. For each element, up to v
  ## representatives are randomly selected for each cluster. The element
  ## is then compared to these representatives using the clustering 
  ## method provided (can be complete, single or average). 
  
  n <- length(x)
  init.idx <- sample(1:n,s)
  other.idx <- setdiff(1:n,init.idx)
  print('Finding initial clusters:')
  init.clusters <- init_clusters(x[init.idx],fun = fun,k=k)
  print('Done finding initial clusters')
  print('Assigning remaining sequences to clusters')
  other.clusters <- sapply(x[other.idx],voting,init.clusters,x[init.idx],v,fun,method)

  final.clusters <- numeric(length=n)
  final.clusters[init.idx] <- init.clusters
  final.clusters[other.idx] <- other.clusters
  names(final.clusters) <- names(x)
  return(final.clusters)
}