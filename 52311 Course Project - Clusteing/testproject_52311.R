getUrl <- function(index){
  #function gets a sequence's index. Extracting it's OEIS' crossrefs
  character_index <- as.character(index)
  index_len <- nchar(character_index)
  to_paste <- c("https://oeis.org/","A",rep(0,(6-index_len)),character_index)
  page_url <- paste(to_paste, collapse = "")
  return(page_url)
}

getLinks <- function(page_url){
  #function returning all links to sequences (all links that are "A" followed by
  #a number)
  page <- read.csv(url(page_url), sep = " ")
  page <- as.matrix(page)
  links <- vector()
  for (i in 1:(length(page))){
    if (is.na(as.double(substr(strsplit(page[i],"A")[[1]][2],1,1)))!=TRUE){
      ref <- as.double(substring(strsplit(page[i], "A")[[1]],1,6))
      links <- cbind(links, ref)}
  }
  links <- links[which(is.na(links) != TRUE)]
  return (unique(links))
}

getCrossrefs <- function(page_url, index){
  #returning the sequences appearing in the line "sequence in context"
  page <- read.csv(url(page_url), sep = "\n")
  page <- as.matrix(page)
  #a loop finding which line in the page contain the words "Sequence in context"
  for (i in 1:(length(page))){
    if (length(grep("Sequence in context", page[i])) > 0){
      crossrefs <- unique(as.double(substr(strsplit(page[i],"A")
                                           [[1]],1,6)))
      crossrefs <- crossrefs[which(is.finite(crossrefs) & crossrefs != index)]
      return(crossrefs)
    }
  }
  return(0)
}



test <- function(model){
  #function recieving a clustering model, returning a sample of 100 sequences'
  #and their matching percentages
  l <- length(model)
  sample <- sample(1:(l-1), 1000)
  significance <- rbind(sample, matrix(data = 0, ncol = 1000, nrow = 5))
  for (i in 1:1000){
    s <- sample[i]
    page <- getUrl(s)
    links <- getLinks(page)
    sample_cluster <- as.double(model[s])
    cluster_size <- length(which(model == sample_cluster))
    links_clusters <- as.double(model[links])
    crossrefs <- getCrossrefs(page, s)
    crossrefs_clusters <- as.double(model[crossrefs])
    percentage1 <- length(which(links_clusters == sample_cluster))/
      length(links_clusters)*100
    percentage2 <- length(which(links_clusters == sample_cluster))/
      cluster_size*100
    matching_crossrefs <- length(which(crossrefs_clusters == model[s]))
    num_crossrefs <- length(crossrefs)
    percentage3 <- matching_crossrefs/num_crossrefs*100
    significance[2,i] <- percentage1
    significance[3,i] <- percentage2
    significance[4,i] <- percentage3
    significance[5,i] <- length(links)
    significance[6,i] <- length(crossrefs)
  }
  return(significance)
}