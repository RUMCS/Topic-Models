library(topicmodels)
##--------------------------------------------------------------------------
## This function reads the data from csv file, cleans it and returns tm text Corpus
## readFromfileName - name of the file to read from
##   year filter
##   year  - optional: year when docs were created, format = YYYY
##--------------------------------------------------------------------------
createCorp <- function(readFromfileName, year){
  cat("@ createCorp \n")
  library(tm)
  ## Load the data
  Posts <- read.delim(file = readFromfileName, header = T, quote = "", sep = "\t")
  
  cat("Read", nrow(Posts), "rows from", readFromfileName, "\n")
  cat("year to process =", year, "\n")
  
  ## the following lines to filter by year
  if(! missing(year)){
    Posts$create_ts <- as.POSIXlt(Posts$create_ts)
    Posts <- subset(Posts, Posts$create_ts$year == (year - 1900))
    cat("Kept", nrow(Posts), "rows\n")
  }
  
  return( doCorpCreation(Posts) )
}

## this is a private function for corpus creation
doCorpCreation <- function(Posts){
 
  Posts <- data.frame(doc_id = Posts$id, question = Posts$title, txt = paste(Posts$title, Posts$body, Posts$answers, sep =" ")) ## added question
  ## Build a corpus
  posts.reader <- readTabular(mapping=list(content="txt", id="doc_id", quest="question")) ## added question
  
  corp <- VCorpus(DataframeSource(Posts), readerControl=list(reader=posts.reader))
  
  ## Transform data
  ## mc.cores = 1 to fix run time error
  
  corp = tm_map(corp, content_transformer(tolower), mc.cores=1) #converting to lower case
  corp = tm_map(corp, removeWords, stopwords('english'), mc.cores=1) # Remove Stopwords
  corp = tm_map(corp, stemDocument, mc.cores=1) # Stemming
  corp = tm_map(corp, stripWhitespace, mc.cores=1) # Eliminate whitespace char
    
  cat ("Created corpus from", length(corp), "documents\n")
  return (corp)
}

## Removes topics from the Document Term Matrix using tf-idf approach
## dtm - Document Term Matrix to process
## threshold - remove words with tf-idf values smaller than threshold
##    this parameter is optional: if threshold is not defined, we will remove 
##    most frequent words, based on the median value of tf-idf (at least 50%)
##    Setting threshold = 0 will eliminate words appearing in every document
removeFrequentWords <- function(dtm, threshold){
  library("slam")
  termTfidf <- tapply(dtm$v/row_sums(dtm)[dtm$i], dtm$j, mean) *
    log2(nDocs(dtm)/col_sums(dtm > 0))
  
  # if no theshold is provided then set the thrshold value to median of tf-idf
  # distribution, removing at least 50% of the words
  if( missing(threshold) ) {
    threshold <- median(termTfidf)
  }
  
  # remove terms which have tf-idf smaller than the median of all the tf-idf values
  # this will shrink the dictionay by approximately 50%
  dtm <- dtm[, termTfidf > threshold]
  dtm <- dtm[row_sums(dtm) > 0,]  #remove docs that have no terms remaining (unlikely event)
  return(dtm)
}

## This function returns frequency of topics for a given LDA model
## dat - Document Term Matrix to feed to the LDA model
## topicCount - Number of topics in the LDA model
getTopicsFrequency <- function(dat, topicCount){
  mdl <- LDA(dat, topicCount) #LDA model
  
  mdl.alpha <- mdl@alpha
  mdl.beta.mean <- mean(mdl@beta)
  mdl.beta.sd <- sd(mdl@beta)
  
  return( list(
    mdl.alpha = mdl.alpha,
    mdl.beta.mean = mdl.beta.mean,
    mdl.beta.sd = mdl.beta.sd,
    topic.frequency = as.vector(table(topics( mdl )))
  ))
}

incCalc <- function(totCount) {
  ## Function to calculate increments, 2..100, 1 101..200 5, 201..300 10 ... capped to 25 after totCount >= 600 to 25
  numRanges <- ceiling(totCount/100) ## this is to ensure that if goes beyond 100's it is taken into account
  vecRanges = matrix(nrow = 1, ncol = numRanges)
  vecCount <- c()
  MRanges <- matrix(nrow = numRanges,ncol = 2)
  ix <- 1
  
  maxCount <- 1
  cat("numRanges=", numRanges, "\n")
  for(i in 1:numRanges){
    indI <- (100*i)-100+1
    indJ <- 100*i
    
    MRanges[i,1] <- indI
    MRanges[i,2] <- indJ
    
    cat("i=", i, "from=", MRanges[i,1], "to ", MRanges[i,2])
    
    cat ("indI=", indI, "indJ=", indJ, "i=", i, "\n")
    
    p <- i
    
    if (p > 1 && p < 7) {p <- 5*(i-1)} else {if (p!=1){p <-25}}
    
    cat("p=",p,"\n")
    
    for (j in seq(from=MRanges[i,1], to=MRanges[i,2], by=p)){
      
      if ( j < totCount ) { ## This is to limit the calculation of intervals up to totCount
        vecCount[ix] <- j
        maxCount <- ix
        cat("j =", j, "\n")
        cat("vecCount[", ix, "]=", vecCount[ix], "\n")
        ix <- ix + 1
      }
    }
    cat("maxCount=", maxCount, "\n")
  }
  return(vecCount)
}
