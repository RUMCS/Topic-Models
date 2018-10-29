library(topicmodels)
##------------------------------------------------------------------------------
## The function takes term document matrix (sparse representation), tdm, and exports it to file
## fileName in the following format "distinctWordCount wordId:WordCount ...", one document per 
## line. For example, "2 15:7 21:4" means that a document contains two distinct words 
## with ids 15 and 21, appearing 7 and 4 times, respectively. 
## This format is used by LDA and HDP-LDA code provided by Blei's team.
##
## The function saves dictionary of words, one word per line to file fileName, suffixed by ".dic".
## Id of a given word is given by word's line number in the file fileName.dic;
## line number count starts from 1.
##
## The function saves document "names", one name per line to file fileName, suffixed by ".doc".
## Id of a given document is equivalent to the document's line number in the file fileName.doc;
## line number count starts from 1.
##
## Usage example: export2lda_c(myTermDocMatrix, "~/foo.txt")
## Note that the existing file will be overwritten
##----------------------------------------------------------------------------
export2lda_c <- function(tdm, fileName){
  # We assume that
  # j - term (a.k.a. word) ids
  # i - doc ids
  # v - word frequency
  
  fid <- file(fileName, "w") # overwrite existing file
  docIds <- unique(tdm$i)
  for(docId in docIds){
    # I am not certain if $j is guaranteed to be consecutive; 
    # hence the inefficient search using which
    wordIndexes <- which( tdm$i == docId)
    cat(length(wordIndexes), file = fid, sep = "", append = T)
    for(wordIndex in  wordIndexes){
      cat(" ", tdm$j[wordIndex], ":", tdm$v[wordIndex], file = fid, sep = "", append = T)
    }
    cat("\n", file = fid, sep = "", append = T)
  }
  close(fid)

  # save dictionary of words
  write(tdm$dimnames$Terms, file = paste(fileName, ".dic", sep=""), sep = "\n")
  
  # save original document names
  write(tdm$dimnames$Docs, file = paste(fileName, ".doc", sep=""), sep = "\n")
}

## This function reads the data from csv file, cleans it and returns tm text Corpus
## readFromfileName - name of the file to read from
## Quarterly Filter (mandatory):
##   year  - mandatory: year when docs were created, format = YYYY
##   month - mandatory: month when docs were created
createCorpQuarter <- function(readFromfileName, year, quarter){
  library(tm)
  ## Load the data
  Posts <- read.delim(file = readFromfileName, header = T, quote = "", sep = "\t")
  
  cat("Read", nrow(Posts), "rows from", readFromfileName, "\n")
  
  if(! missing(year) & ! missing(quarter)){
    Posts$create_ts <- as.POSIXlt(Posts$create_ts)
    Posts <- subset(Posts, quarters(Posts$create_ts) == paste("Q", quarter, sep = "")
                    & Posts$create_ts$year == (year - 1900))
    cat("Kept", nrow(Posts), "rows\n")
  }else{
    stop("Provide year and quarter")
  }
  
  return( doCorpCreation(Posts) )
}


## This function reads the data from csv file, cleans it and returns tm text Corpus
## readFromfileName - name of the file to read from
## Monthly Filter (optional):
##   year  - optional: year when docs were created, format = YYYY
##   month - optional: month when docs were created
createCorp <- function(readFromfileName, year, month){
  library(tm)  
  ## Load the data
  Posts <- read.delim(file = readFromfileName, header = T, quote = "", sep = "\t")
  
  cat("Read", nrow(Posts), "rows from", readFromfileName, "\n")
  
  if(! missing(year) & ! missing(month)){
    Posts$create_ts <- as.POSIXlt(Posts$create_ts)
    Posts <- subset(Posts, Posts$create_ts$mon == (month - 1) 
           & Posts$create_ts$year == (year - 1900))
    cat("Kept", nrow(Posts), "rows\n")
  }
  
  return( doCorpCreation(Posts) )
}

## this is a private function for corpus creation
doCorpCreation <- function(Posts){
  # Concatenate columns, otherwise DataframeSource gets confused
  Posts <- data.frame(paste(Posts$title, Posts$body, Posts$answers, sep =" "))
  
  ## Build a corpus
  corp <- Corpus(DataframeSource(Posts))
  
  ## Transform data
  corp = tm_map(corp, content_transformer(tolower)) #converting to lower case
  corp = tm_map(corp, removeWords, stopwords('english')) # Remove Stopwords
  corp = tm_map(corp, stemDocument) # Stemming
  corp = tm_map(corp, stripWhitespace) # Eliminate whitespace char
  
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
