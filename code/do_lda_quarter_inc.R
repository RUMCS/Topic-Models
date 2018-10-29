##--------------------------------------------------------------
## The script reads posts in delimited format, and generates 
## distribution of LDA topics for a given quarter-year,
## saving the output in in_file_name + "." + year "-q" + quarter ".topic_frequency"
## Usage: Rscript do_lda_quarter_inc.R in_file_name year quarter
##--------------------------------------------------------------
library(tm)
library(foreach)
library(doParallel)

source("utils.R")

args <- commandArgs(trailingOnly = TRUE)
readFrom <- args[1]
year  <-  as.integer(args[2])
quarter <-  as.integer(args[3])

semaphoreFileName <- paste("semaf_", Sys.getpid(), sep='') # append R's process id to semaphore to avoid conflict
if( file.exists(semaphoreFileName) ){ file.remove(semaphoreFileName) }

saveTo <- paste(readFrom, ".", year, "-q", quarter, ".topic_frequency", sep ="")

corp <- createCorpQuarter(readFrom, year, quarter)

## Build a Document-Term Matrix
dtm <- DocumentTermMatrix(corp, control = list(minWordLength = 2)) #keep words of lenght 2 or longer
cat("Before tf-idf: term count =", ncol(dtm), ", doc count =", nrow(dtm), "\n")
dtm <- removeFrequentWords(dtm) #removing based on median tf-idf value
cat("After tf-idf: term count =", ncol(dtm), ", doc count =", nrow(dtm), "\n")
## the following line commented out to avoid sparse error
dtm <- removeSparseTerms(dtm, 1 - (1.1/nrow(dtm)) )  #remove terms appearing only in 1 document
dtm <- dtm[row_sums(dtm) > 0,]  #remove docs that have no terms remaining (unlikely event)
cat("After removing terms appearing only in 1 document: term count =", ncol(dtm), ", doc count =", nrow(dtm), "\n")

#setup parallel backend to use 8 processors
cl<-makeCluster(8)
registerDoParallel(cl)

cat("topicCount\tmdl.alpha\tmdl.beta.mean\tmdl.beta.sd\ttime (sec)\ttopic.frequency\n", sep="\t", append = T
    , file = saveTo) # to file

vecCount <- incCalc(nrow(dtm))

## foreach(topicCount = 2:nrow(dtm) #max = 1 topic per document
cat("lenght(vecCount)=", length(vecCount), "\n")
foreach(ix = 2:length(vecCount) 
        , .packages='topicmodels' #include package
) %do% { #change to %dopar% for multiprocess execution
  
  topicCount <- vecCount[ix]
  
  startRun <- Sys.time()
  
  val <- getTopicsFrequency(dtm, topicCount)
  
  prefix <- paste(topicCount, val$mdl.alpha, val$mdl.beta.mean, val$mdl.beta.sd, difftime(Sys.time(), startRun, units = "secs"), sep="\t")
  
  cat("TopicCount:", topicCount, "\n") #to screen (no screen output is parallel mode)
  
  while (file.exists(semaphoreFileName)==TRUE) {
    Sys.sleep(1); cat("in while sleeping 1 sec\n");
  }
  file.create(semaphoreFileName) ## lock file
  for (i in 1:length(val$topic.frequency)){
    cat(prefix, val$topic.frequency[i]
        , "\n", sep="\t", append = T
        , file = saveTo) # to file
  }
  file.remove(semaphoreFileName)
}
cat("Saved data to", saveTo, "\n")
cat("Done\n")
