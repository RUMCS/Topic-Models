##-------------------------------------------------------------------------------
## The script reads posts in delimited format, and generates 
## distribution of LDA topics for a given month year,
## saving the output in in_file_name + ".year-month.topic_frequency"
## Usage: Rscript do_lda_analysis.R in_file_name year month
##
## This version adds a Semaphore in order to slowdown the process of writing in a parallel processing
## includes refactoring changes (i.e. taking out semaphore code from for loop..)
##------------------------------------------------------------------------------
library(tm)
library(foreach)
library(doParallel)

source("utils.R")

args <- commandArgs(trailingOnly = TRUE)
readFrom <- args[1]
year  <-  as.integer(args[2]) 
month <-  as.integer(args[3])

semaphoreFileName <- paste("semaf_", Sys.getpid(), sep='') # append R's process id to semaphore to avoid conflict
if( file.exists(semaphoreFileName) ){ file.remove(semaphoreFileName) } 

saveTo <- paste(readFrom, ".", year, "-", month, ".topic_frequency", sep ="")

corp <- createCorp(readFrom, year, month)

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

foreach(topicCount = 2:nrow(dtm) #max = 1 topic per document
        , .packages='topicmodels' #include package
) %dopar% { #change to %do% for sequential execution
  startRun <- Sys.time()
  
  val <- getTopicsFrequency(dtm, topicCount)
  
  prefix <- paste(topicCount, val$mdl.alpha, val$mdl.beta.mean, val$mdl.beta.sd, difftime(Sys.time(), startRun, units = "secs"), sep="\t")
  
  cat("TopicCount:", topicCount, "\n") #to screen (no screen output is parallel mode)
  
  while (file.exists(semaphoreFileName)==TRUE) {
        Sys.sleep(1); cat("in while sleeping 1 sec\n");
  }
  file.create(semaphoreFileName) ## lock file
  for (i in 1:length(val$topic.frequency)){
       #cat("in cat for",i,Sys.time(), file.exists(semaphoreFileName), "\n") 
       cat(prefix, val$topic.frequency[i]
            , "\n", sep="\t", append = T
            , file = saveTo) # to file
  }
  file.remove(semaphoreFileName) ## unlock file
  #cat("Semaphore is",file.exists(semaphoreFileName), "\n")
 }
cat("Saved data to", saveTo, "\n")
cat("Done\n")
