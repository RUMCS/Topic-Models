##----------------------------------------------------------------------------
## The script reads posts in delimited format, and generates a list of keyword
## Usage: Rscript my_Script_name.R out_directory/ in_file_name yyyy
## This program includes the probability per keyword 
##        also includes the specificity by keyword by topic  
##----------------------------------------------------------------------------
library(tm)
library(foreach)
library(doParallel)
library(stringr) ## added to handle string functions

memory.limit(24356)
memory.size(max = TRUE)

source("utilsdb2.R")
xURL <- "http://dba.stackexchange.com/questions/"
xSystime <- format(Sys.time(), "%a-%b-%d_%H-%M-%S_%Y")
## following lines commented out to run in PC ##
##args <- commandArgs(trailingOnly = TRUE)
##outDir <- args[1]
##readFrom <- paste(args[1], args[2],sep="")
##year  <-  as.integer(args[3])


## following lines hardcoded to run in PC ##
readFrom <- "nondb2.awk" ##"Posts.xml.filtered.csv.new" ##"db2.awk"
year <- 2015
##
##month <-  as.integer(args[3])

topKeywordCount <- 20 ## SET TO 20 FOR THIS RUN

## corp <- createCorp(readFrom)
corp <- createCorp(readFrom, year) ## added year

# To access original id of a document in VCorp run corp[[ index_of_document_1_to_N ]]$meta$id
# To access original id of a document in DocumentTermMatrix run dtm$dimnames$Docs[ index_of_document_1_to_N ]

cat("DEBUG after calling createCorp!\n")
## Build a Document-Term Matrix
dtm <- DocumentTermMatrix(corp, control = list(minWordLength = 2)) #keep words of lenght 2 or longer
cat("Before tf-idf: term count =", ncol(dtm), ", doc count =", nrow(dtm), "\n")#dtm <- removeFrequentWords(dtm) #removing based on median tf-idf value
cat("After tf-idf: term count =", ncol(dtm), ", doc count =", nrow(dtm), "\n"nlikely event)
cat("After removing terms appearing only in 1 document: term count =", ncol(dtm), ", doc count =", nrow(dtm), "\n")

#setup parallel backend to use 8 processors
cl<-makeCluster(8) ## change to 8 for multiprocessing
registerDoParallel(cl)

foreach(topicCount = c(5) #max = 1 topic per document ## changed topic count=5
        , .packages='topicmodels' #include package
) %do%
{ #change to %dopar% for multi execution
  cat("*RUNNING SINGLE PROCESSING*\n")
  cat("TopicCount:", topicCount, "\n") #to screen (no screen output is parallel mode)  
  mdl <- LDA(dtm, topicCount) #LDA model 
  topic.keyword <- terms(mdl, topKeywordCount)
  mdl.post <- posterior(mdl) #get posterior data
} 

## the following nested loop was added to append the p-value per topic.keyword
Ttopic.keyword <- topic.keyword ## make a copy

## the following nested loop was added to write the p-values of topics by column

outPval <- paste("pvalues", xSystime, "_", readFrom, "_", year, "-", topicCount, ".txt", sep = "")
for (j in 1:topicCount) {
  cat("TOPIC ", j, "\n", append = T, 
      file = outPval)
  cat("KEYWORD,PROBABILITY\n", append = T,
      file = outPval)
  for (i in 1:topKeywordCount){

    Ttopic.keyword[i,j] <- 
      Ttopic.keyword[i,j] <- paste(topic.keyword[i,j], ",", mdl.post$terms[[j, topic.keyword[i,j]]])
    cat(Ttopic.keyword[[i,j]],"\n", append = T, 
        file = outPval)
  }
}

write.csv(Ttopic.keyword, 
          file = paste("pvalues_", xSystime, "_", readFrom, "_", year, "-", topicCount, ".csv", sep = "") )

## code added to support report creation

## Will retrieve all ids and questions
t <- topics(mdl)
len <- length(topics(mdl))
saveTo <- paste("Report-", xSystime, "_", readFrom, "-", year, ".txt")
cat("question_id\tquestion_description\ttopic\t\n", sep="\t", append = T, file = saveTo) # to file
for (i in 1:len){
  ri <- corp[[i]]$meta$id
  rq <- corp[[i]]$meta$quest
  rt <- t[i]
  r <- paste(xURL, ri, "\t", rq, "\t", rt, sep="")
  cat(r, "\n", append = T, file = saveTo) # to file
}

## The following code was added for generating the specificity per keyword

## t is number of terms in dtm
## tC is the index of topicCount
## kC in the index of topKeywordCount
## topKeywordCount is the number of keywords
## dC is the index of number of topics
## countKw is the counter of the keyword that appears in the topics

allDocuments <- nrow(dtm) ## Size in rows of dtm
allKeywords <- ncol(dtm)
countKw <- matrix(0, topKeywordCount, topicCount) ## initialize count keyword matrix with 0
 
## code added to calculate frequencies 
countKw <- array(0,dim=c(topKeywordCount, topicCount)) ## intermediate array to handle the counting of kW in allDocuments
t <- topics(mdl) ## documents and the topics it belongs to
z <- as.matrix(dtm) ## dtm as matrix
dfz <- as.data.frame(z) ## converts z to data frame
tK <- 0 ## array(0,dim=c(topKeywordCount, topicCount))

acumTo <- rep(0,topicCount) ## initialize acum of topics
for(i in 1:topKeywordCount){
     acumKw<-0
     for(j in 1:topicCount){
      xlookup <- topic.keyword[i,j] 
      colZ <- which(colnames(z)==topic.keyword[i,j]) ## finds the column in z that the kW is in
      
      for(k in 1:allDocuments){
        if(dfz[[colZ]][k] > 0) {
         tK <- t[[k]]  ## indek tk is the topic number
         
         if(tK==j){
           if(dfz[[colZ]][k] > 0){  
              acumKw <- acumKw + 1
           }
         }
        } 
      }
      
      if(!is.null(acumKw)){
          countKw[i,j] <- acumKw
          acumKw <-0
      } 
     }  
   
} 

## Calculation of N as the number of documents per topic
for(i in 1:nrow(dtm)){
  topic<-t[[i]]
  acumTo[topic]<-acumTo[topic]+1
}

## the following nested loop was added to write the topics by column
outIdf <- paste("idf_", xSystime, "_", readFrom, "_", year, "-", topicCount, ".txt", sep = "")
for (j in 1:topicCount) {
  cat("TOPIC ", j, "\n", append = T, 
      file = outIdf)
  cat("KEYWORD,PROBABILITY,N,d,idf\n", append = T,
      file = outIdf)
  for (i in 1:topKeywordCount){
 
    idf<-log10(acumTo[j]/countKw[i, j]) 
    cat("i j", i, " ", j, "idf= log(", acumTo[j], "/", countKw[i,j], ")\n")
    Ttopic.keyword[i,j] <- 
      paste(topic.keyword[i,j], ",", mdl.post$terms[[j, topic.keyword[i,j]]], ",",
            acumTo[j], ",", countKw[i, j], ",", idf) ## going to work with a copy of topic.keywords instead
    cat(Ttopic.keyword[[i,j]],"\n", append = T, 
        file = outIdf)
  }
}

cat("Done\n")
