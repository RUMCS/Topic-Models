##------------------------------------------------------------------------------
### This program reads the *.topic.frequency files for each dataset gets the top frequency for x=2..50, and merges them up together 
### for them to be analyzed later with a regression model
##------------------------------------------------------------------------------
readFile <- function(fName) {
  dat <- read.delim(file = fName, header = T, fill=T, sep = "\t",
                    row.names=NULL)
  colnames(dat)<-c(colnames(dat)[-1],"x")
  dat$x<-NULL
  dat<-as.data.frame(dat)
  return(dat)
}

processFile <- function(dat){
  topicCountList <- sort(unique(strtoi(dat$topicCount))) ## making sure the list is ordered
  
  timeframe_type = vector()
  timeframe = vector()
  dataset_name = vector()
  topicCount = vector()
  postFraction = vector()
  topXX = vector()
  documentCount = vector()
  
  topX <- data.frame(timeframe_type = vector(), timeframe = vector(), dataset_name = vector(), topicCount = vector(), postFraction = vector(), topXX = vector(), documentCount = vector())
  
  for(x in 2:50){#top X topics
    for( t in topicCountList ){ ## strtoi(topicCount)
      if (t < x){ ## to avoid writing NA's
        next
      }
      else
      {
        d <- dat [which(dat$topicCount == t), ]
        freq <- sort(d$topic.frequency, decreasing = T)
        documentCount <- sum (freq) #number of docs in corpus
        postFraction = sum(freq[1:x]) / documentCount # this is an F value for a given top X
          if (is.na(postFraction)) {
            cat("====================== postfraction is NA program skips interaction t=", t, "===============\n")
            cat("x=", x, " documentCount=", documentCount, "\n")
            next
          }
        topX <- rbind(topX, data.frame(timeframe_type = Xtimeframe_type[i_dsName], timeframe = Xtimeframe[i_listf], dataset_name = Xdataset_name[i_dsName], topicCount = t, postFraction = postFraction, topXX = x, documentCount = documentCount))
      }
    }
  }
  ## writes output file everytime finishes processing one frequency file
  
  fName <- paste(outPath, "topX", ".csv", sep = "")
  
  if(firstime){
    cat("** WRITING OUTPUT FILE FIRSTIME**", firstime, "\n")
    cat(dir(outPath),"\n")
    ## cat("FILE EXIST NOT NEGATED:", fName, " is ", file.exists(fName), "\n")
    ## cat("FILE EXIST NEGATED:", fName, " is ", !file.exists(fName), "\n")
    write.table(topX, file = fName, sep = ",", row.names=FALSE, col.names=TRUE, append = T)
    firstime <<- FALSE
    print(file.info(fName))
  } else {
    cat("** WRITING OUTPUT FILE SUCCESIVE **", firstime, "\n")
    write.table(topX, file = fName, sep = ",", row.names=FALSE, col.names=FALSE, append = T)
    print(file.info(fName))
  }
} ## EOF


##############################################################
#####                    main                            #####
##############################################################

XSystime <<- format(Sys.time(), "%a-%b-%d_%H-%M-%S_%Y")
dsName <- c("salesforce-m-", "dba-m-", "dba-q-", "android-m-") ## to handle the ds names

### test dsName <- c("test-q-") ## to handle the ds names

Xtimeframe_type <<- c("monthly", "monthly", "quarterly", "monthly")
Xdataset_name <<- c("salesforce", "dba", "dba", "android")

FKdir <- c()
outListf <- c() ### this is the output directory of the files to be plotted, may not be used
ffDir <- c()
i_dsName <<- 1
firstime <<- TRUE ## this flag was placed after first, for so it wrote the hdrs each time ds changed!
for (i_dsName in 1:length(dsName)) {
  cat("DEBUG: Processing dataset=", dsName[i_dsName], " at:", XSystime, "\n")
  
  ## former FKdir[i_dsName] <- file.path("~", "Thesis", "Fitting", paste(dsName[i_dsName],"FK", sep=""), "files-FK")
  ## former ffDir[i_dsName] <- file.path("~", "Thesis", "Fitting", paste(dsName[i_dsName],"FK", sep=""), "files-frequencies")
  
  FKdir[i_dsName] <- file.path("/media", "data", "thesis", "Thesis", "Fitting", paste(dsName[i_dsName],"FK", sep=""), "files-FK")
  ffDir[i_dsName] <- file.path("/media", "data", "thesis", "Thesis", "Fitting", paste(dsName[i_dsName],"FK", sep=""), "files-frequencies")
  
  
  setwd(ffDir[i_dsName])
  cat("DEBUG: ffDir[", i_dsName, "]=", ffDir[i_dsName], "\n")
  
  listf <- list.files(pattern="Posts.xml.csv")
  
  ## former outPath <<- file.path("~", "Thesis", "Fitting","/")
  outPath <<- file.path("/media", "data", "thesis", "Thesis", "Fitting","/")
  ##} ## ends here for debugging purposes next code deactivated:
  ##  for (i_dsName in 1:length(dsName)) {
  
  for (i_listf in 1:length(listf)) {
    readFromfileName <- listf[i_listf]
    cat("listf[", i_listf, "]=", listf[i_listf],"\n")
    Xtimeframe <<- vector()
    Xtimeframe[i_listf] <- substr(gsub("\\.","-",readFromfileName),15,21)
    readFromfileName <- gsub(" ", "", readFromfileName, fixed = TRUE)
    cat("DEBUG: Reading=", readFromfileName, " Xtimeframe=", Xtimeframe[i_listf], "\n")
    dat <- readFile(readFromfileName)
               processFile(dat)
  }
}
cat("** END OF EXECUTION **\n")
