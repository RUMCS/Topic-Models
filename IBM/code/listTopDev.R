####################################################################################
## listTop1Dev.R
## This program reads idf files with format idf_<mmm>_<dd>_<yyyy1>_T<#####> where mmm_dd_yyyy1 is the today's date and T#### 
## is the topic number, retrieves the top 1 topic for each file-yyyy1 and writes it into a file for then being reported to
## IBM.
## The validation wether a value is NA is not exhaustive but partial
## Author. Jorge Lopez
## Date. 7/19/2016
###################################################################################


listf <- list.files(pattern="idf_")
outFileN <- "TopKeysN.txt"
year<-character
keys_year = vector()
key = vector()
topX <- data.frame(year = character(), keys_year = vector())
cat("YEAR,KEYWORD\n", append =T, file = outFileN)
maxKeys <- 30
keyNum <- 5
countK <- 0
i_listf <- 1

for (i_listf in 1:length(listf)) {
  
  #countK <- countK + 1
  #if(countK > keyNum) {
  #  countk <- 0
  #  i_listf <- i_listf + 5
  #}
  
  readFromfileName <- listf[i_listf]
  yyyy <- substr(readFromfileName,17,20)
  
  cat("listf[", i_listf, "]=", listf[i_listf],"\n")
  readFromfileName <- gsub(" ", "", readFromfileName, fixed = TRUE)
  cat("DEBUG: Reading=", readFromfileName, "\n")
  dat <- read.csv(readFromfileName, header = T, sep = ",", row.names = NULL)
  ##i_listf <- i_listf + 1
  
  i_key <- 1
  count <- 0
  topKeys <- 2
  keysXf <- 1
  
   fKeyNotOK <- function(dat, i_key) {
        str1 <- gsub(" ", "", as.character(dat$KEYWORD[i_key]), fixed = TRUE) ## trim keyword spaces
        cat("str1 is ", str1, "\n")
        NotOK <- grepl('^[A-Za-z]+$', str1) & str1 != "ffff"
        cat("NotOK is ", NotOK, "\n")
        keyNotOK <- NotOK
        cat("KeyNotOK is ", keyNotOK, "\n")
        ##NotOK <- !(dat$KEYWORD[i_key]==" ffff ")
        if(NotOK){return(FALSE)} else {return(TRUE)}
                 ##KeyNotOK <- (dat$KEYWORD[i_key]==" db2 ")| +
                 ##(dat$KEYWORD[i_key]==" char(254) ")| +
                 ##(dat$KEYWORD[i_key]==" 0000 ")| +
                 ##(dat$KEYWORD[i_key]==" >db2 ")| +
                 ##(dat$KEYWORD[i_key]==" 11/26/2010 ")| +
                 ##(dat$KEYWORD[i_key]==" >>> ")| +
                 ##(dat$KEYWORD[i_key]==" --| ")| +
                 ##(dat$KEYWORD[i_key]==" node:000 ")| +
                 ##(dat$KEYWORD[i_key]==" >>>> ")| +
                 ##(dat$KEYWORD[i_key]==" 07/03/2008 ")| +
                 ##(dat$KEYWORD[i_key]==" 0.015508423678423 ")| +
                 ##(dat$KEYWORD[i_key]==" ................ ")| +
                 ##(dat$KEYWORD[i_key]==" 0.0376739033577241 ")| +
                 ##(dat$KEYWORD[i_key]==" 0.0310868142555237 ")| +
                 ##(dat$KEYWORD[i_key]==" instance:db2 ")| +
                 ##(dat$KEYWORD[i_key]==" tid:8192 ")| +
                 ##(dat$KEYWORD[i_key]==" database:mmmmm ")| +
                 ##(dat$KEYWORD[i_key]==" >0000 ")| +
                 ##(dat$KEYWORD[i_key]==" 0.015508423678423 ")| +
                 ##(is.na(dat$KEYWORD[i_key]))
     ###return(KeyNotOK)
     }
  
  while(fKeyNotOK(dat, i_key))
  {
    cat("The Keyword found is db2,NA or invalid! \n")
    i_key <- i_key + 1 
    count <- count + 1
    if (count > maxKeys) { 
        cat("There are > maxKeys NA or db2 keywords in this file, file will not be processed\n")
        break
    }
  }
    
  ##  key[i_key] <- as.character(dat$KEYWORD[i_key])
  ##  i_key <- i_key + 1 
  
  cat("$$$$$$$$$$$$ ikey is:", i_key)
  cat(" fkeynotok is:", fKeyNotOK(dat, i_key))
  cat("key[i_key]=",key[i_key],"\n")
  
  if(!fKeyNotOK(dat, i_key)){  
     cat("i_key=", i_key, "\n")
     key[i_key] <- as.character(dat$KEYWORD[i_key])
     ##i_key <- i_key + 1
     ##for (i_key in 1:keysXf) {
       
       cat(yyyy,",", key[i_key], "\n", append =T, file = outFileN)
     ##}
  }
}
