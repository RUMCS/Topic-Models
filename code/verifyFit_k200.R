##----------------------------------------------------------------------------
## This program iterates over the specified datasets, and computes the RMSE of the approximation formulas
##----------------------------------------------------------------------------

#set to True if you want to enable filtering of K > 200, else set to False
data_filter_200 <- T

readFromfileName = "topX.csv.original.zip"

if(data_filter_200){
  model_file_name <- "./models/models_fitting.remove_top_25_percent_and_values_gt_200.rda"
  models_performance_file_name <- "./models/models_performance.remove_top_25_percent_and_values_gt_200.csv"
}else{
  model_file_name <- "./models/models_fitting.remove_top_25_percent.rda"
  models_performance_file_name <- "./models/models_performance.remove_top_25_percent.csv"
}

load(model_file_name)
############################################################################
#####              define approximation formulas                       #####
############################################################################

# The "flexible" model a(X, N) * K^b(X,N)
# K == topic count
# X == the top-X 
# N == number of documents
model_two_param_simple <- function(K, X, N){
  
  val <- data.frame(X = X, N = N)
  a <- predict(dat.lm.intercept, val) # 1.171 + 0.649 * log(X) - 0.116 * log(N)
  b <- predict(dat.lm.slope, val)  # -0.841 + 0.002 * X + 0.00003 * N
  
  exp(a  + b * log(K))
}

model_two_param_complex <- function(K, X, N){
  
  val <- data.frame(X = X, N = N)
  a <- predict(dat.lm.intercept.complex, val) 
  b <- predict(dat.lm.slope.complex, val) 
  
  exp(a  + b * log(K))
}

# The "constrained" model X^-b(X, N) * K^b(X,N)
# K == topic count
# X == the top-X 
# N == number of documents

model_one_param_simple <- function(K, X, N){
  val <- data.frame(X = X, N = N)
  b <- predict(dat.lm.b, val)  #-8.408351E-01 + 2.110667E-03 * X + 3.397901E-05 * N
  
  X^-b * K^b
}

model_one_param_complex <- function(K, X, N){
  val <- data.frame(X = X, N = N)
  b <- predict(dat.lm.b.complex, val)
  
  X^-b * K^b
}


#############################################################################
#####                        helpers                                    #####
#############################################################################

# root-mean-square error 
rmse <- function(actual, expected){
  sqrt(mean((expected-actual)^2))
}

############################################################################
#####                             main                                 #####
############################################################################

dat <- read.csv(unz(readFromfileName, "topX.csv.original"), header = T, sep = ",", row.names = NULL) ## NOT EXCEL
validGroups <- unique(paste(dat$timeframe_type, dat$timeframe, dat$dataset_name, sep = " "))

dat.stats <- data.frame()
for(iG in 1:length(validGroups)) {
  s <- strsplit(validGroups[iG], " ")
  sdf <- as.data.frame(s)
  p1 <- as.character(sdf[1,]) ## timeframe_type
  p2 <- as.character(sdf[2,]) ## timeframe
  p3 <- as.character(sdf[3,]) ## dataset_name
  
  cat("p1=", p1, " p2=", p2, " p3=", p3, " iG=", iG, "\n")
  for(itopXX in 2:50) {
    ds <- dat[dat$timeframe_type ==  p1 & dat$timeframe == p2 & dat$dataset_name == p3 & dat$topXX == itopXX,]
    #cat("p1=", p1, " p2=", p2, " p3=", p3, " itopXX=", itopXX, " iG=", iG, "\n")
 
    if(data_filter_200){ #if True -- keep only values of K <= 200
      ds <- ds[ds$topicCount <= 200, ]
    }
    ds <- ds[ds$topicCount < 0.75 * max(ds$documentCount), ] 
    
    
    #fit indvidual flex model
    ds.lm.flex <- lm(log(postFraction) ~ log(topicCount), data = ds )
    
    #fit indvidual constrained model
    ds.lm.constr <- nls(ds$postFraction ~ ds$topXX^(-b) * ds$topicCount^b, data = ds, start = list( b = -1))
    
    
    F.ds.lm.flex <- exp(predict(ds.lm.flex)) #keep in mind that we are operating on log-transform data here, need to convert it back
    F.ds.lm.constr <- predict(ds.lm.constr)
  
    F.model_two_param_simple <- model_two_param_simple(ds$topicCount, ds$topXX[1], ds$documentCount[1])
    F.model_one_param_simple <- model_one_param_simple(ds$topicCount, ds$topXX[1], ds$documentCount[1])
    F.model_two_param_complex <- model_two_param_complex(ds$topicCount, ds$topXX[1], ds$documentCount[1])
    F.model_one_param_complex <- model_one_param_complex(ds$topicCount, ds$topXX[1], ds$documentCount[1])
    
    #compute rmse and save stats
    dat.stats <- rbind(dat.stats,
      data.frame(
        dt_partition = p1,
        time_frame = p2,
        dataset_name = p3,
        X = ds$topXX[1],
        N = ds$documentCount[1],
        rmse_flex_model_simple = rmse(ds$postFraction, F.model_two_param_simple),
        rmse_constraint_model_simple = rmse(ds$postFraction, F.model_one_param_simple),
        rmse_flex_model_complex = rmse(ds$postFraction, F.model_two_param_complex),
        rmse_constraint_model_complex = rmse(ds$postFraction, F.model_one_param_complex),
        rmse_tailored_flex_lm = rmse(ds$postFraction, F.ds.lm.flex),
        rmse_tailored_constr_lm = rmse(ds$postFraction, F.ds.lm.constr)
      )
    )
  }
}

write.csv(dat.stats, models_performance_file_name, row.names = FALSE)
