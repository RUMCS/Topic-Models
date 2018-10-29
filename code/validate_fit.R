#--------------------------------------------------------------------------
# set to True if you want to enable filtering of K > 200, else set to False
# performs validation using 10 folds
#--------------------------------------------------------------------------
data_filter_200 <- F

# root-mean-square error 
rmse <- function(actual, expected){
  sqrt(mean((expected-actual)^2))
}

# one param complex prediction model
model_one_param_complex <- function(K, X, N){
  val <- data.frame(X = X, N = N)
  b <- predict(dat.lm.b.complex, val)
  
  X^-b * K^b
}


###########################################################################
#####                 get raw data                                    #####
###########################################################################

readFromfileName = "topX.csv.original.zip"
dat.raw <- read.csv(unz(readFromfileName, "topX.csv.original"), header = T, sep = ",", row.names = NULL)

#add key column
dat.raw$data_subset_name <- paste(dat.raw$timeframe_type, dat.raw$timeframe, dat.raw$dataset_name)

###########################################################################
#####                 get the individual fit of b-values              #####
###########################################################################

if(data_filter_200){
  to_fit_file_name <- "to_fit.remove_top_25_percent_and_values_gt_200.csv"
}else{
  to_fit_file_name <- "to_fit.remove_top_25_percent.csv"
}
dat <- read.csv(to_fit_file_name)

#add key column
dat$data_subset_name <- paste(dat$time_interval, dat$time_frame, dat$dataset_name)



###########################################################################
#####                 do run 10-fold validation                       #####
###########################################################################

#split 96 datasets into 10 folds
unique_data_subsets <- unique(dat$data_subset_name)
require(caret)
folds <- createFolds(c(1:length(unique_data_subsets)), k = 10, list = TRUE, returnTrain = FALSE)

rmse_per_fold <- c()
for(f in 1:10){
  train_subsets <- unique_data_subsets[-folds[[f]]]
  test_subsets  <- unique_data_subsets[ folds[[f]]]
  
  
  #let's select the b values associated with the data subsets in the train set
  train_data <- dat[ dat$data_subset_name %in% train_subsets, ]
  
  #complex coanstrained model: compute the model for the value of b
  dat.lm.b.complex <- lm(b ~ X + N + log(X) , data = train_data)
  
  #let us see how the model for b value fits the raw data  for the test set
  dat.raw.test <- dat.raw[dat.raw$data_subset_name %in% test_subsets, ]
  
  #apply DS1 or DS2 filter
  if(data_filter_200){ #if True -- keep only values of K <= 200
    dat.raw.test <- dat.raw.test[dat.raw.test$topicCount <= 200, ]
  }
  dat.raw.test <- dat.raw.test[dat.raw.test$topicCount < 0.75 * dat.raw.test$documentCount, ] 
  
  
  fitted_values <- model_one_param_complex(dat.raw.test$topicCount, dat.raw.test$topXX, dat.raw.test$documentCount)
  #compute and save rmse per fold
  rmse_per_fold <- c(rmse_per_fold, rmse(dat.raw.test$postFraction, fitted_values))
}


###########################################################################
#####                 show summary stats                              #####
###########################################################################

summary(rmse_per_fold)


if(data_filter_200){
  figure_file_name <- "./figures/one_param_complex_10-fold_validation_DS2.pdf"
  rmse.df.ds2 <- data.frame(fiter_name = "Dataset 2", rmse = rmse_per_fold)  
}else{
  figure_file_name <- "./figures/one_param_complex_10-fold_validation_DS1.pdf"
  rmse.df.ds1 <- data.frame(fiter_name = "Dataset 1", rmse = rmse_per_fold)  
}

pdf(figure_file_name)
boxplot(rmse_per_fold, 
        ylab = "RMSE",
        xlab = "Constr. - complex"
        )
grid()
dev.off()

#need to run the script twice with data_filter_200 = T and data_filter_200 = F to get the figure below correctly
pdf("./figures/one_param_complex_10-fold_validation_DS_both.pdf")
rmse.df.ds_both <- rbind(rmse.df.ds1, rmse.df.ds2)
plot(rmse.df.ds_both$fiter_name, rmse.df.ds_both$rmse,
        ylab = "RMSE",
        xlab = "Constr. - complex"
)
grid()
dev.off()
