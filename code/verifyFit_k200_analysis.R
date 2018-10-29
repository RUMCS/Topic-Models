#-------------------------------------------------------------------------
#set to True if you want to enable filtering of K > 200, else set to False
#plots the RMSE boxplots per model and per dataset name
#-------------------------------------------------------------------------
data_filter_200 <- F

if(data_filter_200){
  model_file_name <- "./models/models_fitting.remove_top_25_percent_and_values_gt_200.rda"
  models_performance_file_name <- "./models/models_performance.remove_top_25_percent_and_values_gt_200.csv"
  rmse_per_model_file_name <- "./models/models_performance__remove_top_25_percent_and_values_gt_200.pdf"
  rmse_per_dataset_file_name <- "./models/models_performance__per_dataset__remove_top_25_percent_and_values_gt_200.pdf"
  rmse_per_x_file_name <- "./models/models_performance__per_x__remove_top_25_percent_and_values_gt_200.pdf"
  rmse_per_n_file_name <- "./models/models_performance__per_n__remove_top_25_percent_and_values_gt_200.pdf"
}else{
  model_file_name <- "./models/models_fitting.remove_top_25_percent.rda"
  models_performance_file_name <- "./models/models_performance.remove_top_25_percent.csv"
  rmse_per_model_file_name <- "./models/models_performance__remove_top_25_percent.pdf"
  rmse_per_dataset_file_name <- "./models/models_performance__per_dataset__remove_top_25_percent.pdf"
  rmse_per_x_file_name <- "./models/models_performance__per_x__remove_top_25_percent.pdf"
  rmse_per_n_file_name <- "./models/models_performance__per_n__remove_top_25_percent.pdf"
  
}


dat.stats <- read.csv(models_performance_file_name)

#let's visualize performance of the models
library(tidyr)
dat.stats.long <- gather(dat.stats, model_name, rmse, rmse_flex_model_simple:rmse_tailored_constr_lm, factor_key=TRUE)

#reorder for box-plot
dat.stats.long$model_name <- factor(dat.stats.long$model_name, levels=c("rmse_flex_model_simple", "rmse_flex_model_complex", "rmse_constraint_model_simple", "rmse_constraint_model_complex",  "rmse_tailored_flex_lm", "rmse_tailored_constr_lm" ))

pdf(rmse_per_model_file_name)
boxplot(rmse ~ model_name, 
        data = dat.stats.long,
        log = "y",
        #xlab = "Model Name",
        las = 2,
        ylab = "RMSE",
        par(mar = c(12, 5, 4, 2)+ 0.1),
        names = c("Flex - simple", "Flex - complex", "Constr. - simple", "Constr. - complex", "Ind. fit - flex", "Ind. fit - constr")
)
grid()
dev.off()

pdf(rmse_per_dataset_file_name)
dat.stats$ds_merged_name <- as.factor(paste(dat.stats$dataset_name, "-", dat.stats$dt_partition))
boxplot(rmse_constraint_model_complex ~ ds_merged_name, 
        data = dat.stats,
        log = "y",
        #xlab = "Model Name",
        las = 2,
        ylab = "RMSE",
        par(mar = c(12, 5, 4, 2)+ 0.1)
        #names = c("Flex - simple", "Flex - complex", "Constr. - simple", "Constr. - complex", "Ind. fit - flex", "Ind. fit - constr")
)
grid()
dev.off()

pdf(rmse_per_dataset_file_name)
dat.stats$ds_merged_name <- as.factor(paste(dat.stats$dataset_name, "-", dat.stats$dt_partition))
boxplot(rmse_constraint_model_complex ~ ds_merged_name, 
        data = dat.stats,
        log = "y",
        #xlab = "Model Name",
        main = "Model - constrained, complex",
        las = 2,
        ylab = "RMSE",
        par(mar = c(12, 5, 4, 2)+ 0.1)
        #names = c("Flex - simple", "Flex - complex", "Constr. - simple", "Constr. - complex", "Ind. fit - flex", "Ind. fit - constr")
)
grid()
dev.off()

pdf(rmse_per_x_file_name)
boxplot(rmse_constraint_model_complex ~ X, 
        data = dat.stats,
        log = "y",
        xlab = "X",
        main = "Model - constrained, complex",
        #las = 2,
        ylab = "RMSE"
        #par(mar = c(12, 5, 4, 2)+ 0.1)
        #names = c("Flex - simple", "Flex - complex", "Constr. - simple", "Constr. - complex", "Ind. fit - flex", "Ind. fit - constr")
)
grid()
dev.off()

pdf(rmse_per_n_file_name)
boxplot(rmse_constraint_model_complex ~ N, 
        data = dat.stats,
        log = "y",
        xlab = "N",
        main = "Model - constrained, complex",
        #las = 2,
        ylab = "RMSE"
        #par(mar = c(12, 5, 4, 2)+ 0.1)
        #names = c("Flex - simple", "Flex - complex", "Constr. - simple", "Constr. - complex", "Ind. fit - flex", "Ind. fit - constr")
)
grid()
dev.off()


#summary stats
library(stargazer)
stargazer(dat.stats[,6:11])
