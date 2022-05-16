##-----------------------------------------------------##
##          Part 4: Clusterwise Multiblock             ##
##-----------------------------------------------------##

# Note: Experimental. Unreliable block.

# load library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, zoo, lubridate, tidyverse, fastDummies, rsample, mbclusterwise)

# set the current working directory the same as shell
wd_path = Sys.getenv("PWD") 
if(getwd() != wd_path) setwd(wd_path)

# set reference path 
if(!file.exists(".here")) here::set_here(path= ".")

# load functions
source(here::here("R/utils.R"))

# set system time zone
Sys.setenv(TZ='Europe/Zurich')

# load data
load(here::here('analysis/mbpls_fit_filtered.Rdata'))

# use data prep for Züchter & Mäster
mbdat1 <- züchter_1 # mäster_1

# choosing 75% of the data to be the training data
data_split <- initial_split(mbdat1, prop = .75)

# extracting training data and test data as two separate dataframes
data_train <- training(data_split)
data_test  <- testing(data_split)

# define dependent and independent
dep <- c("lagedesbetriebes", "abschmrngdsstswg", "arbtsblfmngmtdsb")
indep <- setdiff(names(mbdat1), c(dep, "betriebsid", "tvdnummer"))

# X & Y Train
Y_train <- data_train[, dep]
X_train <- data_train[, indep]

block <- c(7,7,6,6) # ref: block_zht 02_multiblock_analysis.Rmd line 333

# Use all data

# multiblock 
# X & Y all
Y_all <- mbdat1[, dep]
X_all <- mbdat1[, indep]

# Multiblock PLS & RA G=1

res1.g1 <- list()
res2.g1 <- list()

for (H in c(1:20)){
  print(paste("G = 1 H=", H, sep=""))
  res1.g1[[H]] <- cw.tenfold(Y = Y_all, X = X_all, blo = block, option = "none", G = 1, H,
                          FOLD = 10, INIT = 20, method = "mbpcaiv", Gamma = NULL, parallel.level = "high");
  gc();
  res2.g1[[H]] <- cw.tenfold(Y = Y_all, X = X_all, blo = block, option = "uniform", G = 1, H,
                          FOLD = 10, INIT = 20, method = "mbpls", Gamma = NULL, parallel.level = "high");
  gc();
}
rm(H)
# 1
res1.g1.cal <- unlist(lapply(1:20, function(x) mean(sqrt(res1.g1[[x]]$sqrmse.cal), na.rm=TRUE)))
res1.g1.val <- unlist(lapply(1:20, function(x) mean(sqrt(res1.g1[[x]]$sqrmse.val), na.rm=TRUE)))
# 2
res2.g1.cal <- unlist(lapply(1:20, function(x) mean(sqrt(res2.g1[[x]]$sqrmse.cal), na.rm=TRUE)))
res2.g1.val <- unlist(lapply(1:20, function(x) mean(sqrt(res2.g1[[x]]$sqrmse.val), na.rm=TRUE)))

rmse.g1.cal <- rbind(res1.g1.cal, res2.g1.cal)
rmse.g1.val <- rbind(res1.g1.val, res2.g1.val)

rownames(rmse.g1.cal) <- rownames(rmse.g1.val) <- c("mbpcaiv", "mbpls")
colnames(rmse.g1.cal) <- colnames(rmse.g1.val) <- paste("H", 1:20, sep = "=")

save(res1.g1, res1.g1.cal, res1.g1.val, 
     res2.g1, res2.g1.cal, res2.g1.val, 
     rmse.g1.cal, rmse.g1.val, 
     file = here::here("analysis/mbfit_cvg1_none.RData"))

rmse.g1.cal
rmse.g1.val

# Multiblock PLS & RA G=2

res1.g2 <- list()
res2.g2 <- list()

for (H in c(1:20)){
  print(paste("G = 2 H=", H, sep=""))
  res1.g2[[H]] <- cw.tenfold(Y = Y_all, X = X_all, blo = block, option = "none", G = 2, H,
                             FOLD = 10, INIT = 20, method = "mbpcaiv", Gamma = NULL, parallel.level = "high");
  gc();
  res2.g2[[H]] <- cw.tenfold(Y = Y_all, X = X_all, blo = block, option = "none", G = 2, H,
                             FOLD = 10, INIT = 20, method = "mbpls", Gamma = NULL, parallel.level = "high");
  gc();
}
rm(H)
# 1
res1.g2.cal <- unlist(lapply(1:20, function(x) mean(sqrt(res1.g2[[x]]$sqrmse.cal), na.rm=TRUE)))
res1.g2.val <- unlist(lapply(1:20, function(x) mean(sqrt(res1.g2[[x]]$sqrmse.val), na.rm=TRUE)))
# 2
res2.g2.cal <- unlist(lapply(1:20, function(x) mean(sqrt(res2.g2[[x]]$sqrmse.cal), na.rm=TRUE)))
res2.g2.val <- unlist(lapply(1:20, function(x) mean(sqrt(res2.g2[[x]]$sqrmse.val), na.rm=TRUE)))

rmse.g2.cal <- rbind(res1.g2.cal, res2.g2.cal)
rmse.g2.val <- rbind(res1.g2.val, res2.g2.val)

rownames(rmse.g2.cal) <- rownames(rmse.g2.val) <-  c("mbpcaiv", "mbpls")
colnames(rmse.g2.cal) <- colnames(rmse.g2.val) <- paste("H", 1:20, sep = "=")

save(res1.g2, res1.g2.cal, res1.g2.val,  
     res2.g2, res2.g2.cal, res2.g2.val, 
     rmse.g2.cal, rmse.g2.val, 
     file = here::here("analysis/mbfit_cvg2_none.RData"))

rmse.g2.cal
rmse.g2.val