library(minfi)
library(glmnet)
library(doParallel)

GRFunnorm <- readRDS("rlt/all/QC/GRFunnormNoSNP.rds")
idxRmAllFunnorm <- readRDS("rlt/all/QC/idxRmAllFunnorm.rds")
#phenoRemain <- readRDS("rlt/all/QC/phenoRemain220621_V2.rds")
#phenoNewRemain <- readRDS("rlt/all/QC/phenoNewRemain.rds")
phenoRemain <- readRDS("data/phenotype/phenoRemain221011.rds")

#  801605    338

Beta <- getBeta(GRFunnorm)[!idxRmAllFunnorm,]
sum(colnames(Beta) == phenoRemain$ship_ID) == nrow(phenoRemain)

rm(GRFunnorm)

## there are two samples with missing GA
## 336 samples left
idxNA <- is.na(phenoRemain$gest_age_weeks_birth_r)

tBeta <- t(Beta[,!idxNA])
GA <- phenoRemain$gest_age_weeks_birth_r[!idxNA]

nSample <- length(GA)


# keep only CpGs
probeType <- stringr::str_sub(colnames(tBeta), 1, 2)

flagCg <- probeType == "cg"

tBeta <- tBeta[,flagCg]

load("data/NEST/beta.rds")

flag <- colnames(tBeta) %in% rownames(beta_nest)
tBeta <- tBeta[,flag]

t_beta_nest <- t(beta_nest)
idx <- match(colnames(tBeta), colnames(t_beta_nest))
t_beta_nest <- t_beta_nest[,idx]

sum(colnames(tBeta) == colnames(t_beta_nest)) == ncol(tBeta)


set.seed(777)

idx <- list()
for(i in 1:10){
  if(i < 10) {
    idx[[i]] <- (34*(i-1) + 1) : (34 * i)
  } else {
    idx[[i]] <- 307 : 336
  }
}

GA.min <- NULL
GA.1se <- NULL
for(i in 1:10){
  print(i)
  gaModel <- cv.glmnet(x = tBeta[-idx[[i]],], y = GA[-idx[[i]]], alpha = 0.5,  parallel=TRUE)
  GA.min <- c(GA.min, predict(gaModel, newx = tBeta[idx[[i]],], s = "lambda.min"))
  GA.1se <- c(GA.1se, predict(gaModel, newx = tBeta[idx[[i]],], s = "lambda.1se"))
}

cor.test(GA.min, GA)
cor.test(GA.1se, GA)

save(GA.min, GA.1se, file = "rlt/nestCV10.RData")


set.seed(123)
GA.min <- NULL
GA.1se <- NULL
for(i in 1:nSample){
  print(i)
  gaModel <- cv.glmnet(x = tBeta[-i,], y = GA[-i], alpha = 0.5,  parallel=FALSE)
  GA.min <- c(GA.min, predict(gaModel, newx = tBeta[i,], s = "lambda.min"))
  GA.1se <- c(GA.1se, predict(gaModel, newx = tBeta[i,], s = "lambda.1se"))
}

cor.test(GA.min, GA)
cor.test(GA.1se, GA)

save(GA.min, GA.1se, file = "rlt/nestLOOCV_V2.RData")
