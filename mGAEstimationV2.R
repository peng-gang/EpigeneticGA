library(minfi)
library(glmnet)
library(wateRmelon)

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


set.seed(123)
gaModel <- cv.glmnet(x = tBeta, y = GA, alpha = 0.5)

tmp_coeffs <- coef(gaModel, s = "lambda.1se")
coefGA.1se <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)


tmp_coeffs <- coef(gaModel, s = "lambda.min")
coefGA.min <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

GA.min <- predict(gaModel, newx = tBeta[], s = "lambda.min")
GA.1se <- predict(gaModel, newx = tBeta[], s = "lambda.1se")

cor.test(GA, GA.min)
cor.test(GA, GA.1se)


GA.min.nest <- predict(gaModel, newx = t_beta_nest, s = "lambda.min")
GA.1se.nest <- predict(gaModel, newx = t_beta_nest, s = "lambda.1se")

cor.test(pheno_nest$GEST_AGE_WKS, GA.min.nest)
cor.test(pheno_nest$GEST_AGE_WKS, GA.1se.nest)


write.csv(coefGA.1se, file = "rlt/all/mAge/coef_mGA_1se_V2.csv", row.names = FALSE)
write.csv(coefGA.min, file = "rlt/all/mAge/coef_mGA_min_V2.csv", row.names = FALSE)

saveRDS(gaModel, file = "rlt/all/mAge/gaModelV2.rds")






## Other Clocks ----
### Horvath clock ----
ageHorvaths <- wateRmelon::agep(Beta)
ageHorvaths_nest <- wateRmelon::agep(beta_nest)

### PedBE ----
## 94 total, 94 found
anti.trafo= function(x,adult.age=20) {
  ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) 
}

datClock=read.csv("data/PedBE/datcoefInteresting94.csv") 

selectCpGsClock=is.element(rownames(Beta),
                           as.character(datClock[,1][-1]))
sum(selectCpGsClock)
datMethClock0=data.frame(t(Beta[selectCpGsClock,]))

datMethClock= data.frame(datMethClock0[as.character(datClock[,1][-1])])
sum(colnames(datMethClock) == as.character(datClock[,1][-1])) == ncol(datMethClock)

PedBE_age <- as.numeric(
  anti.trafo(
    datClock[1,2]+as.numeric(as.matrix(datMethClock)%*%as.numeric(datClock[,2][-1]))
  )
)


# 94
selectCpGsClock=is.element(rownames(beta_nest),
                           as.character(datClock[,1][-1]))
sum(selectCpGsClock)
datMethClock0=data.frame(t(beta_nest[selectCpGsClock,]))

datMethClock= data.frame(datMethClock0[as.character(datClock[,1][-1])])
sum(colnames(datMethClock) == as.character(datClock[,1][-1])) == ncol(datMethClock)

PedBE_age_nest <- as.numeric(
  anti.trafo(
    datClock[1,2]+as.numeric(as.matrix(datMethClock)%*%as.numeric(datClock[,2][-1]))
  )
)




### EPIC GA ----
coefEPICGA <- read.csv("data/mAgeGA/coefEPICGA.csv", header = TRUE)
# 176 total, 170 CpG found
idx <- match(coefEPICGA$cpgs[-1], rownames(Beta))
sum(!is.na(idx))
sum(rownames(Beta)[idx] == coefEPICGA$cpgs[-1], na.rm = TRUE) == sum(!is.na(idx))
betaSel <- Beta[idx,]
rownames(betaSel) == coefEPICGA$cpgs[-1]
betaSel <- as.matrix(betaSel)
betaSel[is.na(betaSel)] <- 0
EPICGA <- (t(betaSel) %*% coefEPICGA$s0[-1] + coefEPICGA$s0[1])/7

# 176 total, 87 CpG found
idx <- match(coefEPICGA$cpgs[-1], rownames(beta_nest))
sum(!is.na(idx))
sum(rownames(beta_nest)[idx] == coefEPICGA$cpgs[-1], na.rm = TRUE) == sum(!is.na(idx))
betaSel <- beta_nest[idx,]
rownames(betaSel) == coefEPICGA$cpgs[-1]
betaSel <- as.matrix(betaSel)
betaSel[is.na(betaSel)] <- 0
EPICGA_nest <- (t(betaSel) %*% coefEPICGA$s0[-1] + coefEPICGA$s0[1])/7



### Knight Clock ----
coefKnight<- read.csv("data/mAgeGA/coefKnight.csv", header = TRUE)
# 148 total, 141 found
idx <- match(coefKnight$CpGmarker[-1], rownames(Beta))
sum(!is.na(idx))
sum(rownames(Beta)[idx] == coefKnight$CpGmarker[-1], na.rm = TRUE) == sum(!is.na(idx))
betaSel <- Beta[idx,]
rownames(betaSel) == coefKnight$CpGmarker[-1]
betaSel <- as.matrix(betaSel)
betaSel[is.na(betaSel)] <- 0
KnightGA <- t(betaSel) %*% coefKnight$CoefficientTraining[-1] + coefKnight$CoefficientTraining[1]


idx <- match(coefKnight$CpGmarker[-1], rownames(beta_nest))
sum(!is.na(idx))
sum(rownames(beta_nest)[idx] == coefKnight$CpGmarker[-1], na.rm = TRUE) == sum(!is.na(idx))
betaSel <- beta_nest[idx,]
rownames(betaSel) == coefKnight$CpGmarker[-1]
betaSel <- as.matrix(betaSel)
betaSel[is.na(betaSel)] <- 0
KnightGA_nest <- t(betaSel) %*% coefKnight$CoefficientTraining[-1] + coefKnight$CoefficientTraining[1]



### Bohlin Clock ----
# https://github.com/JonBohlin/predictGA
library(predictGA)

tmp_coeffs <- coef(UL.mod.cv, s = "lambda.1se")
coefBohlin <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)


cpgs<-extractSites(type="se")
# allcpgs<-extractSites(type="all")
# numsamples<-10
# mlmatr<-matrix(NA, ncol=length(allcpgs), nrow=numsamples)
# mlmatr<-data.frame(mlmatr)
# 
# for(i in cpgs)
#   mlmatr[,i]<-runif(numsamples, min=0, max=1)
# 
# mypred<-predictGA(mlmatr)
# head(mypred)


## The CpGs needed for prediction can not contain NA's


mlmatr<-matrix(0, ncol=1000, nrow=ncol(Beta))
mlmatr<-data.frame(mlmatr)

betaSel <- matrix(0, ncol = length(cpgs), nrow = ncol(Beta))
colnames(betaSel) <- cpgs

## 85 out of 96
idx <- match(cpgs, rownames(Beta))
sum(!is.na(idx))
sum(!is.na(idx))
tmp <- Beta[idx,]
rownames(tmp) == cpgs
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp[is.na(tmp)] <- 0
colnames(tmp) == cpgs
colnames(tmp) <- cpgs

rownames(mlmatr) <- rownames(tmp)
mlmatr <- cbind(mlmatr, tmp)

BohlinGA <- predictGA(tmp, transp = FALSE)




mlmatr<-matrix(0, ncol=1000, nrow=ncol(beta_nest))
mlmatr<-data.frame(mlmatr)

betaSel <- matrix(0, ncol = length(cpgs), nrow = ncol(beta_nest))
colnames(betaSel) <- cpgs

## 96 out of 96
idx <- match(cpgs, rownames(beta_nest))
sum(!is.na(idx))
sum(!is.na(idx))
tmp <- beta_nest[idx,]
rownames(tmp) == cpgs
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp[is.na(tmp)] <- 0
colnames(tmp) == cpgs
colnames(tmp) <- cpgs

rownames(mlmatr) <- rownames(tmp)
mlmatr <- cbind(mlmatr, tmp)

BohlinGA_nest <- predictGA(tmp, transp = FALSE)



idxNA <- is.na(phenoRemain$gest_age_weeks_birth_r)

mGA <- data.frame(
  ShipID = phenoRemain$ship_ID[!idxNA],
  GA = phenoRemain$gest_age_weeks_birth_r[!idxNA],
  HorvathAge = ageHorvaths$horvath.age[!idxNA],
  PedBEAge = PedBE_age[!idxNA],
  KnightAge = KnightGA[!idxNA],
  BohlinAge = BohlinGA$GA[!idxNA]/7,
  EPICAge = EPICGA[!idxNA],
  GAMin = GA.min,
  GA1se = GA.1se,
  stringsAsFactors = FALSE
)

mGA_nest <- data.frame(
  id = pheno_nest$nest_id,
  GA = pheno_nest$GEST_AGE_WKS,
  HorvathAge = ageHorvaths_nest$horvath.age,
  PedBEAge = PedBE_age_nest,
  KnightAge = KnightGA_nest,
  BohlinAge = BohlinGA_nest$GA/7,
  EPICAge = EPICGA_nest,
  GAMin = GA.min.nest,
  GA1se = GA.1se.nest,
  stringsAsFactors = FALSE
)


save(mGA, mGA_nest, file = "rlt/all/mAge/mGAV2.RData")


