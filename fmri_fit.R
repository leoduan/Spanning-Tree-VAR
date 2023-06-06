rm(list=ls())

setwd("~/git/treed_var/")
source("treed_var_class.R")

load(file = "HCP_desikan_XY.RDa")

Y<- XY_fmri$Y_concat
X<- XY_fmri$X_concat


model<- TreedVAR$new()

d<- 5
m<- 10
p_star<- 5
lambda<- 0.5


model$init(Y,X,m,p_star)

dim(model$X)
dim(model$L_I_minus_X2)

model$runMCMC(1000)


model$sampleB() #slow
model$sampleTrees()
model$updateA_T()
model$sampleCbar()
model$sampleR()
model$sampleEta()
model$eta_r <- kronecker(model$r, model$eta)
model$sampleSigma2()#slow

model$sampleZ()#slow
model$sampleW()#slow


model$Z %*% t(model$W)
