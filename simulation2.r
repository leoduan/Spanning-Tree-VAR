setwd("~/git/treed_var/")
# rm(list=ls())
source("treed_var_class.R")

args <- commandArgs(trailingOnly = TRUE)
T <- as.numeric(args[1])
p <- as.numeric(args[2])


# T<- 1200
d<- 10
# p<- 100

m<- 10
p_star<- 5
loglambda<- -10



sparsify_by_graph<- function(C_bar, A){
    
      C_bar_star<- C_bar
      idx1 <- 1
      idx2 <- p
      for (k in 1:d) {
        C_bar_star[idx1:idx2, ] <- C_bar[idx1:idx2, ]* A
        idx1 <- idx1 + p
        idx2 <- idx2 + p
      }
      C_bar_star
    
}

#generate a graph with tree-rank =2
findMST = function(S) {
  G <- graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE, diag = FALSE)
  G_mst <- mst(graph = G)
  A_T_ <- as.matrix(get.adjacency(G_mst))
  A_T_
}

A0_S<- matrix(0,p,p)

for(l in 1:2){
    S<- matrix(runif(p*p),p)
    S<- S+t(S)
    A0_S <-  A0_S + findMST(S)
    }

A0 = (A0_S>0)
#

diag(A0)<- 1

d0 = 3

Cbar0_0<- matrix(rnorm(d0*p*p),p*d0)
zeros<- matrix(0, (d-d0)*p, p )

Cbar0 <- sparsify_by_graph( rbind(Cbar0_0,zeros), A0)

Sigma2_0<- 0.5**abs(outer(c(1:p),c(1:p),"-"))

rho = sqrt(sum((Cbar0)**2)/sum((Sigma2_0)**2))/0.1

Sigma2_0 = Sigma2_0 * rho

L_Sigma2_0<- chol(Sigma2_0)



X<- matrix(rnorm( (T-d)*p*d),(T-d))

XCbar <- X%*%Cbar0

Y<- XCbar + matrix(rnorm(length(XCbar)),ncol=p)%*%L_Sigma2_0





model_full<- TreedVAR$new()
  
model_full$init(Y,X,m,p_star, loglambda_ =  - 50 )

model_full$runMCMC(1000,burnin=800, silent=TRUE)

model_no_tree<- TreedVAR$new()
  
model_no_tree$init(Y,X,m,p_star, loglambda_ =  - 50 )

model_no_tree$runMCMC(1000,burnin=800, useTrees=FALSE, silent=TRUE)

model_no_gdp<- TreedVAR$new()
  
model_no_gdp$init(Y,X,m,p_star, loglambda_ =  - 50 )

model_no_gdp$runMCMC(1000,burnin=800, useGDP=FALSE, silent=TRUE)



source("var_lasso.R")

lasso<- lassoVAR$new()
lasso$fit_lasso(Y,X,lam1_rate=1E-3, lam2_rate=0, num_iter=2000,lr=1E-4, silent=TRUE)

enets<- lassoVAR$new()
enets$fit_lasso(Y,X,lam1_rate=5E-4, lam2_rate=5E-4, num_iter=2000,lr=1E-4, silent=TRUE)

# dim(model$trace_Cbar)

# relative error
# sum ((enets$Cbar- Cbar0)**2)/ sum(Cbar0**2)

# relative error
# sum ((lasso$Cbar- Cbar0)**2)/ sum(Cbar0**2)

res<- list( model_full=model_full, model_no_tree = model_no_tree, 
           model_no_gdp = model_no_gdp,
           lasso = lasso,
           enets = enets,
           Cbar0=Cbar0, A0 = A0)

filename = paste("sim_results/sim_lowtreerank_T",T,"p",p, sep="_") 

save(res,file=filename)
