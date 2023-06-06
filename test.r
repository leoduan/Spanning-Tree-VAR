setwd("~/git/treed_var/")
rm(list=ls())
source("treed_var_class.R")

T<- 1200
d<- 5
p<- 100

m<- 15
p_star<- 5
loglambda<- -10



Cbar0<- matrix(rnorm(d*p*p),p*d)

Cbar0[abs(Cbar0)< quantile(abs(Cbar0),0.99)]<- 0


X<- matrix(rnorm( (T-d)*p*d),(T-d))

XCbar <- X%*%Cbar0


Y<- XCbar + rnorm(length(XCbar))


model<- TreedVAR$new()
  
model$init(Y,X,m,p_star, loglambda_ =  - 50 )

model$runMCMC(100)

sum(model$A_T)


((p-1)*2)*(m-1) +p


for(k in 2:m){
  cover_by_other_trees <- Reduce("+", model$tree_list[1:(k-1)]) > 0
  print(sum((!cover_by_other_trees) & (model$tree_list[[k]]>0)))
}


# for(k in 1:m){
  # cover_by_other_trees <- Reduce("+", model$tree_list[-k]) > 0
  # print(sum((!cover_by_other_trees) & (model$tree_list[[k]]>0)))
# }


image(model$tree_list[[10]])

mean(model$Y_minus_XCbar^2)


# max(abs(model$Cbar))

# plot(abs(model$Cbar[,1]))
