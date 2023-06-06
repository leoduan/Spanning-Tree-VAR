lassoVAR <- setRefClass("lassoVAR",
  fields = list(
    Y = "matrix",
    X = "matrix",
    Cbar = "matrix",
    X2 = "matrix",
    XY = "matrix",
    Y2 = "matrix",
    dp = "integer",
    p = "integer",
    d = "integer",
    T = "integer",
    n = "integer"
  ),
  methods = list(
    elastic_net_loss = function(C,lam1=0.1,lam2=0.0){
        part1 = (sum(diag(Y2)) - 2 * sum((C)*XY) + sum((C%*%t(C))*X2))/2
        # part1 = sum((Y-X%*%C)**2)/2
        part1 + sum(abs(C))*lam1 + sum(C**2)/2*lam2
    },
    softThreshold = function(C, lam){
        absC = abs(C)
        sign(C)*(absC>lam)*(absC-lam)
    },   
    fit_lasso = function(Y_, X_, lam1_rate=0.1 , lam2_rate=0.1, num_iter = 1E5, lr = 1E-8, silent=FALSE) {
        
      Y <<- Y_
      X <<- X_
      p <<- ncol(Y)
      d <<- as.integer(ncol(X) / p)
      T <<- nrow(Y) + d
      
      n <<- dim(X)[2]
        
      lam1 <- lam1_rate* T*p
      lam2 <- lam2_rate* T*p
        
      X2 <<- t(X)%*%X
      XY <<- t(X)%*%Y
      Y2 <<- t(Y)%*%Y 

      dp <<- d * p

     Cbar <<- solve(X2 + diag(1, dp), XY)
        
     for(epoch in 1:num_iter){
        Cbar<<- softThreshold(Cbar+ lr* (XY-X2%*%Cbar -lam2*Cbar), lam1*lr)

        if(epoch %% 100 ==0){
            loss = elastic_net_loss(Cbar,lam1=lam1,lam2=lam2)
            if(!silent){print(loss)}
        }
     }
    }
  )
)
