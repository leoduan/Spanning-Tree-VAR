suppressMessages(require("igraph"))

TreedVAR <- setRefClass("TreedVAR",
  fields = list(
    Y = "matrix",
    X = "matrix",
    T = "integer",
    d = "integer",
    p = "integer",
    m = "integer",
    p_star = "integer",
    A_T = "matrix",
    dp = "integer",
    loglambda = "numeric",
    Cbar = "matrix",
    Y_minus_XCbar = "matrix",
    W = "matrix",
    Z = "matrix",
    B = "matrix",
    ZW = "matrix",
    eta_r = "matrix",
    eta = "matrix",
    r = "numeric",
    L_I_minus_X2 = "matrix",
    I_minus_X2 = "matrix",
    sigma2 = "numeric",
    m_tilde = "numeric",
    tree_list = "list",
    trace_Cbar = "array",
    trace_A_T = "array",
    trace_eta = "array"
  ),
  methods = list(
    rgumbel = function(n) {
      -log(-log(runif(n)))
    },
    findMST = function(logS) {
      G <- graph_from_adjacency_matrix(-logS, mode = "undirected", weighted = TRUE, diag = FALSE)
      G_mst <- mst(graph = G)
      A_T_ <- as.matrix(get.adjacency(G_mst))
      A_T_
    },
    drawT_fast = function(logS) {
      n <- nrow(logS)
      gumbelMat <- matrix(0, n, n)
      gumbelMat[lower.tri(gumbelMat, diag = FALSE)] <- rgumbel((n) / 2 * (n - 1))
      gumbelMat <- gumbelMat + t(gumbelMat)
      A_T_ <- findMST(logS + gumbelMat)
      A_T_
    },
    init = function(Y_, X_, m_ = 10, p_star_ = 5, loglambda_ = -10) {
      Y <<- Y_
      X <<- X_
      p <<- ncol(Y)
      d <<- as.integer(ncol(X) / p)
      T <<- nrow(Y) + d

      m <<- as.integer(m_)
      p_star <<- as.integer(p_star_)

      W <<- matrix(rnorm(p * p_star), p)
      Z <<- matrix(0, T - d, p_star)
      ZW <<- Z %*% t(W)

      loglambda <<- loglambda_
      # initialize m_tilde
      X2 <- t(X) %*% X

      X2_spec_norm <- max(eigen(X2)$values)
      m_tilde <<- X2_spec_norm * (1.0+1E-3)
      dp <<- d * p

      I_minus_X2 <<- diag(1, dp) * m_tilde - X2
      L_I_minus_X2 <<- t(chol(I_minus_X2))


      tree_list <<- lapply(c(1:m), function(x) drawT_fast(matrix(0, p, p)))
      A_T <<- matrix(1, p, p)



      Cbar <<- solve(X2 + diag(1, dp), t(X) %*% Y)
      sigma2 <<- 0.9
      Y_minus_XCbar<<- Y- X%*%Cbar

      u0 <- matrix(rnorm(d * p * p), d * p)
      U <- L_I_minus_X2 %*% u0 * sqrt(sigma2) + (I_minus_X2) %*% Cbar


      YZW <- Y - ZW
      B <<- t(X) %*% (YZW) + U

      eta <<- matrix(10, p, p)
      r <<- 0.1^c(1:d) # rep(1,d)

      eta_r <<- kronecker(r, eta)
    },
    sumOverLags = function(X) {
      s <- matrix(0, p, p)
      idx1 <- 1
      idx2 <- p
      for (k in 1:d) {
        s <- s + X[idx1:idx2, ]
        idx1 <- idx1 + p
        idx2 <- idx2 + p
      }
      s
    },
    sumOverGraph = function(X) {
      s <- rep(0, d)
      idx1 <- 1
      idx2 <- p
      for (k in 1:d) {
        s[k] <- sum(X[idx1:idx2, ])
        idx1 <- idx1 + p
        idx2 <- idx2 + p
      }
      s
    },
    sampleB = function() {
      u0 <- matrix(rnorm(d * p * p), d * p)
      U <- L_I_minus_X2 %*% u0 * sqrt(sigma2) + (I_minus_X2) %*% Cbar
      YZW <- Y - ZW
      B <<- t(X) %*% (YZW) + U
    },
    sampleCbar = function() {
      var_cbar <- 1 / (m_tilde + 1 / eta_r)
      m_cbar <- B * var_cbar
      A_T_rep <- (kronecker(matrix(1, d), A_T))
      Cbar <<- (m_cbar + sqrt(var_cbar) * sqrt(sigma2) * rnorm(d * p * p)) * A_T_rep
      Y_minus_XCbar<<- Y- X%*%Cbar
    },
    sampleTrees = function() {
      var_cbar <- 1 / (m_tilde + 1 / eta_r)
      part1 <- log(var_cbar) / 2 - log(eta_r) / 2
      part2 <- B^2 * var_cbar / 2 / sigma2

      logW0 <- sumOverLags(part1 + part2)

      # add the upper and lower triangular parts (note: it's ok that the diagonals are added twice, because they will not be used for sampling trees.)

      logW <- logW0 + t(logW0)

      logWLambda <- logW + loglambda



      for (l in 1:m) {
        cover_by_other_trees <- Reduce("+", tree_list[-l]) > 0
        tree_list[[l]] <<- drawT_fast((1 - cover_by_other_trees) * logWLambda)
      }
    },
    updateA_T = function() {
      A_T <<- (Reduce("+", tree_list) > 0) * 1
      diag(A_T) <<- 1
    },
    sampleR = function() {
      a_k <- p*p
      b_k <- 0.1^c(1:d) * p*p

      par1 <- sum(A_T) / 2 + a_k
      par2 <- sumOverGraph(Cbar^2 / (kronecker(matrix(1, d, 1), eta)) / sigma2 / 2) + b_k
      r <<- 1 / rgamma(d, par1, rate = par2)
    },
    sampleEta = function() {
      alpha_eta <- 3
      gamma_eta <- 1E-3

      par1 <- alpha_eta + 1
      par2 <- 1 / (1 / gamma_eta + 1 / eta)
      beta <- matrix(rgamma(p * p, par1, 1), p) * par2

      par1 <- d * A_T / 2 + alpha_eta
      par2 <- sumOverLags(Cbar^2 / (kronecker(r, matrix(1, p, p))) / sigma2 / 2) + beta
      eta <<- matrix(1 / rgamma(p * p, c(par1)), p, p) * par2
    },
    sampleZ = function() {
      W2_I <- t(W) %*% W / sigma2 + diag(1, p_star)
      W2_I_inv <- solve(W2_I)
      par2 <- W2_I_inv
      par1 <- Y_minus_XCbar %*% (W %*% W2_I_inv / sigma2)

      z0 <- matrix(rnorm((T - d) * p_star), T - d)
      Z <<- z0 %*% chol(par2) + par1
    },
    sampleW = function() {
      gamma_W <- 1E-3
      Z2_I <- t(Z) %*% Z / sigma2 + diag(1 / gamma_W, p_star)
      Z2_I_inv <- solve(Z2_I)
      par2 <- Z2_I_inv
      par1 <- t(Y_minus_XCbar) %*% Z %*% Z2_I_inv / sigma2

      w0 <- matrix(rnorm(p * p_star), p)
      W <<- w0 %*% chol(par2) + par1
    },
    sampleSigma2 = function() {
      alpha_sigma <- 2
      beta_sigma <- 0.1

      par1 <- ((T - d) * p + sum(A_T) * d) / 2 + alpha_sigma
      par2 <- sum((Y_minus_XCbar - ZW)^2) / 2 + sum(Cbar^2 / eta_r) / 2 + beta_sigma
      sigma2 <<- 1 / rgamma(1, par1, rate = par2)
    },
    runMCMC = function(num_iter = 100, burnin=50, useTrees=TRUE, useGDP=TRUE, silent=FALSE) {
        
    if(!silent){
      pb <- txtProgressBar(
        min = 0, # Minimum value of the progress bar
        max = num_iter, # Maximum value of the progress bar
        style = 3, # Progress bar style (also available style = 1 and style = 2)
        width = 50, # Progress bar width. Defaults to getOption("width")
        char = "="
      ) # Character used to create the bar
        }
      trace_Cbar <<- array(0, dim = c(dp, p, num_iter-burnin))
      trace_A_T <<- array(0, dim = c(p, p, num_iter-burnin))
      trace_eta <<- array(0, dim = c(p, p, num_iter-burnin))
        
      if(!useGDP){
          eta <<- matrix(1,p,p)
      }
        
      if(!useTrees){
          A_T <<- matrix(1,p,p)
      }
        

      for (step in 1:num_iter) {
        sampleB()
          
          
        if(useTrees){
            sampleTrees()
            updateA_T()
        }
        sampleCbar()
        sampleR()
          
        if(useGDP){
            sampleEta()
        }
          
        eta_r <<- kronecker(r, eta)
        sampleSigma2()

        sampleZ()
        sampleW()
        ZW <<- Z %*% t(W)

        if(step>burnin){
            trace_Cbar[, , step - burnin] <<- Cbar
            trace_A_T[, , step - burnin] <<- A_T
            trace_eta[, , step - burnin] <<- eta
        }
          
        if(!silent){
         setTxtProgressBar(pb, step)
        }
        # print(sum(model$A_T))
      }
    }
  )
)