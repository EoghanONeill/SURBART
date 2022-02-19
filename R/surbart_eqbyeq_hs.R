# get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2)
  }

  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))

  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)

  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}




#' @title Seemingly Unrelated Regression Bayesian Additive Regression Trees implemented using MCMC with equation-by-equation tree draws and Horseshoe prior on covariances.
#'
#' @description Seemingly Unrelated Regression Bayesian Additive Regression Trees implemented using MCMC with equation-by-equation tree drawsand Horseshoe prior on covariances. Implementation is edited from Huber et al (2021).
#' @import dbarts
#' @import truncnorm
#' @import LaplacesDemon
#' @import MASS
#' @import stochvol
#' @import invgamma
#' @param x.train The training covariate data for all training observations. Matrix or list. If one matrix (not a list), then the same set of covariates are used for each outcome. Each element of the list is a covariate matrix that  corresponds to a different outcome variable. Number of rows equal to the number of observations. Number of columns equal to the number of covariates.
#' @param x.test The test covariate data for all test observations. Matrix or list. If one matrix (not a list), then the same set of covariates are used for each outcome. Each element of the list is a covariate matrix that  corresponds to a different outcome variable. Number of rows equal to the number of observations. Number of columns equal to the number of covariates.
#' @param y The training data list of vectors of outcomes. The length of the list should equal the number of types of outcomes. Each element of the list should be a vector of length equal to the number of units.
#' @param num_outcomes The number of outcome variables.
#' @param num_obs The number of observations per outcome.
#' @param num_test_obs THe number of test observations per outcome.
#' @param n.iter Number of iterations excluding burnin.
#' @param n.burnin Number of burnin iterations.
#' @param n.trees (dbarts control option) A positive integer giving the number of trees used in the sum-of-trees formulation. Assuming each outcome variable is modelled using the same number of trees. Each outcome is modelled using a distinct sum of trees (as opposed to shared trees).
#' @param n.chains (dbarts control option) A positive integer detailing the number of independent chains for the dbarts sampler to use (more than one chain is unlikely to improve speed because only one sample for each call to dbarts).
#' @param n.threads  (dbarts control option) A positive integer controlling how many threads will be used for various internal calculations, as well as the number of chains. Internal calculations are highly optimized so that single-threaded performance tends to be superior unless the number of observations is very large (>10k), so that it is often not necessary to have the number of threads exceed the number of chains.
#' @param printEvery (dbarts control option)If verbose is TRUE, every printEvery potential samples (after thinning) will issue a verbal statement. Must be a positive integer.
#' @param printCutoffs (dbarts control option) A non-negative integer specifying how many of the decision rules for a variable are printed in verbose mode
#' @param rngKind (dbarts control option) Random number generator kind, as used in set.seed. For type "default", the built-in generator will be used if possible. Otherwise, will attempt to match the built-in generator’s type. Success depends on the number of threads.
#' @param rngNormalKind (dbarts control option) Random number generator normal kind, as used in set.seed. For type "default", the built-in generator will be used if possible. Otherwise, will attempt to match the built-in generator’s type. Success depends on the number of threads and the rngKind
#' @param rngSeed (dbarts control option) Random number generator seed, as used in set.seed. If the sampler is running single-threaded or has one chain, the behavior will be as any other sequential algorithm. If the sampler is multithreaded, the seed will be used to create an additional pRNG object, which in turn will be used sequentially seed the threadspecific pRNGs. If equal to NA, the clock will be used to seed pRNGs when applicable.
#' @param updateState (dbarts control option) Logical setting the default behavior for many sampler methods with regards to the immediate updating of the cached state of the object. A current, cached state is only useful when saving/loading the sampler.
#' @param tree.prior (dbarts option) An expression of the form dbarts:::cgm or dbarts:::cgm(power,base) setting the tree prior used in fitting.
#' @param node.prior (dbarts option) An expression of the form dbarts:::normal or dbarts:::normal(k) that sets the prior used on the averages within nodes.
#' @param resid.prior (dbarts option) An expression of the form dbarts:::chisq or dbarts:::chisq(df,quant) that sets the prior used on the residual/error variance
#' @param proposal.probs (dbarts option) Named numeric vector or NULL, optionally specifying the proposal rules and their probabilities. Elements should be "birth_death", "change", and "swap" to control tree change proposals, and "birth" to give the relative frequency of birth/death in the "birth_death" step.
#' @param sigmadbarts (dbarts option) A positive numeric estimate of the residual standard deviation. If NA, a linear model is used with all of the predictors to obtain one.
#' @param print.opt Print every print.opt number of Gibbs samples.
#' @param quiet Does not show progress bar if TRUE
#' @param sv If TRUE use stochastic voloatility prior
#' @param outcome_draws If TRUE, output draws of the outcome.
#' @export
#' @return The following objects are returned:
#' \item{mutrain_draws}{Matrix of MCMC draws of expected values for training observations. Number of rows equal to the number of training observations multiplied by the number of outcomes. The rows are ordered by beginning with all N (units') observations for the first outcome variable, then all N for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{mutest_draws}{Matrix of MCMC draws of expected values for test observations. Number of rows equal to the number of test observations multiplied by the number of outcomes. The rows are ordered by beginning with all Ntest (units') observations for the first outcome variable, then all Ntest for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{ytrain_draws}{Matrix of MCMC draws of outcomes for training observations. Number of rows equal to the number of observations multiplied by the number of outcomes. The rows should be ordered by beginning with all N (units') observations for the first outcome variable, then all N for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{ytest_draws}{Matrix of MCMC draws of outcomes for training observations. Number of rows equal to the number of observations multiplied by the number of outcomes. The rows should be ordered by beginning with all N (units') observations for the first outcome variable, then all N for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{Sigma_draws}{3 dimensional array of MCMC draws of the covariance matrix for the outcome-specific error terms. The numbers of rows and columns equal are equal to the number of outcome variables. The number of slices is }
#' @examples
#'
#' library(foreign)
#' library(systemfit)
#'
#' hsb2 <- read.dta("https://stats.idre.ucla.edu/stat/stata/notes/hsb2.dta")
#'
#' train_inds <- sample(1:nrow(hsb2),size = 180, replace = FALSE)
#' test_inds <- (1:nrow(hsb2))[-train_inds]
#'
#' hsb2train <- hsb2[train_inds,]
#' hsb2test <- hsb2[test_inds,]
#'
#'
#' r1 <- read~female + as.numeric(ses) + socst
#' r2 <- math~female + as.numeric(ses) + science
#'
#' fitsur <- systemfit(list(readreg = r1, mathreg = r2), data=hsb2train)
#'
#' linSURpreds <- predict(fitsur, hsb2test)
#'
#' sqrt(mean((linSURpreds$readreg.pred-hsb2test$read)^2))
#'
#' sqrt(mean((linSURpreds$mathreg.pred-hsb2test$math)^2))
#'
#'
#'
#' xtrain <- list()
#' xtrain[[1]] <- cbind(hsb2train$female, as.numeric(hsb2train$ses), hsb2train$socst)
#' xtrain[[2]] <- cbind(hsb2train$female, as.numeric(hsb2train$ses), hsb2train$science)
#'
#' xtest <- list()
#' xtest[[1]] <- cbind(hsb2test$female, as.numeric(hsb2test$ses), hsb2test$socst)
#' xtest[[2]] <- cbind(hsb2test$female, as.numeric(hsb2test$ses), hsb2test$science)
#'
#' ylist <- list()
#' ylist[[1]] <- hsb2train$read
#' ylist[[2]] <- hsb2train$math
#'
#'
#'
#' surbartres <- surbart_eqbyeq_hs(x.train = xtrain, #either one matrix or list
#'                              x.test = xtest, #either one matrix or list
#'                              y = ylist,
#'                              2,
#'                              nrow(hsb2train),
#'                              nrow(hsb2test),
#'                              n.iter = 1000,
#'                              n.burnin = 100)
#'
#'
#' readpredsurbart <- apply(surbartres$mutest_draws[,,1],2,mean)
#'
#' mathpredsurbart <- apply(surbartres$mutest_draws[,,2],2,mean)
#'
#' sqrt(mean((readpredsurbart  -hsb2test$read)^2))
#'
#' sqrt(mean((mathpredsurbart- hsb2test$math)^2))
#'
#'
#'
#'
#' @export
surbart_eqbyeq_hs <- function(x.train, #either one matrix or list
                           x.test, #either one matrix or list
                           y,
                           num_outcomes,
                           num_obs,
                           num_test_obs,
                           n.iter=1000,
                           n.burnin=100,
                           n.trees = 50L,
                           n.burn = 0L,
                           n.samples = 1L,
                           n.thin = 1L,
                           n.chains = 1,
                           n.threads = guessNumCores(),
                           printEvery = 100L,
                           printCutoffs = 0L,
                           rngKind = "default",
                           rngNormalKind = "default",
                           rngSeed = NA_integer_,
                           updateState = FALSE,
                           tree.prior = dbarts:::cgm,
                           node.prior = dbarts:::normal,
                           resid.prior = dbarts:::chisq,
                           proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
                           sigmadbarts = NA_real_,
                           print.opt = 100,
                           quiet = FALSE,
                           sv=FALSE,
                           outcome_draws = FALSE){




  # if(is.vector(x.train) | is.factor(x.train)| is.data.frame(x.train)) x.train = as.matrix(x.train)
  # if(is.vector(x.test) | is.factor(x.test)| is.data.frame(x.test)) x.test = as.matrix(x.test)

  if((!is.matrix(x.train))&(!is.list(x.train))) stop("argument x.train must be a double matrix or a list of matrices.")
  if((!is.matrix(x.test))&(!is.list(x.test))) stop("argument x.test must be a double matrix or a list of matrices.")



  if(is.matrix(x.train)){
    if(!(is.matrix(x.test))){stop("Both x.train and x.test must be matrices, or both must be lists.")}
    print("input is a matrix. Using same covariates for all outcomes.")
    Xlist <- lapply(1:num_outcomes, function(j) x.train)
    Xtestlist <- lapply(1:num_outcomes, function(j) x.test)
    if(nrow(x.train) != num_obs) stop("number of rows of x.train must equal num_obs.")
    if(nrow(x.test) != num_test_obs) stop("number of rows of x.test must equal num_test_obs.")

    if(nrow(x.train) != length(y[[1]])) stop("number of rows of x.train must equal length of each element of y.train")
    if((ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")


  }


  if(is.list(x.train)){
    print("input is a list. Using separate covariate matrix for each outcome.")
    if(!(is.matrix(x.train[[1]]))){stop("Elements of the x.train list must be matrices.")}

    if(!(is.list(x.test))){stop("Both x.train and x.test must be matrices, or both must be lists.")}
    if(!(is.matrix(x.test[[1]]))){stop("Elements of the x.test list must be matrices.")}

    if(length(x.train) != num_outcomes ){ stop("Number of elements of covariate matrix list x.train must equal num_outcomes.") }
    Xlist <- x.train
    Xtestlist <- x.test

    if(nrow(x.train[[1]]) != num_obs) stop("number of rows of each element of x.train must equal num_obs.")
    if(nrow(x.test[[1]]) != num_test_obs) stop("number of rows of x.test must equal num_test_obs.")

    if((ncol(x.test[[1]])!=ncol(x.train[[1]]))) stop("Each element of input list (of matrices) x.test must have the same number of columns as each element of x.train list.")


    if(nrow(x.train[[1]]) != length(y[[1]])) stop("number of rows in each element of x.train must equal length of each element of y.train")

  }

  if(length(y) != num_outcomes) stop("length of list y must equal number of outcomes.")
  if(length(y[[1]]) != num_obs) stop("length of each vector element of the list y must equal number of observations.")


  n <- length(y)




  # draw = list(
  #   Z.mat = array(NA, dim = c(n, n.iter)),
  #   Z.matcens = array(NA, dim = c(n0, n.iter)),
  #   #Z.matuncens = array(NA, dim = c(n1, n.iter)),
  #   Z.matcensbelow = array(NA, dim = c(n_censbelow, n.iter)),
  #   Z.matcensabove = array(NA, dim = c(n_censabove, n.iter)),
  #   mu = array(NA, dim = c(n, n.iter)),#,
  #   mucens = array(NA, dim = c(n0, n.iter)),#,
  #   muuncens = array(NA, dim = c(n1, n.iter)),#,
  #   mucensbelow = array(NA, dim = c(n_censbelow, n.iter)),#,
  #   mucensabove = array(NA, dim = c(n_censabove, n.iter)),#,
  #   ystar = array(NA, dim = c(n, n.iter)),#,
  #   ystarcens = array(NA, dim = c(n0, n.iter)),#,
  #   ystaruncens = array(NA, dim = c(n1, n.iter)),#,
  #   ystarcensbelow = array(NA, dim = c(n_censbelow, n.iter)),#,
  #   ystarcensabove = array(NA, dim = c(n_censabove, n.iter)),#,
  #   test.mu =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.y_nocensoring =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.y_withcensoring =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.probcensbelow =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.probcensabove =  array(NA, dim = c(ntest, n.iter)),
  #   sigma = rep(NA, n.iter)
  # )



  Ymu <- rep(NA, num_outcomes)
  Ysd <- rep(NA, num_outcomes)


  # standardize outcomes
  for(i in 1:num_outcomes){
    Ymu[i] <- mean(y[[i]])
    Ysd[i] <- sd(y[[i]])
    y[[i]] <- (y[[i]]-Ymu[i])/Ysd[i]
  }

  ymat <- matrix(NA, nrow = num_obs, ncol = num_outcomes)

  for(i in 1:num_outcomes){
    ymat[,i] <- y[[i]]
  }

  XAOLStemp <- matrix(NA, nrow = num_obs, ncol = num_outcomes)

  for(mm in 1:num_outcomes){
    A_OLS <- solve(crossprod(Xlist[[i]]))%*%crossprod(Xlist[[i]],y[[i]])
    XAOLStemp[,mm] <- Xlist[[i]]%*%A_OLS
  }
  # Sig_OLS <- crossprod(Y-X%*%A_OLS)/T
  Sig_OLS <- crossprod(ymat-XAOLStemp)/num_obs

  # Sig_t <- array(0,dim=c(num_obs,num_outcomes,num_outcomes))
  # for(tt in 1:num_obs){
  #   Sig_t[tt,,] <- Sig_OLS
  # }

  # covariance related objects
  eta <- matrix(NA,num_obs,num_outcomes)
  H <- matrix(-3,num_obs,num_outcomes)
  A0_draw <- diag(num_outcomes)

  # stochastic volatility using stochvol package
  sv_draw <- list()
  sv_latent <- list()
  for (mm in seq_len(num_outcomes)){
    sv_draw[[mm]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent[[mm]] <- rep(0,num_obs)
  }

  # construct priors for SV
  sv_priors <- list()
  if(sv){
    for(mm in 1:num_outcomes){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 10),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 10),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
  }else{
    for(mm in 1:num_outcomes){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_constant(0),
        phi = sv_constant(1-1e-12),
        sigma2 = sv_constant(1e-12),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
  }

  # A_prior <- matrix(0,K,num_outcomes)
  # theta_A <- matrix(1,K,num_outcomes)
  theta_A0 <- matrix(1,num_outcomes,num_outcomes)



  #Set tree sampler parameters

  control <- dbartsControl(updateState = updateState, verbose = FALSE,  keepTrainingFits = TRUE,
                           keepTrees = TRUE,
                           n.trees = n.trees,
                           n.burn = n.burn,
                           n.samples = n.samples,
                           n.thin = n.thin,
                           n.chains = n.chains,
                           n.threads = n.threads,
                           printEvery = printEvery,
                           printCutoffs = printCutoffs,
                           rngKind = rngKind,
                           rngNormalKind = rngNormalKind,
                           rngSeed = rngSeed)


  # print(colnames(Xmat.train))

  # print("begin dbarts")


  #take initial tree draws

  sampler.list <- list()
  svdraw.list <- list()

  if(nrow(Xtestlist[[1]])==0){


    for (jj in 1:num_outcomes){

      Xmat.train <- data.frame(y = y[[jj]], x = Xlist[[jj]] )

      sampler.list[[jj]] <- dbarts(y ~ .,
                                   data = Xmat.train,
                                   #test = x.test,
                                   control = control,
                                   tree.prior = tree.prior,
                                   node.prior = node.prior,
                                   resid.prior = resid.prior,
                                   proposal.probs = proposal.probs,
                                   sigma = sqrt(Sig_OLS[jj,jj]))
    }


  }else{
    for (jj in 1:num_outcomes){

      Xmat.train <- data.frame(y = y[[jj]], x = Xlist[[jj]] )
      Xmat.test <- data.frame(x = Xtestlist[[jj]] )

      sampler.list[[jj]] <- dbarts(y ~ .,
                                   data = Xmat.train,
                                   test = Xmat.test,
                                   control = control,
                                   tree.prior = tree.prior,
                                   node.prior = node.prior,
                                   resid.prior = resid.prior,
                                   proposal.probs = proposal.probs,
                                   sigma = sqrt(Sig_OLS[jj,jj])
      )
    }

  }

  sampler.run <- list()
  sigma.mat <- matrix(NA, num_outcomes, 1)

  # -----------------------------------------------------------------------------
  # Initialize HS prior on covariances (if any)
  lambda.A0 <- 1
  nu.A0 <- 1
  tau.A0 <- 1
  zeta.A0 <- 1
  prior.cov <- rep(1, num_outcomes*(num_outcomes-1)/2)


  # ------------------------------------------------------------------
  preds.train <- matrix(NA, num_obs, num_outcomes)

  preds.test <- matrix(NA, num_test_obs, num_outcomes)

  # sigma.mat <- matrix(NA, num_outcomes, num_outcomes)

  # count.mat <- list() #stored as list to allow for varying numbers of splitting variables per outcome
  # require separate counts for each tree


  Y_store <- array(NA,dim=c(n.iter,num_obs,num_outcomes))
  mu_store <- array(NA,dim=c(n.iter,num_obs,num_outcomes))

  Ytest_store <- array(NA,dim=c(n.iter,num_test_obs,num_outcomes))
  mutest_store <- array(NA,dim=c(n.iter,num_test_obs,num_outcomes))

  Sigma_store <- array(NA,dim=c(num_outcomes,num_outcomes, n.iter))



  # eta <- matrix(NA, nrow = num_obs, ncol =  num_outcomes)

  # start Gibbs sampler
  # show progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = n.iter+n.burnin, style = 3)
    start  <- Sys.time()
  }



  #loop through the Gibbs sampler iterations
  for(iter in 1:(n.iter+n.burnin)){




    for (mm in 1:num_outcomes){
      if (mm > 1){
        Z_mm <- eta[,1:(mm-1), drop=F]
        A0_mm <- A0_draw[mm,1:(mm-1)]
        sampler.list[[mm]]$setResponse(ymat[,mm] - Z_mm%*%A0_mm)
      }

      #Draw all trees for outcome mm in iteration iter
      rep_mm <- sampler.list[[mm]]$run(0L, 1L) # construct BART sample using dbarts (V. Dorie)

      sampler.run[[mm]] <- rep_mm

      #line commented out because sampling from inverse wishart, not the prior from BAVART paper
      sigma.mat[mm,] <- rep_mm$sigma

      if (any(is.na(rep_mm$train))) break
      eta[,mm] <- ymat[,mm] - rep_mm$train[,1]
      # A_draw[,mm] <- X.ginv%*%rep_mm$train
      preds.train[,mm] <-  rep_mm$train[,1]


      # print("rep_mm$test[,1] =")
      # print(rep_mm$test[,1])
      preds.test[,mm] <-  rep_mm$test[,1]

      # count.mat[,mm] <- rep_mm$varcount
      if (mm > 1){
        norm_mm <- as.numeric(exp(-.5*sv_latent[[mm]]) * 1/sigma.mat[mm,])
        u_mm <- eta[,1:(mm-1),drop=F]*norm_mm
        eta_mm <- eta[,mm]*norm_mm
        if (mm == 2) v0.inv <- 1/theta_A0[mm,1] else v0.inv <- diag(1/theta_A0[mm,1:(mm-1)])
        V.cov <- solve(crossprod(u_mm) + v0.inv)
        mu.cov <- V.cov %*% crossprod(u_mm, eta_mm)
        mu.cov.draw <- mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov))
        A0_draw[mm,1:(mm-1)] <- mu.cov.draw
      }


    } # end loop over mm



    shock_norm <- eta %*% t(solve(A0_draw))
    if(sv){
      for (mm in seq_len(num_outcomes)){
        svdraw_mm <- svsample_general_cpp(shock_norm[,mm]/sigma.mat[mm], startpara = sv_draw[[mm]], startlatent = sv_latent[[mm]], priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        sv_latent[[mm]] <- svdraw_mm$latent

        normalizer <- as.numeric(exp(-.5*svdraw_mm$latent))
        weights.new <- as.numeric(exp(-svdraw_mm$latent))

        dat <- dbartsData(formula = ymat[,mm]~Xlist[[mm]],weights=weights.new)
        sampler.list[[mm]]$setData(dat)
        # H[,mm] <- log(sigma.mat[mm]^2) + svdraw_mm$latent
      }
    }#else{
    #   H[,mm] <- log(sigma.mat[mm]^2)
    # }

    # for(tt in seq_len(num_obs)){
    #   S_tmp <- exp(H[tt,])
      # S.t <- t(A0_draw)%*%crossprod(diag(S_tmp),(A0_draw))
      # Sig_t[tt,,] <- S.t
    # }

    # 2) updating shrinkage priors
    # hierarchical prior values
    hs_draw <- get.hs(bdraw=A0_draw[lower.tri(A0_draw)],lambda.hs=lambda.A0,nu.hs=nu.A0,tau.hs=tau.A0,zeta.hs=zeta.A0)
    lambda.A0 <- hs_draw$lambda
    nu.A0 <- hs_draw$nu
    tau.A0 <- hs_draw$tau
    zeta.A0 <- hs_draw$zeta
    prior.cov <- hs_draw$psi
    theta_A0[lower.tri(theta_A0)] <- prior.cov
    theta_A0[theta_A0>10] <- 10
    theta_A0[theta_A0<1e-8] <- 1e-8

    #now save training and test draws if past burn-in


    if(iter>n.burnin){
      iter_min_burnin <- iter-n.burnin

      #save mu draws
      for(mm in 1:num_outcomes){
        mu_store[iter_min_burnin,1:num_obs,mm] <- preds.train[,mm]*Ysd[mm] + Ymu[mm]
        mutest_store[iter_min_burnin,1:num_test_obs,mm] <- preds.test[,mm]*Ysd[mm] + Ymu[mm]
      }
      #save y draws




      #save sigma matrix draws
      #if mm==1, then diagonal with sigma given by independent value from dbarts initial value? or save next draw?






      Sigma_store[,,iter_min_burnin] <- (A0_draw%*%diag(as.vector(sigma.mat)^2))%*%t(A0_draw)


      if(outcome_draws==TRUE){
        #draw each individual's vector of observations from a multivariate normal
        for(obs_ind in 1:num_obs){
          temp_sample <- mvrnorm(n = 1,
                                 mu = preds.train[obs_ind,],
                                 Sigma = Sigma_store[,,iter_min_burnin])

          #Note: element-wise vector multiplication
          Y_store[iter_min_burnin,obs_ind , 1:num_outcomes] <- temp_sample*Ysd + Ymu

        }

        for(obs_ind in 1:num_test_obs){
          temp_sample <- mvrnorm(n = 1,
                                 mu = preds.test[obs_ind,],
                                 Sigma = Sigma_store[,,iter_min_burnin])

          #Note: element-wise vector  multiplication
          Ytest_store[iter_min_burnin,obs_ind , 1:num_outcomes] <- temp_sample*Ysd + Ymu

        }
      }


    }












    if(!quiet) setTxtProgressBar(pb, iter)


#
#     if(iter %% print.opt == 0){
#       print(paste("Gibbs Iteration", iter))
#       # print(c(sigma2.alpha, sigma2.beta))
#     }


  }#end iterations of Giibs sampler


  ret_list <- list()

  ret_list$mu_draws <- mu_store
  ret_list$mutest_draws <- mutest_store

  if(outcome_draws==TRUE){
    ret_list$Ytrain_draws <- Y_store
    ret_list$Ytest_draws <- Ytest_store
  }

  ret_list$Sigma_draws <- Sigma_store

  return(ret_list)



}
