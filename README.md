
# SURBART

<!-- badges: start -->
<!-- badges: end -->

The goal of SURBART is to provide an implementation of Seemingly Unrelated Regression Bayesian Additive Regression Trees, as described by Chakraborty (2016). The implementation currently does not support variable numbers of trees.



The package contains both an implementation as described by Chakraborty (2016) and implementations based on equation-by-equation estimation as described by Huber and Rossini (2021) and Huber et al. (2020). The code is mostly based on the Mixed Frequency BAVART implementation of Huber et al (2020) available at https://github.com/mpfarrho/mf-bavart .

The function surbart_original() is an implementation based on Chakraborty (2016).

The function surbart_eqbyeq_hs() is based on the equation-by-equation implementaiton with horseshoe prior for the covariances from Huber and Rossini (2021). The implementation is an edited version of that available at https://github.com/mpfarrho/mf-bavart  , except it allows for different covariate matrices for each outcome, but not mixed frequency data.

The function surbart_eqbyeq() is based on the equation-by-equation implementaiton from Huber and Rossini (2021), with the inverse Wishart covariance matrix prior from Chakraborty (2016).

Chakraborty, S. (2016). Bayesian additive regression tree for seemingly unrelated regression with automatic tree selection. In Handbook of Statistics (Vol. 35, pp. 229-251). Elsevier.

Huber, F., & Rossini, L. (2020). Inference in Bayesian additive vector autoregressive tree models. arXiv preprint arXiv:2006.16333.

Huber, F., Koop, G., Onorante, L., Pfarrhofer, M., & Schreiner, J. (2020). Nowcasting in a pandemic using non-parametric mixed frequency VARs. Journal of Econometrics.


## Installation

You can install the development version of SURBART like so:

``` r
library(devtools)
install_github("EoghanONeill/SURBART")
```

## Example

This is a basic example:

``` r
library(SURBART)
## basic example code

library(foreign)
library(systemfit)

hsb2 <- read.dta("https://stats.idre.ucla.edu/stat/stata/notes/hsb2.dta")

train_inds <- sample(1:nrow(hsb2),size = 180, replace = FALSE)
test_inds <- (1:nrow(hsb2))[-train_inds]

hsb2train <- hsb2[train_inds,]
hsb2test <- hsb2[test_inds,]


r1 <- read~female + as.numeric(ses) + socst
r2 <- math~female + as.numeric(ses) + science

fitsur <- systemfit(list(readreg = r1, mathreg = r2), data=hsb2train)

linSURpreds <- predict(fitsur, hsb2test)

sqrt(mean((linSURpreds$readreg.pred-hsb2test$read)^2))

sqrt(mean((linSURpreds$mathreg.pred-hsb2test$math)^2))



xtrain <- list()
xtrain[[1]] <- cbind(hsb2train$female, as.numeric(hsb2train$ses), hsb2train$socst)
xtrain[[2]] <- cbind(hsb2train$female, as.numeric(hsb2train$ses), hsb2train$science)

xtest <- list()
xtest[[1]] <- cbind(hsb2test$female, as.numeric(hsb2test$ses), hsb2test$socst)
xtest[[2]] <- cbind(hsb2test$female, as.numeric(hsb2test$ses), hsb2test$science)

ylist <- list()
ylist[[1]] <- hsb2train$read
ylist[[2]] <- hsb2train$math



surbartres_eqbyeq <- surbart_eqbyeq(x.train = xtrain, #either one matrix or list
               x.test = xtest, #either one matrix or list
               y = ylist,
               2,
               nrow(hsb2train),
               nrow(hsb2test),
               n.iter = 1000,
               n.burnin = 100)


readpredsurbart_eqbyeq <- apply(surbartres_eqbyeq$mutest_draws[,,1],2,mean)

mathpredsurbart_eqbyeq <- apply(surbartres_eqbyeq$mutest_draws[,,2],2,mean)

sqrt(mean((readpredsurbart_eqbyeq  -hsb2test$read)^2))

sqrt(mean((mathpredsurbart_eqbyeq- hsb2test$math)^2))



surbartres_eqbyeq_hs <- surbart_eqbyeq(x.train = xtrain, #either one matrix or list
               x.test = xtest, #either one matrix or list
               y = ylist,
               2,
               nrow(hsb2train),
               nrow(hsb2test),
               n.iter = 1000,
               n.burnin = 100)


readpredsurbart_eqbyeq_hs <- apply(surbartres_eqbyeq_hs$mutest_draws[,,1],2,mean)

mathpredsurbart_eqbyeq_hs <- apply(surbartres_eqbyeq_hs$mutest_draws[,,2],2,mean)

sqrt(mean((readpredsurbart_eqbyeq_hs  -hsb2test$read)^2))

sqrt(mean((mathpredsurbart_eqbyeq_hs- hsb2test$math)^2))




surbartres_original <- surbart_original(x.train = xtrain, #either one matrix or list
               x.test = xtest, #either one matrix or list
               y = ylist,
               2,
               nrow(hsb2train),
               nrow(hsb2test),
               n.iter = 1000,
               n.burnin = 100)


readpredsurbart_original <- apply(surbartres_original$mutest_draws[,,1],2,mean)

mathpredsurbart_original <- apply(surbartres_original$mutest_draws[,,2],2,mean)

sqrt(mean((readpredsurbart_original - hsb2test$read)^2))

sqrt(mean((mathpredsurbart_original - hsb2test$math)^2))



```

