Implementation of blasso
================
SAK LEE
April 29, 2018

Required packages
-----------------

**lars** has a data set. **mvtnorm** is used for calculating the multivariate normal cdf. **magrittr** is just for readibility of code.

``` r
library(lars)
library(magrittr)
library(mvtnorm)
```

Data description.
-----------------

**diabets** consists of 3 lists: 1. regular data matrix X 2. y values 3. extend version of data matrix(intersection).

``` r
data("diabetes")
diabetes %>% dim
```

    ## [1] 442   3

``` r
# Take data matrix and y values.
X <- diabetes[,1]
y <- diabetes[,2]

dim(X); length(y)
```

    ## [1] 442  10

    ## [1] 442

``` r
colnames(X)
```

    ##  [1] "age" "sex" "bmi" "map" "tc"  "ldl" "hdl" "tch" "ltg" "glu"

``` r
n <- length(y)
```

### Normalization & set params.

Data matrix and y values are decentered and set the paramerters for sigma 1 & 2, and tau.

``` r
# standardized
X <- scale(X)
y <- scale(y)

sigma1 <- 1
sigma2 <- 0.492
tau <- 4.25
```

### Enumerate the models.

Since we have 10 indep. variables, we can actually write down the all the possible 1024 models. I stored this information in matrix m.

``` r
# expand full possible models
l <- rep(list(0:1), dim(X)[2])
m <- expand.grid(l)
m %>% dim
```

    ## [1] 1024   10

``` r
head(m)
```

    ##   Var1 Var2 Var3 Var4 Var5 Var6 Var7 Var8 Var9 Var10
    ## 1    0    0    0    0    0    0    0    0    0     0
    ## 2    1    0    0    0    0    0    0    0    0     0
    ## 3    0    1    0    0    0    0    0    0    0     0
    ## 4    1    1    0    0    0    0    0    0    0     0
    ## 5    0    0    1    0    0    0    0    0    0     0
    ## 6    1    0    1    0    0    0    0    0    0     0

The elements of matrix m indicates the inclusion of the corresponding variable in the given model. So we can calculate the each number of variables in the possible model by applying the **rowSums** function to matrix m.

``` r
k_gamma <- m %>% rowSums()
k_gamma %>% head
```

    ## [1] 0 1 1 2 1 2

Calculate the weights for marginal distribution.
------------------------------------------------

First, we need to calculate the weights vector of each models to get the marginal likelihood for a particular model(equation 6 in the paper). **w\_gamma** will store the calculated weights.

``` r
lik <- log(dnorm(y, mean = 0, sd = sigma1)) %>% sum %>% exp

w_gamma <- rep(0, 2^dim(X)[2])

ptm <- proc.time() ## start of clock and then to end
for (i in 1:length(w_gamma)){
  if(k_gamma[i] == 0){
    w_gamma[i] <- 1
  } else{
    X_gamma <- X[,colnames(X)[which(m[i,] == 1)]]
    z <- expand.grid(rep(list(0:1), k_gamma[i]))
    result <- 0
    
    for (j in 1:dim(z)[1]){
      sz <- as.numeric(z[j,])
      temp <- solve(t(X_gamma) %*% X_gamma)
      mu <- temp %*% (t(X_gamma) %*% y - tau * sigma1 * sz) %>% as.numeric()
      sig <- sigma1 * temp
      
      lower_v <- rep(-Inf, length(sz))
      lower_v[sz == 1] <- 0
      
      upper_v <- rep(0, length(sz))
      upper_v[sz == 1] <- Inf
      result <- result + pmvnorm(lower = lower_v, 
                                 upper = upper_v,
                                 mean = mu, sigma = sig) /
                         dmvnorm(rep(0, k_gamma[i]), mean = mu, sigma = sig)
    }
    
    w_gamma[i] <- result
  }
  # print(paste(round(i / length(w_gamma) * 100, 2), "%"))
}
proc.time()-ptm ## end of clock.
```

    ##    user  system elapsed 
    ##  521.17    0.68  535.33

According to the equation 6, we can get a marginal likelihood for each models in **w\_gamma** vector. Since the paper assume the bernoulli prior for the model space, the posterior model distribution can be obtained by just normalizing the **w\_gamma** vector.

``` r
marginal_gamma <- w_gamma * (tau / (2*sigma1))^k_gamma * lik
posterior_gamma <- marginal_gamma / sum(marginal_gamma)
sum(posterior_gamma)
```

    ## [1] 1

``` r
posterior_gamma %>% head()
```

    ## [1] 3.405080e-46 9.529242e-44 1.076044e-46 2.290711e-44 6.275001e-15
    ## [6] 5.175410e-15

Calculating the inclusion prob.
-------------------------------

``` r
colnX <- X %>% colnames()
inclusion_p <- rep(0, 10)
for (i in 1:10){
  inclusion_p[i] <- posterior_gamma[which(as.numeric(m[,i]) == 1)] %>% sum  
}
data.frame(variableName = colnX, prob = inclusion_p)
```

    ##    variableName      prob
    ## 1           age 0.2050831
    ## 2           sex 0.8835613
    ## 3           bmi 0.9999995
    ## 4           map 0.9885209
    ## 5            tc 0.5649279
    ## 6           ldl 0.4501651
    ## 7           hdl 0.7522674
    ## 8           tch 0.4254341
    ## 9           ltg 0.9981532
    ## 10          glu 0.2639510
