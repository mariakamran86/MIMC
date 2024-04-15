#                                    ********************************************************
#                                        Multiple Indicator Multiple Cause (MIMIC) model 
#                                                  with Discrete Indicators
#
#                                    *********************************************************

# Preliminary----
options(java.parameters="-Xmx5000m")
set.seed(300)
library(maxLik)
library(dplyr)
library(haven)
library(dummies)
library(xtable)


# I.   Data----

WS <- read_dta("C:/Users/path")#Loading data
WS <- WS[complete.cases(WS),]# to remove observations with na



# setting the data                                        

y1   <- as.numeric(WS$y1)
y2 <- as.numeric(WS$y2)
x1 <- as.numeric(WS$x1)
x2 <- as.numeric(WS$x2)
x3 <- as.numeric(WS$x3)
z1 <- as.numeric(WS$z1)
z2 <- as.numeric(WS$z2)
d1 <- as.numeric(WS$d1)
d2 <- as.numeric(WS$d2)
d3 <- as.numeric(WS$d3)


# II.  Initial Values for Maximum likelihood----

# In this step, initial values for the maximum likelihood estimation using a ordinary
# stage least square method is estimated. In this case since it is multi-choice-multi-indicators
# 3SLS is used. Other options include 2sls and GMM.

library(systemfit)


eq1 <- y1 ~ x1 + x2 + x3 + d1 + d2 + d3
eq2 <- y2 ~ x1 + x2 + x3 + d1 + d2 + d3
eq3 <-  x2 ~ x1  + x3 + d1 + d2 + d3 + z1 + z2


system <- list(eq1=eq1,
               eq2=eq2,
               eq3=eq3)


aa <-
  systemfit(
    system,
    "3SLS",
    data = WS,
    inst =  ~ x1 + x3 + d1 + d2 + d3 + z1 + z2,
    pooled = TRUE,
    methodResidCov = "noDfCor",
    residCovWeighted = TRUE,
    method3sls = "GMM"
  )

summary(aa)
logLik(aa)



names(aa$coefficients)[names(aa$coefficients) == "eq1_(Intercept)"] <- "b0"
names(aa$coefficients)[names(aa$coefficients) == "eq1_x1"]          <- "b1"
names(aa$coefficients)[names(aa$coefficients) == "eq1_x2" ]         <- "b2"
names(aa$coefficients)[names(aa$coefficients) == "eq1_x3"]          <- "b3"
names(aa$coefficients)[names(aa$coefficients) == "eq1_d1"]          <- "bl1"
names(aa$coefficients)[names(aa$coefficients) == "eq1_d2"  ]        <- "bl2"
names(aa$coefficients)[names(aa$coefficients) == "eq3_d3"]          <- "bl3"


names(aa$coefficients)[names(aa$coefficients) == "eq2_(Intercept)"] <- "c0"
names(aa$coefficients)[names(aa$coefficients) == "eq2_x1"]          <- "c1"
names(aa$coefficients)[names(aa$coefficients) == "eq2_x2" ]         <- "c2"
names(aa$coefficients)[names(aa$coefficients) == "eq2_x3"]          <- "c3"
names(aa$coefficients)[names(aa$coefficients) == "eq2_d1"]          <- "cl1"
names(aa$coefficients)[names(aa$coefficients) == "eq2_d2"  ]        <- "cl2"
names(aa$coefficients)[names(aa$coefficients) == "eq2_d3"]          <- "cl3"


names(aa$coefficients)[names(aa$coefficients) == "eq3_(Intercept)"] <- "g0"
names(aa$coefficients)[names(aa$coefficients) == "eq3_x1"]          <- "g1"
names(aa$coefficients)[names(aa$coefficients) == "eq3_x3"]          <- "g3"
names(aa$coefficients)[names(aa$coefficients) == "eq3_z1"]          <- " gz1"
names(aa$coefficients)[names(aa$coefficients) == "eq3_z2" ]         <- "gz2"
names(aa$coefficients)[names(aa$coefficients) == "eq3_d1"]          <- "gl1"
names(aa$coefficients)[names(aa$coefficients) == "eq3_d2"  ]        <- "gl2"
names(aa$coefficients)[names(aa$coefficients) == "eq3_d3"]          <- "gl3"

# III. Maximum Likelihood ----

maxLik_run <- function(r, w, q, a, k, s, f) {
  tryCatch({
    
    llr <- function(param) {
      
      #Param----
      b0  <-  param[1]
      b1  <-  param[2]
      b2  <-  param[3]
      b3  <-  param[4]
      bl1  <-  param[5]
      bl2  <-  param[6]
      bl3  <-  param[7]
      
      c0  <-  param[8]
      c1  <-  param[9]
      c2  <-  param[10]
      c3  <-  param[11]
      cl1 <-  param[12]
      cl2 <-  param[13]
      cl3 <-  param[14]
      
      g0 <- param[15]
      g1 <- param[16]
      g3 <- param[17]
      gz1 <- param[18]
      gz2 <- param[19]
      gl1  <-  param[20]
      gl2  <-  param[21]
      gl3  <-  param[22]
      
      q1 <- param[23]
      q2 <- param[24]
      rho1 <- param[25]
      rho2 <- param[26]
      rho3 <- param[27]
      gamma <- param[28]
  
  #ML Function----
  sum(
    log((exp(q1) / (exp(q1) + exp(q2) + 1))
        *( ((pnorm(b0+b1*x1+ b2*x2+ b3*x3+ bl1*d1 + bl2*d2+ bl3*d3))^y1)
           * ((pnorm(c0+c1*x1+ c2*x2+ c3*x3+ cl1*d1 + cl2*d2+ cl3*d3))^y2)
           *((pnorm(g0+g1*x1+ g3*x3+gl1*d1+ gl2*d2 + gl3*d3 +gz1*z1+gz2*z2))^x2)
           * ((1-pnorm(b0+b1*x1+ b2*x2+ b3*x3+ bl1*d1 + bl2*d2+ bl3*d3))^(1-y1))
           *((1-pnorm(c0+c1*x1+ c2*x2+ c3*x3+ cl1*d1 + cl2*d2+ cl3*d3))^(1-y2))
           * ((1-pnorm(g0+g1*x1+ g3*x3+gl1*d1+ gl2*d2 + gl3*d3 +gz1*z1+gz2*z2))^(1-x2)))
        +   (exp(q2)/(exp(q1)+exp(q2)+1)) 
        *(
          ((pnorm(b0+b1*x1+ b2*x2+ b3*x3+ bl1*d1 + bl2*d2+ bl3*d3+rho1))^y1)
          * ((pnorm(c0+c1*x1+ c2*x2+ c3*x3+ cl1*d1 + cl2*d2+ cl3*d3+rho2))^y2)
          *((pnorm(g0+g1*x1+ g3*x3+gl1*d1+ gl2*d2 + gl3*d3 +gz1*z1+gz2*z2+rho3))^x2)
          * ((1-pnorm(b0+b1*x1+ b2*x2+ b3*x3+ bl1*d1 + bl2*d2+ bl3*d3+rho1))^(1-y1))
          *((1-pnorm(c0+c1*x1+ c2*x2+ c3*x3+ cl1*d1 + cl2*d2+ cl3*d3+rho2))^(1-y2))
          * ((1-pnorm(g0+g1*x1+ g3*x3+gl1*d1+ gl2*d2 + gl3*d3 +gz1*z1+gz2*z2+rho3))^(1-x2)) )
        + (1-(exp(q1)/(exp(q1)+exp(q2)+1))-(exp(q2)/(exp(q1)+exp(q2)+1)))
        *
    (
      ((pnorm(b0+b1*x1+ b2*x2+ b3*x3+ bl1*d1 + bl2*d2+ bl3*d3+rho1*(exp(gamma)/(exp(gamma)+1))))^y1)
      * ((pnorm(c0+c1*x1+ c2*x2+ c3*x3+ cl1*d1 + cl2*d2+ cl3*d3+rho2*(exp(gamma)/(exp(gamma)+1))))^y2)
      *((pnorm(g0+g1*x1+ g3*x3+gl1*d1+ gl2*d2 + gl3*d3 +gz1*z1+gz2*z2*rh3**(exp(gamma)/(exp(gamma)+1))))^x2)
      * ((1-pnorm(b0+b1*x1+ b2*x2+ b3*x3+ bl1*d1 + bl2*d2+ bl3*d3+rho1*(exp(gamma)/(exp(gamma)+1))))^(1-y1))
      *((1-pnorm(c0+c1*x1+ c2*x2+ c3*x3+ cl1*d1 + cl2*d2+ cl3*d3+rho2*(exp(gamma)/(exp(gamma)+1))))^(1-y2))
      * ((1-pnorm(g0+g1*x1+ g3*x3+gl1*d1+ gl2*d2 + gl3*d3 +gz1*z1+gz2*z2+rho3*(exp(gamma)/(exp(gamma)+1))))^(1-x2)) )))
}
  
    return(maxLik(
      llr,
      start = c(
        aa$coefficients,
        rho1 = r,
        rho2 = w,
        rho3 = q,
        rho4 = a,
        q1 = k,
        q2 = s,
        gamma = f
      ),
      method = "bfgs"
      ,
      control = list(
        tol = -1,
        reltol = 1e-12,
        gradtol = 1e-12
      )
    ))
  },
error = function(e)
  return(NA))  
  }

# Assuming that residual of all three equations are mutually uncorrelated, 
# rhoi accounts for the compound coefficient for error variance of i equation.
# Here qi is the associated probability of each specification. 
# gamma here is a trend of the error variance
# choosing these values is an iterative process done with the following function.


st <-
  expand.grid(
    rho1 = c(0.09),
    rho2 = c(1.251976),
    rho3 = c(0.5),
    q1 = c(0.1),
    q2 = c(-3.42),
    gamma = c(-3.26)
  )

maxLik_list <-
  Map(maxLik_run,
      st$rho1,
      st$rho2,
      st$rho3,
      st$q1,
      st$q2,
      st$gamma) 

# Exporting results to latex

latex <- function(dta) {
  hmn <- as.data.frame(summary(dta)[6])
  hmn <- hmn[!grepl("*l", row.names(hmn)),-c(3 , 4)]
  dd  <- xtable(hmn[1:2], digits = 3)
  return(dd)
}


#III.A Residual Correlation----

rcor<-function(dta){
  rho1 <- dta$estimate["rho1"]
  rho2 <- dta$estimate["rho2"]
  rho3 <- dta$estimate["rho3"]
  q1   <- dta$estimate["q1"]
  q2   <- dta$estimate["q2"]
  gamma<- dta$estimate["gamma"]
  
  c_gamma<-exp(gamma)/(exp(gamma)+1)
  
  c_q1<-exp(q1)/(1+exp(q1)+exp(q2))
  c_q2<-exp(q2)/(1+exp(q1)+exp(q2))
  c_q3<-1-c_q1-c_q2
  
  S2<-c_q1+c_q3*(c_gamma^2)-(c_q1+c_q3*c_gamma)^2
  
  rho12<-rho1*rho2*S2/sqrt((rho1^2*S2+1)*(rho2^2*S2+1))
  rho13<-rho1*rho3*S2/sqrt((rho1^2*S2+1)*(rho3^2*S2+1))

  rho23<-rho2*rho3*S2/sqrt((rho2^2*S2+1)*(rho3^2*S2+1))

  
  result<-rbind(rho12,rho13,rho23,c_q1,c_q2,c_q3,c_gamma)
  
  return(result)
  
}

b<-as.matrix(dta$estimate)
V<-vcov(dta)


pload_3 <- function(b, V,i,j) {
  r1 <- b[paste0("rho", i),1]
  r2 <- b[paste0("rho", j),1]
  c1 <- b["gamma",1]
  p1 <- b["q1",1]
  p2 <- b["q2",1]
  
  b0 <- rbind(1, r1, r2, p1, c1, p2)
  
  q <- V[c("p1","p2"),c("p1","p2")]
  w <- rbind(V[paste0("rho", i), c("p1","p2")], V[paste0("rho", j), c("p1","p2")])
  h <- rbind(cbind(V[paste0("rho", i),paste0("rho", i)],V[paste0("rho", i),paste0("rho", j)]),
             cbind(V[paste0("rho", j),paste0("rho", i)],V[paste0("rho", j),paste0("rho", j)]))
  V0 <- rbind(cbind(h, w), cbind(w, q))
  
  return(V0)
}


pload_g<-function(i,j){
  r1<-dta$gradient[paste0("rho",i)]
  r2<-dta$gradient[paste0("rho",j)]
  q1<-dta$gradient["q1"]
  q2<-dta$gradient["q2"]
  gamma<-dta$gradient["gamma"]
  
  gamma1<-exp(gamma)/(exp(gamma)+1)
  p1<-exp(q1)/(1+exp(q1)+exp(q2))
  p2<-exp(q2)/(1+exp(q1)+exp(q2))
  p3<-1-p1-p2
  S2<-p1+p3*gamma1^2-(p1+p3)^2
  grho<-r1*r2*S2/sqrt((r1^2*S2+1)*(r2^2*S2+1))
  G<-rbind(grho,grho,p1,p2)
  return(G)
  
}


#Standard Error for residual correlation rho_12
V0<-pload_3(b,V,1,2)
G<-pload_g(1,2)

rcor(dta)["rho12",]/(t(G)%*%V0%*%G)

#Standard Error for residual correlation rho_13
V0<-pload_3(b,V,1,3)
G<-pload_g(1,3)

rcor(dta)["rho13",]/(t(G)%*%V0%*%G)


#Standard Error for residual correlation rho_23
V0<-pload_3(b,V,2,3)
G<-pload_g(2,3)

rcor(dta)["rho23",]/(t(G)%*%V0%*%G)

#IV Bootstraping the estimates----

library(boot)
library(parallel)

gmm_beta <- function(formula, data, indices) {
  d <- data[indices,]
  fit <- maxLik_run(rho1=r,rho2=w,rho3=q,rho4=a,q1=k,q2=s,
                    gamma=f)
  #return(coef(fit))
  return(logLik(fit))
}

gmm_beta_boot <-
  boot(
    data = WS,
    statistic = gmm_beta,
    R = 1342,
    formula = list(r1, r2, ...),
    parallel = "snow"
  )

boot1000 <-
  as.data.frame(cbind(gmm_beta_boot$t0, confint(gmm_beta_boot, type = "bca")))

bootlik <- as.data.frame(gmm_beta_boot$t)

#gmm_beta_jack<-jack.after.boot(gmm_beta_boot, useJ = FALSE, stinf = FALSE,index = 7)