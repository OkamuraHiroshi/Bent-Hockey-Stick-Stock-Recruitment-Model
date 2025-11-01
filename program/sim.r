## Simulation


##

logit <- function(x) log(x/(1-x))
ilogit <- function(x) 1/(1+exp(-x))
pop_update <- function(p,x,eps,GL,g=0.01,model) {
  if (model==0) res <- ifelse(x < p[,2], p[,1]*x, p[,1]*p[,2]) 
  if (model==1) res <- p[,1]*(x+sqrt(p[,2]^2+g^2/4)-sqrt((x-p[,2])^2+g^2/4))
  if (model==2) {
    a <- p[,1]
    b <- p[,2]
    pred <- function(a,b,x) sum(GL$xr_x*GL$wt*a*b*(x/b)^(1-(x/b)^(1/GL$xp_x)))/(2*GL$xr_x)
    pred <- Vectorize(pred, c("a","b","x"))
    res <- ifelse(x < b, pred(a,b,x), a*b)
  }
  if (model==3) res <- p[,1]*x/(1+x/p[,2])
  if (model==4) res <- p[,1]*x*exp(-x/p[,2])

  log(res)+eps
}
Pop_simulator <- function(parms=NULL, res, model=1, Sim=1000, seed0=1){
  set.seed(seed0)
  require(statmod)
  dat <- res$dat
  obj <- res$obj
  sdrep <- res$sdrep
  model <- res$model
  g <- res$g
  x_int <- res$x_int
  n_g <- res$n_g
  GL <- res$GL
  M <- res$M

  if (is.null(parms)){
    mu_p <- sdrep$par.fixed
    Sigma_p <- sdrep$cov.fixed
    parms <- rmvnorm(Sim, mu_p, Sigma_p)
  }
  pars <- as.data.frame(parms)
  Sim <- nrow(pars)
  nY <- dat$Y
  nT <- length(dat$CPUE)
  n0 <- n_last <- matrix(NA, nrow=nY, ncol=Sim)
  n <- u <- matrix(NA, nrow=nT, ncol=Sim)
  log_a <- pars$log_a
  log_b <- pars$log_b
  a <- exp(log_a)
  b <- exp(log_b)
  base_parms <- cbind(a,b)
  sigma <- exp(pars$log_sigma)
  eps <- matrix(rnorm(nY*Sim,0,sigma),nrow=nY,ncol=Sim,byrow=TRUE)
  YEAR <- dat$YEAR+1
  START <- dat$START
  END <- dat$END
  WEEK <- dat$WEEK
  U <- summary(sdrep)[rownames(summary(sdrep))=="U",]
  logit_U <- cbind(logit(U[,1]), (1/(U[,1]*(1-U[,1])))*U[,2])
  set.seed(seed0)
  
  for (i in 1:nT){
      u[i,] <- ilogit(rnorm(Sim, logit_U[i,1], logit_U[i,2]))
      if (START[i]==1) {
        if (i==1) n0[YEAR[i],] <- pars$tilde_n0 else n0[YEAR[i],] <- pop_update(base_parms,exp(n_last[YEAR[i]-1,]),eps[YEAR[i]-1,],GL,g=g,model=model)
        n[i,] <- n0[YEAR[i],]-M*7*WEEK[i]
      }
      if (START[i]==0) n[i,] <- n[i-1,]-M*7*(WEEK[i]-WEEK[i-1])+log(1-u[i-1,])
      if (END[i]==1) n_last[YEAR[i],] <- n[i,]-M*7*(26-WEEK[i])+log(1-u[i,])
  }
  n0_new <- pop_update(base_parms,exp(n_last[nY,]),eps[nY,],GL,g=g,model=model)
  list(parms=pars, base_parms=base_parms, M=M, WEEK=WEEK, n0=n0, n_last=n_last, n0_new=n0_new, n=n, u=u)
}

#

Simulated_data <- function(res, Sim=1000, p=0.1, seed0=1){
  dat <- res$dat

  YEAR <- dat$YEAR+1
  START <- dat$START
  END <- dat$END
  WEEK <- dat$WEEK
  M <- res$M
  model <- res$model
  sdrep <- res$sdrep
  
  set.seed(seed0)
  
  pops <- Pop_simulator(NULL, res, model, Sim, seed0)
  gen_par <- pops$parms
  nY <- dat$Y
  est_n0 <- summary(sdrep)[rownames(summary(sdrep))=="n0",]
  
  W <- sapply(1:Sim, function(i) 1/nY*sum(dnorm(pops$n0[,i],est_n0[,1],est_n0[,2],log=TRUE)))
  W[is.na(W)] <- min(W, na.rm=TRUE)-abs(min(W, na.rm=TRUE))
  W <- W-max(W)
  P <- exp(W)/sum(exp(W))
  m <- round(p*Sim)
  id <- sample(Sim,m,prob=P)
  nt <- length(dat$CPUE)
  
  parms <- gen_par[id,]
  base_parms <- pops$base_parms[id,]
  n <- pops$n[,id]
  u <- pops$u[,id]
  n0 <- pops$n0[,id]
  n_last <- pops$n_last[,id]
  n0_new <- pops$n0_new[id]
  
  catch <- exp(n)*u
    
  log_tilde_q <- parms[,"log_tilde_q"]
  tau <- exp(parms[,"log_tau"])
  eta <- exp(parms[,"log_eta"])
  z <- matrix(rnorm(nt*m), nrow=nt, ncol=m)
  z2 <- matrix(rnorm(nY*m), nrow=nY, ncol=m)
  log_cpue <- matrix(NA, nrow=nt, ncol=m)
  log_q <- matrix(NA, nrow=nY, ncol=m)
  log_q[1,] <- rnorm(m, log_tilde_q, eta)
  
  for (i in 1:nt){
    if (i>1 & START[i]==1) log_q[YEAR[i],] <- log_q[YEAR[i]-1,]+eta*z2[YEAR[i]-1,]
    log_cpue[i,] <- log_q[YEAR[i],]+n[i,]+tau*z[i,]
  }
  
  tcat <- sapply(1:m, function(i) tapply(catch[,i], YEAR, sum))
  
  list(m=m, parms=parms, base_parms=base_parms, model=model, M=M, YEAR=YEAR-1, WEEK=WEEK, START=START, END=END, n0=n0, n_last=n_last, n0_new=n0_new, n=n, u=u, cpue=exp(log_cpue), catch=catch, tcat=tcat, q=exp(log_q))
}

##

make_dat <- function(x,i) data.frame(Year=x$YEAR, Week=x$WEEK,Cat=x$catch[,i], CR=x$cpue[,i])
make_par <- function(x,i){
  list(log_a=x$parms[i,"log_a"],
       log_b=x$parms[i,"log_b"],
       tilde_n0=x$parms[i,"tilde_n0"],
       log_tilde_q=x$parms[i,"log_tilde_q"],
       log_sigma=x$parms[i,"log_sigma"],
       log_tau=x$parms[i,"log_tau"],
       log_eta=x$parms[i,"log_eta"],
       log_q=log(x$q[,i]),
       n0=x$n0[,i]
      )
}
sim2est <- function(Sim_dat,model=1,a_init=NULL,b_init=NULL,a_range=NULL,b_range=NULL,mod="BHS",start=1){
  Res_sim <- list()
  m <- Sim_dat$m
  last_n <- last(Sim_dat$n)
  last_n <- ifelse(is.na(last_n),-Inf,last_n)
  a_init0 <- a_init
  
  if (model==2) X = rbind(c(0,5),c(0,6),c(0,7),c(0,8),c(0,9),c(0,10),c(0,20),c(0,4),c(0,3),c(0,2),c(0,1),c(0,0.5),c(0,0.1))
  if (model==1) X = rbind(c(0,0.01),c(0,0.1),c(0,0.5))
  if (model==0) X = rbind(c(0,0.01))
  if (model>=3) X = rbind(c(0,1))
  
  nr <- nrow(X)
  
  temp_res <- list(conv_diag=rep(TRUE,5))
  
  for (i in start:m){
    print(i)
    for (kk in 1:nr){
    if (mod=="HS"){
      if (model==0) if (i==7 | i==21 | i==23 | i==31 | i==67 | i==68) a_init=4 else {if (i==15 | i==48 | i==85 | i==95) a_init=2 else {if (i==13 | i==29 | i==57 | i==93 | i==97 | i==98) a_init=6 else {if (i==59) a_init=1 else a_init=a_init0}}}
      if (model==1) if (i==89) a_init=4 else a_init=a_init0
      if (model==2) {
        if (i==32 | i==45 | i==52 | i==64 | i==88 | i==89) a_init=2 else {if (i==65 | i==86) a_init=4 else a_init=a_init0}
      }
      if (model==3) if (i==6 | i==7 | i==23 | i==30 | i==48 | i==62 | i==77) a_init=2 else {if (i==32 | i==76) a_init=4 else {if (i==45 | i==53 | i==54) a_init=6 else {if (i==51) a_init=8 else a_init=a_init0}}}
      if (model==4) if (i==1 | i==5 | i==6 | i==10 | i==11 | i==14 | i==15 | i==17 | i==30 | i==33 | i==42 | i==77 | i==80 | i==89) a_init=4 else {if (i==25 | i==63 | i==95) a_init=6 else {if (i==38 | i==55 | i==75 | i==92) a_init=8 else a_init=a_init0}}
    }
    if (mod=="MR"){ }
    if (mod=="BHS"){
      if (model==0) if (i %in% c(3,10,13,26,42,53,57,62,77,78,90,92,97,100)) a_init=2 else {if (i %in% c(11,72)) a_init=4 else {if (i %in% c(32,35,51,58,70,76,79,94,95)) a_init=6 else {if (i %in% c(6,16,19,21,39)) a_init=8 else {if (i==15) a_init=1 else a_init=a_init0}}}}
      if (model==1) if (i==41) a_init=1 else {if (i==51 | i==68 | i==84 | i==90) a_init=4 else a_init=a_init0}
      if (model==2) if (i==16 | i==45 | i==72 | i==76 | i==94) a_init=2 else {if (i==81 | i==92) a_init=4 else a_init=a_init0}
      if (model==3) if (i %in% c(1,52,53,71,72,84,92,98)) a_init=2 else {if (i==50) a_init=4 else {if (i==36) a_init=6 else a_init=a_init0}}
    #  if (model==3) if (i %in% c(12,34,61,89)) a_init=2 else a_init=1.2
      if (model==4) if (i %in% c(7,31,33,43,44,77)) a_init=4 else {if (i %in% c(11,18,22)) a_init=6 else {if (i %in% c(39,41,58)) a_init=8 else {if (i==1) a_init=10 else a_init=a_init0}}}
    }
    if (mod=="BH"){
      if (model==0) if (i %in% c(4,23,34,41,47,51,67,68,69,79,86)) a_init=2 else {if (i %in% c(43,72)) a_init=4 else {if (i %in% c(15,32)) a_init=6 else {if (i %in% c(28,44,63,76,81,90)) a_init=8 else {if (i==100) a_init=1 else a_init=a_init0}}}}
      if (model==1) if (i %in% c(2,36,76)) a_init=4 else {if (i==90) a_init=2 else a_init=a_init0}
      if (model==2) if (i %in% c(14,33,40,43,52,54,58)) a_init=2 else {if (i %in% c(27)) a_init=1 else {if (i %in% c(29,98,99)) a_init=4 else {if (i %in% c(44)) a_init=6 else a_init=a_init0}}}
      if (model==4) if (i %in% c(3,4,7,9,22,43,57,65,75,86)) a_init=4 else {if (i==60) a_init=6 else a_init=a_init0}
    }
    if (mod=="RI"){
      if (model==0) if (i %in% c(22,55,83,84)) a_init=4 else {if (i %in% c(17,25,41,48)) a_init=6 else {if (i %in% c(33,66,80,85)) a_init=8 else {if (i %in% c(75)) a_init=10 else {if (i %in% c(52,76,77,82)) a_init=12 else {if (i %in% c(26)) a_init=14 else a_init=a_init0}}}}}
      if (model==1) if (i %in% c(41,47,73,76,85)) a_init=8 else {a_init=a_init0}
      if (model==2) if (i==25) a_init=4 else {if (i==77) a_init=2 else a_init=a_init0}
      if (model==3) if (i %in% c(9,55,75)) a_init=2 else a_init=a_init0
      if (model==4) if (i %in% c(100)) a_init=8 else {a_init=a_init0}
    }
    if (exp(last_n[i]) > 0){
      dat_sim <- make_dat(Sim_dat,i)
      par_sim <- make_par(Sim_dat,i)
      if (kk==1 | (kk > 1 & temp_res$conv_diag[5]!=TRUE)) temp_res <- leslie(dat_sim, par_sim, model=model, do_compile=FALSE, pre_process=FALSE, x_int=X[kk,], g0=X[kk,2],a_range=a_range,b_range=b_range,a_init=a_init,b_init=b_init)
    } else {
      temp1 <- matrix(0, nrow=4, ncol=2)
      rownames(temp1) <- c("b","n0_new","log_r","q")
      temp1[4,1] <- 0.000001
      colnames(temp1) <- c("Estimate","Std. Error")
      conv_diag <- c("mod_conv"=NA, "pdHess"=NA, "vcov"=NA, "convergence"=NA)
      temp_res <- list(sdrep=temp1, conv_diag=conv_diag)
    }
    }
    Res_sim[[i]] <- temp_res
  }
  return(Res_sim)
}
