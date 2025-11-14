### State-Space Leslie Removal Method with the BHS SR function

## Population dynamics model

leslie <- function(
  dat,
  parms=NULL,
  model=2,
  do_compile=TRUE,
  tmb_file="leslie_r3",
  M=0.0133,     # Roa-Ureta & Arkhipkin (2007)
  g0=0.01,
  x_int=c(0,5),
  n_g=10,
  a_range=c(1,10),
  b_range=c(1,10),
  a_init=NULL,
  b_init=NULL,
  max_P=0.5,
  maxit=1000,
  maxeval=1000,
  silent=TRUE,
  pre_process=TRUE,
  phase_a=FALSE,
  phase_b=FALSE,
  no_est=FALSE,
  make_plot=TRUE,
  optimizer="nlminb",
  method="BFGS"
){
  require(TMB)
  require(tidyverse)
  require(statmod)

  # Gauss-Legendre integration
  
  g_le <- gauss.quad(n_g,kind="legendre")
  
  nodes <- g_le$nodes
  wt <- g_le$weights
  
  #

  if (pre_process){
    dat <- dat %>% mutate(Cat=Catch/Mass, CR=CPUE/Mass)
    TC <- dat %>% group_by(Year) %>% summarize(TC=sum(Cat))
    Years <- unique(dat$Year)
    leslie_res <- lapply(Years, function(i) lm(CR~cumsum(Cat),data=subset(dat, Year==i)))
    lm_res <- sapply(1:length(leslie_res), function(i) c(-leslie_res[[i]]$coef[2], -leslie_res[[i]]$coef[1]/leslie_res[[i]]$coef[2]))
    lm_res <- rbind(lm_res, TC$TC)
    lm_res <- as.data.frame(t(lm_res))
    names(lm_res) <- c("q","N0","Catch")
    rownames(lm_res) <- Years
  }
  
  if (do_compile){
    compile(paste0(tmb_file,".cpp"))
    dyn.load(dynlib(tmb_file))
  }
  
  model_name <- c("HS","MR","BHS","BH","RI")
  
  dat1 <- dat %>% mutate(Y=Year-min(Year)) %>% group_by(Year) %>% mutate(W=Week-min(Week), START=(W==0), END=(W==max(W)))
  dat1c <- dat %>% group_by(Year) %>% summarize(CC=sum(Cat))
  nY <- max(dat1$Y)+1
  dat_bhs <- list(CPUE=dat1$CR, CAT=dat1$Cat, WEEK=dat1$Week, YEAR=dat1$Y, START=as.numeric(dat1$START), END=as.numeric(dat1$END), Y=nY, M=M, x_lo=x_int[1], x_up=x_int[2], nodes=nodes, wt=wt, g=g0, model=model)
  if (is.null(parms)){
    par_bhs <- list(
      log_a=2,
      log_b=2,
      tilde_n0=max(mean(log(lm_res$N0)), max(log(dat1c$CC*1.5))),
      log_tilde_q=mean(log(lm_res$q)),
      log_sigma=log(0.3),
      log_tau=log(0.3),
      log_eta=log(0.3),
      log_q=log(lm_res$q),
      n0=pmax(log(lm_res$N0), log(dat1c$CC*1.5))
    )
  } else par_bhs <- parms
  
  if (!is.null(a_init)) par_bhs$log_a <- a_init
  if (!is.null(b_init)) par_bhs$log_b <- b_init

  if (par_bhs$log_b < min(par_bhs$n0)) par_bhs$log_b <- 1.1*min(par_bhs$n0)
  if (par_bhs$log_b > max_P*max(par_bhs$n0)) par_bhs$log_b <- max_P*(max(par_bhs$n0)+0.1*par_bhs$log_b)
 
  lo <- c(a_range[1],b_range[1],rep(-20,5))
  up <- c(a_range[2],b_range[2],rep(20,5))
  
  maps <- list()
  if (phase_a) maps$log_a <- factor(NA)
  if (phase_b) maps$log_b <- factor(NA)
  obj <- MakeADFun(dat_bhs, par_bhs, map=maps, random=c("log_q","n0"), DLL=tmb_file, silent=silent)
  kk <- 1
  dd <- c(seq(0.01,0.5,len=5),seq(-0.5,-0.01,len=5))
  while(is.na(as.numeric(obj$fn())) & kk <= 10){
    par_bhs$log_b <- par_bhs$log_b+dd[kk]
    obj <- MakeADFun(dat_bhs, par_bhs, map=maps, random=c("log_q","n0"), DLL=tmb_file, silent=silent)
    kk <- kk + 1
  }
  if (no_est) mod_bhs <- list(par=obj$par, convergence=1) else {
    if (optimizer=="optim") {
        if (method=="L-BFGS-B") mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit), lower=lo, upper=up) else mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit))
    }
    if (optimizer=="nlminb") mod_bhs <- nlminb(obj$par, obj$fn, obj$gr, lower=lo, upper=up, control=list(iter.max=maxit, eval.max=maxeval))
    if (phase_a){
      par_bhs <- obj$env$parList(mod_bhs$par)
      maps$log_a <- NULL
      obj <- MakeADFun(dat_bhs, par_bhs, map=maps, random=c("log_q","n0"), DLL=tmb_file, silent=silent) 
      if (optimizer=="optim") {
        if (method=="L-BFGS-B") mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit), lower=lo, upper=up) else mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit))
      }
      if (optimizer=="nlminb") mod_bhs <- nlminb(obj$par, obj$fn, obj$gr, lower=lo, upper=up, control=list(iter.max=maxit, eval.max=maxeval))
    }
    if (phase_b){
      par_bhs <- obj$env$parList(mod_bhs$par)
      maps$log_b <- NULL
      obj <- MakeADFun(dat_bhs, par_bhs, map=maps, random=c("log_q","n0"), DLL=tmb_file, silent=silent) 
      if (optimizer=="optim") {
        if (method=="L-BFGS-B") mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit), lower=lo, upper=up) else mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit))
      }
      if (optimizer=="nlminb") mod_bhs <- nlminb(obj$par, obj$fn, obj$gr, lower=lo, upper=up, control=list(iter.max=maxit, eval.max=maxeval))
    }
  }
  
  sdrep <- sdreport(obj)
  
  conv_diag <- c("mod_conv" = (mod_bhs$convergence==0), "grad" = (max(abs(sdrep$gradient.fixed)) < 10^(-2)), "pdHess" = (sdrep$pdHess==TRUE), "vcov" = (all(is.na(sdrep$cov.fixed))==FALSE))
  conv_diag <- c(conv_diag, "convergence"=all(conv_diag))
  
  if (rev(mod_bhs$par)[1] < -10){
    par_bhs$log_eta <- -10
    maps$log_eta <- factor(NA)
    lo <- lo[1:(length(lo)-1)]
    up <- up[1:(length(up)-1)]    
    obj <- MakeADFun(dat_bhs, par_bhs, map=maps, random=c("log_q","n0"), DLL=tmb_file, silent=silent)
    if (optimizer=="optim") {
      if (method=="L-BFGS-B") mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit), lower=lo, upper=up) else mod_bhs <- optim(obj$par, obj$fn, obj$gr, method=method, control=list(maxit=maxit))
    }
    if (optimizer=="nlminb") mod_bhs <- nlminb(obj$par, obj$fn, obj$gr, lower=lo, upper=up, control=list(iter.max=maxit, eval.max=maxeval))
    
    sdrep <- sdreport(obj)
    conv_diag <- c("mod_conv" = (mod_bhs$convergence==0), "grad" = (max(abs(sdrep$gradient.fixed)) < 10^(-2)), "pdHess" = (sdrep$pdHess==TRUE), "vcov" = (all(is.na(sdrep$cov.fixed))==FALSE))
    conv_diag <- c(conv_diag, "convergence"=all(conv_diag))
  }

  pars <- obj$env$parList(mod_bhs$par)
  pars$a <- obj$report()$a
  pars$b <- obj$report()$b
  
  xm_x <- 0.5*(x_int[2]+x_int[1])
  xr_x <- 0.5*(x_int[2]-x_int[1])
  dx_gi_x <- xr_x*nodes
  xp_x <- xm_x+dx_gi_x
  
  GL <- list(xr_x=xr_x, xp_x=xp_x, wt=wt)

  N_E <- summary(sdrep)[rownames(summary(sdrep))=="N_E",1]
  N_S <- summary(sdrep)[rownames(summary(sdrep))=="N_S",1]
  
  ts_out <- data.frame(S=N_E, R=N_S)
  
  if (make_plot){
    X <- seq(0,round(max(N_E)*1.1))
    pred_out <- data.frame(ES = X, ER = pred_R(X, pars, GL, g0, model=model))

    p1 <- ggplot(ts_out, aes(x=S,y=R))+geom_point()+labs(x="S", y="R")+theme_bw()+geom_line(data=pred_out, aes(x=ES,y=ER),color="blue")
  } else p1 <- NULL
  
  list(model=model, model_name=model_name[model+1], M=M, x_int=x_int, n_g=n_g, GL=GL, no_est=no_est, g=g0, phase_b=phase_b, conv_diag=conv_diag, p=p1, obj=obj, mod=mod_bhs, pars=pars, dat=dat_bhs, sdrep=sdrep, ts_out=ts_out)
}

## Prediction

pred_R <- function(S,p,GL,g=0.01,model=0) {
   a <- p$a
   b <- p$b
   
   if (model==0) res <- ifelse(S < b, a*S, a*b) 
   if (model==1) res <- a*(S+sqrt(b^2+g^2/4)-sqrt((S-b)^2+g^2/4))
   if (model==2){
     pred <- function(S) sum(GL$xr_x*GL$wt*a*b*(S/b)^(1-(S/b)^(1/GL$xp_x)))/(2*GL$xr_x)
     pred <- Vectorize(pred, "S")
     res <- ifelse(S < b, pred(S), a*b)
   }
   if (model==3) res <- a*S/(1+S/b)
   if (model==4) res <- a*S*exp(-S/b)
   
   return(res)
}

SR_plot <- function(res, msy_res=NULL){
 ts_out <- res$ts_out
 model <- res$model
 pars <- res$par
 GL <- res$GL
 g <- res$g
 
 X <- seq(0,round(max(ts_out$S)*1.1))
 pred_out <- data.frame(ES = X, ER = pred_R(X, pars, GL, g=g, model=model))

 p1 <- ggplot(ts_out, aes(x=S,y=R))+geom_point()+labs(x="S", y="R")+theme_bw()+geom_line(data=pred_out, aes(x=ES,y=ER),color="blue")
 
 p1
}

## MSY estimation

n_est <- function(n, p, M, GL, g, U=0.5, n_g=50, stochastic=TRUE, model=2){
   require(statmod)
   g_he <- gauss.quad(n_g,kind="hermite")

   if (stochastic) sigma <- exp(p$log_sigma) else sigma <- rep(0, n_g)
   
   S <- function(z) exp(n-M*7*26+log(1-U)+z)
   integrand <- function(z) pred_R(S(z), p, GL, g, model=model)
   n_pred <- sum(g_he$weights*integrand(sqrt(2)*sigma*g_he$nodes))/sqrt(pi)
   
   return(log(n_pred))
}

n_finder <- function(x, res, U=0.5, st=TRUE){
  model <- res$model
  p <- res$par
  GL <- res$GL
  g <- res$g
  M <- res$M
 
 (x - n_est(x, p, M, GL, g, U=U, stochastic=st, model=model))^2
}

sy <- function(res, U=0.5, st=TRUE, n_range=c(2,12)) {
  n_eq <- optimize(n_finder, n_range, res=res, U=U, st=st)
  res <- log(U)+n_eq$minimum
  attr(res, "obj") <- n_eq$objective
  res
}

msy_est <- function(
  res,
  ST=TRUE,
  n_range=c(2,12),
  U_range=c(0.01,0.99)
){
  model <- res$model
  M <- res$M
  p <- res$par
  
  U_msy <- optimize(sy, U_range, maximum=TRUE, res=res, st=ST, n_range=n_range)
  n_msy <- optimize(n_finder, n_range, res=res, U=U_msy$maximum, st=ST)
  S_msy <- exp(n_msy$minimum)*exp(-M*7*26)*(1-U_msy$maximum)

  list(model=model, model_name=res$model_name, p=p, U_msy=U_msy, n_msy=n_msy, S_msy=S_msy, MSY=U_msy$maximum*exp(n_msy$minimum))
}

##

sim_est <- function(res, U=0.5, d=1, Sim=3000, Y=100, seed0=1){
  set.seed(seed0)
  model <- res$model
  p <- res$par
  M <- res$M
  GL <- res$GL
  g <- res$g
  n0 <- p$log_a+p$log_b
  sigma <- exp(p$log_sigma)
  
  n_vec <- matrix(NA, nrow=Y, ncol=Sim)

  z <- matrix(rnorm(Y*Sim), nrow=Y, ncol=Sim)

  n_vec[1,] <- n0+sigma*z[1,]
  
  for (i in 1:(Y-1)){
    n_vec[i+1,] <- log(pred_R(exp(n_vec[i,]-M*7*26+log(1-U)),p,GL,g,model))+sigma*z[i+1,]
  }
  
  exp(n_vec+log(U)-d*sigma^2/2)
}

##

msy_sim_est <- function(
  Sim_dat,
  res,
  ST=TRUE,
  n_range=c(5,10),
  U_range=c(0.01,0.99)
){
  model <- Sim_dat$model
  M <- res$M
  GL <- res$GL
  g <- res$g
  
  m <- Sim_dat$m
  
  res1 <- list(model=model, M=M, GL=GL, g=g)
  
  Res <- list()
  for (i in 1:m){
    p <- as.list(Sim_dat$parms[i,])
    p$a <- exp(p$log_a)
    p$b <- exp(p$log_b)
    res1$par <- p
  
    U_msy <- optimize(sy, U_range, maximum=TRUE, res=res1, st=ST, n_range=n_range)
    n_msy <- optimize(n_finder, n_range, res=res1, U=U_msy$maximum, st=ST)
    S_msy <- exp(n_msy$minimum)*exp(-M*7*26)*(1-U_msy$maximum)

    Res[[i]] <- list(model=model, model_name=res$model_name, p=p, U_msy=U_msy, n_msy=n_msy, S_msy=S_msy, MSY=U_msy$maximum*exp(n_msy$minimum))
  }
  
  Res
}

## Convergence Comparison

comp_table <- function(Res0, Res1){
  table(sapply(1:100, function(i) Res0[[i]]$conv_diag[4]), sapply(1:100, function(i) Res1[[i]]$conv_diag[4]))
}
