## Simulation_bhs
##

library(TMB)
library(mvtnorm)
library(tidyverse)

compile("leslie.cpp")
dyn.load(dynlib("leslie"))

## data

dat <- read.csv("squid.csv")
dat <- dat %>% mutate(Cat=Catch/Mass, CR=CPUE/Mass)
TC <- dat %>% group_by(Year) %>% summarize(TC=sum(Cat))
res <- lapply(1987:2000, function(i) lm(CR~cumsum(Cat),data=subset(dat, Year==i)))
lm_res <- sapply(1:length(res), function(i) c(-res[[i]]$coef[2], -res[[i]]$coef[1]/res[[i]]$coef[2]))
lm_res <- rbind(lm_res, TC$TC)
lm_res <- as.data.frame(t(lm_res))
names(lm_res) <- c("q","N0","Catch")
rownames(lm_res) <- 1987:2000

source("leslie.r")
source("sim.r")

res0 <- leslie(dat, model=0, do_compile=FALSE, pre_process=FALSE)

Sim_dat0 <- Simulated_data(res0,Sim=1000)

Res_sim02 <- sim2est(Sim_dat0, model=2, mod="HS")
Res_sim01 <- sim2est(Sim_dat0, model=1, mod="HS")
Res_sim00 <- sim2est(Sim_dat0, model=0, mod="HS")
Res_sim03 <- sim2est(Sim_dat0, model=3, mod="HS")
Res_sim04 <- sim2est(Sim_dat0, model=4, mod="HS")

if (0){
res1 <- leslie(dat, model=1, do_compile=FALSE, pre_process=FALSE, g0=0.01)

Sim_dat1 <- Simulated_data(res1,Sim=1000)

Res_sim12 <- sim2est(Sim_dat1, model=2, mod="MR")
Res_sim11 <- sim2est(Sim_dat1, model=1, mod="MR")
Res_sim10 <- sim2est(Sim_dat1, model=0, mod="MR")
Res_sim13 <- sim2est(Sim_dat1, model=3, mod="MR")
Res_sim14 <- sim2est(Sim_dat1, model=4, mod="MR")
}

res2 <- leslie(dat, model=2, do_compile=FALSE, pre_process=FALSE)

Sim_dat2 <- Simulated_data(res2,Sim=1000,seed0=3)

Res_sim21 <- sim2est(Sim_dat2, model=1, mod="BHS")
Res_sim22 <- sim2est(Sim_dat2, model=2, mod="BHS")
Res_sim20 <- sim2est(Sim_dat2, model=0, mod="BHS")
Res_sim23 <- sim2est(Sim_dat2, model=3, mod="BHS")
Res_sim24 <- sim2est(Sim_dat2, model=4, mod="BHS")

res3 <- leslie(dat, model=3, do_compile=FALSE, pre_process=FALSE)

Sim_dat3 <- Simulated_data(res3,Sim=1000,seed0=5)

Res_sim31 <- sim2est(Sim_dat3, model=1, mod="BH")
Res_sim30 <- sim2est(Sim_dat3, model=0, mod="BH")
Res_sim32 <- sim2est(Sim_dat3, model=2, mod="BH")
Res_sim33 <- sim2est(Sim_dat3, model=3, mod="BH")
Res_sim34 <- sim2est(Sim_dat3, model=4, mod="BH")

res4 <- leslie(dat, model=4, do_compile=FALSE, pre_process=FALSE)

Sim_dat4 <- Simulated_data(res4,Sim=1000,seed0=7)

Res_sim41 <- sim2est(Sim_dat4, model=1, mod="RI")
Res_sim42 <- sim2est(Sim_dat4, model=2, mod="RI")
Res_sim40 <- sim2est(Sim_dat4, model=0, mod="RI")
Res_sim43 <- sim2est(Sim_dat4, model=3, mod="RI")
Res_sim44 <- sim2est(Sim_dat4, model=4, mod="RI")

if (0){
save(res1, Sim_dat1, Res_sim10, Res_sim11, Res_sim12, Res_sim13, Res_sim14, file="Res_s1.rda")
save(res2, Sim_dat2, Res_sim20, Res_sim21, Res_sim22, Res_sim23, Res_sim24, file="Res_s2.rda")
save(res3, Sim_dat3, Res_sim30, Res_sim31, Res_sim32, Res_sim33, Res_sim34, file="Res_s3.rda")
save(res4, Sim_dat4, Res_sim40, Res_sim41, Res_sim42, Res_sim43, Res_sim44, file="Res_s4.rda")
}

p0 <- cowplot::plot_grid(res0$p, res2$p, res3$p, res4$p, labels=c("   Hockey-Stick","Bent Hockey-Stick","   Beverton-Holt","       Ricker"))

##

U_range <- seq(0.01,0.99,by=0.01)

Sim1 <- sapply(U_range, function(u) mean(sim_est(res0, U=u)[100,]))
Ana1 <- sapply(U_range, function(u) exp(sy(res0, U=u)))
Ana10 <- sapply(U_range, function(u) exp(sy(res0, U=u, st=FALSE)))

dat1 <- data.frame(Model="Hockey-Stick", Type=rep(c("Simulation","Stochastic","Deterministic"),each=length(U_range)),U=rep(U_range, 3), SY=c(Sim1, Ana1, Ana10))

Sim2 <- sapply(U_range, function(u) mean(sim_est(res2, U=u)[100,]))
Ana2 <- sapply(U_range, function(u) exp(sy(res2, U=u)))
Ana20 <- sapply(U_range, function(u) exp(sy(res2, U=u, st=FALSE)))
dat1 <- rbind(dat1, 
  data.frame(Model="Bent Hockey-Stick", Type=rep(c("Simulation","Stochastic","Deterministic"),each=length(U_range)),U=rep(U_range, 3), SY=c(Sim2, Ana2, Ana20))
)

Sim3 <- sapply(U_range, function(u) mean(sim_est(res3, U=u)[100,]))
Ana3 <- sapply(U_range, function(u) exp(sy(res3, U=u)))
Ana30 <- sapply(U_range, function(u) exp(sy(res3, U=u, st=FALSE)))
dat1 <- rbind(dat1, 
  data.frame(Model="Beverton-Holt", Type=rep(c("Simulation","Stochastic","Deterministic"),each=length(U_range)),U=rep(U_range, 3), SY=c(Sim3, Ana3, Ana30))
)

Sim4 <- sapply(U_range, function(u) mean(sim_est(res4, U=u)[100,]))
Ana4 <- sapply(U_range, function(u) exp(sy(res4, U=u)))
Ana40 <- sapply(U_range, function(u) exp(sy(res4, U=u, st=FALSE)))
dat1 <- rbind(dat1, 
  data.frame(Model="Ricker", Type=rep(c("Simulation","Stochastic","Deterministic"),each=length(U_range)),U=rep(U_range, 3), SY=c(Sim4, Ana4, Ana40))
)

dat1$Model <- factor(dat1$Model, levels=c("Hockey-Stick","Bent Hockey-Stick","Beverton-Holt","Ricker"))
dat1$Type <- factor(dat1$Type, levels=c("Stochastic","Deterministic","Simulation"))

pSY <- ggplot()+geom_point(data = subset(dat1, Type == "Deterministic"), aes(x = U, y = SY, color = Type, shape = Type, alpha=Type), size = 2.5) + geom_point(data = subset(dat1, Type %in% c("Stochastic", "Simulation")), aes(x = U, y = SY, color = Type, shape = Type,alpha=Type), size = 2.5) + theme_bw() + facet_wrap(~Model, scale="free") + scale_alpha_manual(values = c("Stochastic" = 1, "Deterministic" = 1, "Simulation" = 1))+labs(x="Fishing Rate", y="Sustainable Yield")+
scale_color_manual(values = c("Stochastic" = "#F8766D",
                                "Deterministic" = "#7CAE00",
                                "Simulation" = "#00BFC4")) +
  scale_shape_manual(values = c("Stochastic" = 16,  # 丸
                                "Deterministic" = 17, # 三角
                                "Simulation" = 15))   # 四角

# load(file="Res_s1.rda")

par(mfrow=c(4,2))

true_MSY1 <- msy_sim_est(Sim_dat0, res0)

Res_MSY1 <- NULL
for (i in 1:100){
   print(i)
   temp0 <- msy_est(Res_sim00[[i]])
   temp1 <- msy_est(Res_sim01[[i]])
   temp2 <- msy_est(Res_sim02[[i]])
   temp3 <- msy_est(Res_sim03[[i]])
   temp4 <- msy_est(Res_sim04[[i]])

   U_msy = c(as.numeric(temp0$U_msy$maximum),as.numeric(temp1$U_msy$maximum)
,as.numeric(temp2$U_msy$maximum),as.numeric(temp3$U_msy$maximum),as.numeric(temp4$U_msy$maximum))
   n_msy = c(as.numeric(temp0$n_msy$minimum),as.numeric(temp1$n_msy$minimum)
,as.numeric(temp2$n_msy$minimum),as.numeric(temp3$n_msy$minimum)
,as.numeric(temp4$n_msy$minimum))
   S_msy = c(temp0$S_msy,temp1$S_msy,temp2$S_msy,temp3$S_msy,temp4$S_msy)
   MSY = c(temp0$MSY,temp1$MSY,temp2$MSY,temp3$S_msy,temp4$S_msy)
   Conv = c(Res_sim00[[i]]$conv_diag[5], Res_sim01[[i]]$conv_diag[5], Res_sim02[[i]]$conv_diag[5], Res_sim03[[i]]$conv_diag[5], Res_sim04[[i]]$conv_diag[5])

   Res_MSY1 <- rbind(Res_MSY1, data.frame(TrueModel="HS", Model=c("HS","MR","BHS","BH","RI"),U_msy=U_msy, n_msy=n_msy, S_msy=S_msy, MSY=MSY, Conv=Conv))
}

true_Umsy <- sapply(1:100, function(i) true_MSY1[[i]]$U_msy$maximum)
true_Umsy <- rep(true_Umsy, each=5)
true_Smsy <- sapply(1:100, function(i) true_MSY1[[i]]$S_msy)
true_Smsy <- rep(true_Smsy, each=5)

Res_MSY1 <- Res_MSY1 %>% mutate(BU_msy = (U_msy-true_Umsy)/true_Umsy, BS_msy = (S_msy-true_Smsy)/true_Smsy)

boxplot(BU_msy~Model,data=Res_MSY1,main="True=HS")
abline(h=0,col="red",lty=2)
boxplot(BS_msy~Model,data=Res_MSY1,main="True=HS")
abline(h=0,col="red",lty=2)

true_MSY2 <- msy_sim_est(Sim_dat2, res2)

Res_MSY2 <- NULL
for (i in 1:100){
   print(i)
   temp0 <- msy_est(Res_sim20[[i]])
   temp1 <- msy_est(Res_sim21[[i]])
   temp2 <- msy_est(Res_sim22[[i]])
   temp3 <- msy_est(Res_sim23[[i]])
   temp4 <- msy_est(Res_sim24[[i]])
      
   U_msy = c(as.numeric(temp0$U_msy$maximum),as.numeric(temp1$U_msy$maximum)
,as.numeric(temp2$U_msy$maximum),as.numeric(temp3$U_msy$maximum),as.numeric(temp4$U_msy$maximum))
   n_msy = c(as.numeric(temp0$n_msy$minimum),as.numeric(temp1$n_msy$minimum)
,as.numeric(temp2$n_msy$minimum),as.numeric(temp3$n_msy$minimum)
,as.numeric(temp4$n_msy$minimum))
   S_msy = c(temp0$S_msy,temp1$S_msy,temp2$S_msy,temp3$S_msy,temp4$S_msy)
   MSY = c(temp0$MSY,temp1$MSY,temp2$MSY,temp3$S_msy,temp4$S_msy)
   Conv = c(Res_sim20[[i]]$conv_diag[5], Res_sim21[[i]]$conv_diag[5], Res_sim22[[i]]$conv_diag[5], Res_sim23[[i]]$conv_diag[5], Res_sim24[[i]]$conv_diag[5])

   Res_MSY2 <- rbind(Res_MSY2, data.frame(TrueModel="BHS", Model=c("HS","MR","BHS","BH","RI"),U_msy=U_msy, n_msy=n_msy, S_msy=S_msy, MSY=MSY, Conv=Conv))
}

true_Umsy <- sapply(1:100, function(i) true_MSY2[[i]]$U_msy$maximum)
true_Umsy <- rep(true_Umsy, each=5)
true_Smsy <- sapply(1:100, function(i) true_MSY2[[i]]$S_msy)
true_Smsy <- rep(true_Smsy, each=5)

Res_MSY2 <- Res_MSY2 %>% mutate(BU_msy = (U_msy-true_Umsy)/true_Umsy, BS_msy = (S_msy-true_Smsy)/true_Smsy)

boxplot(BU_msy~Model,data=Res_MSY2,ylab=expression(U[msy]),main="TRUE=BHS")
abline(h=0,col="red",lty=2)
boxplot(BS_msy~Model,data=Res_MSY2,ylab=expression(S[msy]),main="TRUE=BHS")
abline(h=0,col="red",lty=2)

true_MSY3 <- msy_sim_est(Sim_dat3, res3)

Res_MSY3 <- NULL
for (i in 1:100){
   print(i)
   temp0 <- msy_est(Res_sim30[[i]])
   temp1 <- msy_est(Res_sim31[[i]])
   temp2 <- msy_est(Res_sim32[[i]])
   temp3 <- msy_est(Res_sim33[[i]])
   temp4 <- msy_est(Res_sim34[[i]])
      
   U_msy = c(as.numeric(temp0$U_msy$maximum),as.numeric(temp1$U_msy$maximum)
,as.numeric(temp2$U_msy$maximum),as.numeric(temp3$U_msy$maximum),as.numeric(temp4$U_msy$maximum))
   n_msy = c(as.numeric(temp0$n_msy$minimum),as.numeric(temp1$n_msy$minimum)
,as.numeric(temp2$n_msy$minimum),as.numeric(temp3$n_msy$minimum)
,as.numeric(temp4$n_msy$minimum))
   S_msy = c(temp0$S_msy,temp1$S_msy,temp2$S_msy,temp3$S_msy,temp4$S_msy)
   MSY = c(temp0$MSY,temp1$MSY,temp2$MSY,temp3$S_msy,temp4$S_msy)
   Conv = c(Res_sim30[[i]]$conv_diag[5], Res_sim31[[i]]$conv_diag[5], Res_sim32[[i]]$conv_diag[5], Res_sim33[[i]]$conv_diag[5], Res_sim34[[i]]$conv_diag[5])

   Res_MSY3 <- rbind(Res_MSY3, data.frame(TrueModel="BH", Model=c("HS","MR","BHS","BH","RI"),U_msy=U_msy, n_msy=n_msy, S_msy=S_msy, MSY=MSY, Conv=Conv))
}

true_Umsy <- sapply(1:100, function(i) true_MSY3[[i]]$U_msy$maximum)
true_Umsy <- rep(true_Umsy, each=5)
true_Smsy <- sapply(1:100, function(i) true_MSY3[[i]]$S_msy)
true_Smsy <- rep(true_Smsy, each=5)

Res_MSY3 <- Res_MSY3 %>% mutate(BU_msy = (U_msy-true_Umsy)/true_Umsy, BS_msy = (S_msy-true_Smsy)/true_Smsy)

boxplot(BU_msy~Model,data=Res_MSY3,ylab=expression(U[msy]),main="TRUE=BH")
abline(h=0,col="red",lty=2)
boxplot(BS_msy~Model,data=Res_MSY3,ylab=expression(S[msy]),main="TRUE=BH")
abline(h=0,col="red",lty=2)

true_MSY4 <- msy_sim_est(Sim_dat4, res4)

Res_MSY4 <- NULL
for (i in 1:100){
   print(i)
   temp0 <- msy_est(Res_sim40[[i]])
   temp1 <- msy_est(Res_sim41[[i]])
   temp2 <- msy_est(Res_sim42[[i]])
   temp3 <- msy_est(Res_sim43[[i]])
   temp4 <- msy_est(Res_sim44[[i]])
      
   U_msy = c(as.numeric(temp0$U_msy$maximum),as.numeric(temp1$U_msy$maximum)
,as.numeric(temp2$U_msy$maximum),as.numeric(temp3$U_msy$maximum),as.numeric(temp4$U_msy$maximum))
   n_msy = c(as.numeric(temp0$n_msy$minimum),as.numeric(temp1$n_msy$minimum)
,as.numeric(temp2$n_msy$minimum),as.numeric(temp3$n_msy$minimum)
,as.numeric(temp4$n_msy$minimum))
   S_msy = c(temp0$S_msy,temp1$S_msy,temp2$S_msy,temp3$S_msy,temp4$S_msy)
   MSY = c(temp0$MSY,temp1$MSY,temp2$MSY,temp3$S_msy,temp4$S_msy)
   Conv = c(Res_sim40[[i]]$conv_diag[5], Res_sim41[[i]]$conv_diag[5], Res_sim42[[i]]$conv_diag[5], Res_sim43[[i]]$conv_diag[5], Res_sim44[[i]]$conv_diag[5])

   Res_MSY4 <- rbind(Res_MSY4, data.frame(TrueModel="RI", Model=c("HS","MR","BHS","BH","RI"),U_msy=U_msy, n_msy=n_msy, S_msy=S_msy, MSY=MSY, Conv=Conv))
}

true_Umsy <- sapply(1:100, function(i) true_MSY4[[i]]$U_msy$maximum)
true_Umsy <- rep(true_Umsy, each=5)
true_Smsy <- sapply(1:100, function(i) true_MSY4[[i]]$S_msy)
true_Smsy <- rep(true_Smsy, each=5)

Res_MSY4 <- Res_MSY4 %>% mutate(BU_msy = (U_msy-true_Umsy)/true_Umsy, BS_msy = (S_msy-true_Smsy)/true_Smsy)

boxplot(BU_msy~Model,data=Res_MSY4,ylab=expression(U[msy]),main="TRUE=RI")
abline(h=0,col="red",lty=2)
boxplot(BS_msy~Model,data=Res_MSY4,ylab=expression(S[msy]),main="TRUE=RI")
abline(h=0,col="red",lty=2)

Res_MSYs <- rbind(Res_MSY1, Res_MSY2, Res_MSY3, Res_MSY4)

Res_MSYs$Model <- factor(Res_MSYs$Model, levels=c("HS","MR","BHS","BH","RI"))
Res_MSYs$TrueModel <- factor(Res_MSYs$TrueModel, levels=c("HS","BHS","BH","RI"))

pU <- ggplot(Res_MSYs,aes(x=Model, y=BU_msy))+geom_boxplot(fill="skyblue")+geom_hline(yintercept=0, linetype="dashed", colour="red")+theme_bw()+labs(y=expression(U[MSY]))+facet_wrap(~TrueModel, scales="free")
pS <- ggplot(Res_MSYs,aes(x=Model, y=BS_msy))+geom_boxplot(fill="skyblue")+geom_hline(yintercept=0, linetype="dashed", colour="red")+theme_bw()+labs(y=expression(S[MSY]))+facet_wrap(~TrueModel, scales="free")

pU1 <- ggplot(subset(Res_MSYs, Conv),aes(x=Model, y=BU_msy))+geom_boxplot(fill="skyblue")+geom_hline(yintercept=0, linetype="dashed", colour="red")+theme_bw()+labs(y=expression(U[MSY]))+facet_wrap(~TrueModel, scales="free")
pS1 <- ggplot(subset(Res_MSYs, Conv),aes(x=Model, y=BS_msy))+geom_boxplot(fill="skyblue")+geom_hline(yintercept=0, linetype="dashed", colour="red")+theme_bw()+labs(y=expression(S[MSY]))+facet_wrap(~TrueModel, scales="free")
