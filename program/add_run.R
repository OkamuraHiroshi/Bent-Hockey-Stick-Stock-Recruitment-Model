## Additional Run

km_range <- c(0.1,0.5,1,3,5,7,10,15,20,30)

Res2 <- list()
for (i in 1:length(km_range)){
  Res2[[i]] <- leslie(dat, model=2, do_compile=FALSE, pre_process=FALSE, x_int=c(0,km_range[i]))
}

S_X <- seq(0,500)

dat_XX <- NULL

for (i in 1:length(km_range)){
  dat_XX <- rbind(dat_XX,
  data.frame(km=km_range[i], S=S_X, R=pred_R(S_X, Res2[[i]]$par, res2$GL, model=2)))
}

dat_XX$km <- factor(dat_XX$km)

pp0 <- ggplot(subset(dat_XX, !(km %in% c(0.1,3))), aes(x=S, y=R, color=km))+geom_line()+geom_line(data=subset(dat_XX,km==0.1),aes(x=S, y=R), color=gray(0.3), linetype="dotted", linewidth=1)+geom_line(data=subset(dat_XX,km==3),aes(x=S, y=R), color=gray(0.3), linetype="dotted", linewidth=1)+annotate("text", x=300,   y=650, label="km=0.1",color=gray(0.3))+annotate("text", x=0,   y=4000, label="km=3",color=gray(0.3))+theme_bw()

U_msy <- sapply(1:length(km_range), function(i) msy_est(Res2[[i]])$U_msy$maximum)
S_msy <- sapply(1:length(km_range), function(i) msy_est(Res2[[i]])$S_msy)
Obj <- sapply(1:length(km_range), function(i) Res2[[i]]$mod$objective)
Conv <- sapply(1:length(km_range), function(i) Res2[[i]]$conv_diag[5])

dat_MSY <- data.frame(km=km_range, Convergence=Conv, Obj=Obj-min(Obj), U_msy=U_msy, S_msy=S_msy)

pp1 <- ggplot(dat_MSY, aes(x=km, y=U_msy, color=Convergence, shape=Convergence))+geom_point(size=3)+geom_line(data=subset(dat_MSY,!(km %in% c(0.1,3))),aes(x=km,y=U_msy))+labs(y=expression(F[MSY]))+  
theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
        )+theme_bw()
pp2 <- ggplot(dat_MSY, aes(x=km, y=S_msy, color=Convergence, shape=Convergence))+geom_point(size=3)+geom_line(data=subset(dat_MSY,!(km %in% c(0.1,3))),aes(x=km,y=S_msy))+labs(y=expression(S[MSY]))+ 
theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
        )+theme_bw()
        
pp3 <- ggplot(subset(dat_MSY, km>0.1), aes(x=km, y=Obj, color=Convergence, shape=Convergence))+geom_point(size=3)+geom_line(data=subset(dat_MSY,!(km %in% c(0.1,3))),aes(x=km,y=Obj))+labs(y="Difference in Negative Log-Likelihood")+theme_bw()

###

maxY <- max(dat$Year)

Error <- NULL

n <- 4
for (k in 1:length(km_range)){
res21 <- Res2[[k]]
error <- NULL
for (i in 1:n){
res20 <- leslie(subset(dat, Year <= maxY-i), model=2, do_compile=FALSE, pre_process=TRUE, x_int=c(0,km_range[k]))

error <- c(error, log((log(pred_R(res20$ts_out[14-i,1],res20$par,res20$GL,model=2))-log(res21$ts_out[14-(i-1),2]))^2))
}
Error <- cbind(Error, error)
}
colnames(Error) <- km_range
RF <- exp(colMeans(Error))

dat_RF <- data.frame(km=km_range, RF=RF, Convergence=Conv)
        
pp4 <- ggplot(dat_RF, aes(x=km, y=RF, color=Convergence, shape=Convergence))+geom_point(size=3)+geom_line(data=subset(dat_RF,!(km %in% c(0.1,3))),aes(x=km,y=RF))+labs(y="Retrospective Forecasting Error")+theme_bw()

