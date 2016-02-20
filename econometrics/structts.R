
#### Example for R - Structural time series model (linear gaussian)
#### Source: Koopman, S.J. and Durbin J. (2001). Time Series Analysis by State Space Methods. Oxford: Oxford University Press
### Code: http://radhakrishna.typepad.com/TimeSeries_Analysis_State_Space_Methods.pdf

data.1 <- log(read.table("UKdriversKSI.txt",skip=1))
colnames(data.1) <- "logKSI"
data.1 <- ts(data.1, start = c(1969),frequency=12)

model <-SSModel(log(drivers)~SSMtrend(1,Q=list(NA))+
SSMseasonal(period=12,sea.type= 'trigonometric' ,Q=NA),data=Seatbelts,H=NA)

ownupdatefn <- function(pars,model,...){
model$H[] <- exp(pars[1])
diag(model$Q[,,1])<- exp(c(pars[2],rep(pars[3],11)))
model
}

fit<-fitSSM(inits=
log(c(var(log(Seatbelts[, 'drivers' ])),0.001,0.0001)),
model=model,
            updatefn=ownupdatefn,method= 'BFGS' )

res<- data.frame(V = fit$model$H[,,1],
W = diag(fit$model$Q[,,1])[1],
                 Ws = diag(fit$model$Q[,,1])[2])

out <- KFS(fit$model,filtering=c( 'state' ), smoothing=c( 'state' , 'disturbance' , 'mean' ))
level.sm <- out$alphahat[,1]

seas.sm <- rowSums(out$alphahat[,2:12])

epshat <- out$epshat

par(mfrow=c(3,1))

temp <- ts( cbind(data.1,level.sm) ,start = 1969,frequency =12)

plot.ts(temp,plot.type="single", ylim=c(7,8),xlab="",ylab = "log KSI",
        col = c("darkgrey","blue"),lwd=c(1,2),cex.lab=1.3,main = "i")

temp<- ts( seas.sm,start = 1969,frequency = 12)
plot(temp,xlab="",ylab = "log KSI",col = "blue",lwd=1,main="ii")
abline(h=0,col="grey")

temp <- ts( epshat ,start = 1969,frequency =12)
plot(temp,xlab="",ylab = "log KSI",col = "blue",lwd=1,main="iii")
abline(h=0,col="grey")
