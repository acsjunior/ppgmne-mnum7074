## Exemplo: Intervalo de confiança ---------------------------------------------
## Autor: Prof. Wagner Hugo Bonat LEG/UFPR -------------------------------------
## CE085 - Estatística Inferencial ---------------------------------------------


## Simulando dados
set.seed(123) ; y <- rexp(20, rate=1)

grid.theta <- seq(0.68, 2, l=100)

# Verosssimilhança
Lik <- function(theta,y){
  Lt <- prod(dexp(y, rate = theta, log = FALSE))
  return(Lt)
}

L.theta <- sapply(grid.theta, Lik, y=y)
LR.theta <- L.theta/Lik(theta=1/mean(y),y=y)

# log-verossimilhança
ll <- function(theta,y){
  llt <- sum(dexp(y,rate=theta,log=TRUE))
  return(llt)
}

ll.theta <- sapply(grid.theta, ll, y=y)

# Deviance
dev.theta <- 2*(ll(1/mean(y),y=y)-ll.theta)

# Deviance aproximada
dev.aprox <- function(theta,theta.hat,y){
  dev <- (theta - theta.hat)^2*(length(y)/(theta.hat^2))
  return(dev)
}

devA.theta <- sapply(grid.theta, dev.aprox, theta.hat = 1/mean(y), y=y)

par(mfrow=c(1,2))
plot(LR.theta ~ grid.theta,type="l",ylab=expression(LR(theta)),xlab=expression(theta))
abline(h=0.146)

plot(dev.theta ~ grid.theta,type="l",ylab=expression(D(theta)),xlab=expression(theta))
lines(devA.theta ~ grid.theta,type="l",col="red",lty=1);abline(h=3.84)
legend("top",legend=c("Deviance"," Aproximação Quadrática"),col=c("black","red"), lty=c(1,1), bty="n")

theta.est <- round(1/mean(y), dig=2)
theta.est

Ic_min <- theta.est - 1.96*sqrt( (theta.est^2)/20  )
Ic_max <- theta.est + 1.96*sqrt( (theta.est^2)/20  )
c(Ic_min, Ic_max)


## Intervalo Deviance

ICdevExp <- function(theta, theta.hat, y, nivel = 0.95){
  n <- length(y)
  dv <- 2*n*( log( theta.hat/theta) + mean(y)*(theta- theta.hat))
  return(dv - qchisq(nivel,df=1))
}

require(rootSolve)
uniroot.all(ICdevExp, interval = c(0,10), theta.hat = 1/mean(y), y = y)

THETA <- 1
set.seed(12)
ic <- matrix(NA, ncol=2, nrow=100)
for(i in 1:100){
  y <- rexp(70, rate=THETA)
  est <- 1/mean(y)
  ic[i,] <- uniroot.all(ICdevExp, int=c(0,5), theta.hat=est, y=y)
}
mean(apply(ic, 1, function(x) (x[1] < THETA & x[2] > THETA)))

plot(c(0.4,1.8)~c(1,100),type="n",ylab=expression(theta),xlab="Ensaio")
abline(h=1)
for(i in 1:100){
  arrows(c(i),ic[i,1],c(i),ic
         [i,2],code=3,angle=90,length=0.03,
         col=ifelse(ic[i,1] > 1 | ic[i,2] < 1, "darkred","lightgray"))
  # ,      lty=ifelse(ic[i,1] > 1 | ic[i,2] < 1, 1, 2))
}