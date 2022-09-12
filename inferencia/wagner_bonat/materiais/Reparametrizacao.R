## Exemplo: Reparametrização ---------------------------------------------------
## Autor: Prof. Wagner Hugo Bonat LEG/UFPR -------------------------------------
## CE085 - Estatística Inferencial ---------------------------------------------

## Dados
dados <- c(0.19,1.68,2.81,0.59,1.18)
theta.est <- length(dados)/sum(dados^2)
mu.est <- sum(dados^2)/length(dados)
c(theta.est, mu.est)

## Função de verossimilhança
# Verossimilhança para theta
L.theta <- function(theta, y){
  n <- length(y)
  return((2 * theta)^n * prod(y) * exp(-theta*sum(y^2)))
}

# Verossimilhança para mu 
L.mu <- function(mu, y){
  n <- length(y)
  return((2*mu^-1)^n * prod(y) * exp( -(1/mu)*sum(y^2)))
}

## Graficamente
LLMAX <- L.theta(theta.est, y=dados)
LLMIN <- L.theta(0.045, y=dados)
require(rootSolve)
grid.theta <- seq(0.045,0.9,l=100)
grid.mu <- seq(1/0.90, 1/0.045,l=100)

ll.mu <- apply(as.matrix(grid.mu), 1, L.mu,y=dados)
ll.theta <- apply(as.matrix(grid.theta),1,L.theta, y=dados)

par(mfrow = c(2,2), mar=c(2.6, 2.8, 1.2, 0.5), mgp = c(1.6, 0.6, 0))
plot(ll.theta ~ grid.theta, type="l", ylab=expression(L(theta)), xlab=expression(theta))
segments(theta.est, LLMIN, theta.est, LLMAX)
plot(ll.mu ~ grid.mu,type="l", ylab=expression(L(mu)), xlab=expression(mu))
segments(mu.est, LLMIN, mu.est, LLMAX)

plot(log(ll.theta) ~ grid.theta, type="l", ylab=expression(l(theta)), xlab=expression(theta))
segments(theta.est, log(LLMIN), theta.est, log(LLMAX))
plot(log(ll.mu) ~ grid.mu,type="l", ylab=expression(l(mu)), xlab=expression(mu))
segments(mu.est, log(LLMIN), mu.est, log(LLMAX))


## Deviance
dev.theta <- function(theta, theta.hat, y, desloca=0){
  saida <- 2*(length(y)*log(theta.hat/theta) - 
                (theta.hat - theta)*sum(y^2))
  return(saida - desloca)
}

dev.mu <- function(mu, mu.hat, y, desloca=0){
  saida <- 2*(length(y)*log(mu/mu.hat) + 
                ((1/mu)-mu.hat^-1)*sum(y^2))
  return(saida - desloca)
}

dev.app.theta <- function(theta, theta.hat, y){
  return((theta - theta.hat)^2 * (length(y)/theta.hat^2))
}

dev.app.mu.obs <- function(mu, mu.hat, y){
  Io <- -((length(y)/mu.hat^2) - (2*sum(y^2)/mu.hat^3))
  return((mu - mu.hat)^2 * Io)
}

dev.app.mu.esp <- function(mu, mu.hat, y){
  Ie <- length(y)/(mu.hat^2)
  return((mu - mu.hat)^2 *Ie)
}


## Intervalo assintótico
ic.assintotico <- function(EMV, Io, alpha){
  return(EMV + c(-1, 1) * qnorm(1-alpha/2)*sqrt(Io^(-1)))
}

n <- length(dados)
ic.theta <- ic.assintotico(EMV = theta.est, 
                           Io = (n / theta.est^2), alpha=0.05)
ic.mu.obs <- ic.assintotico(EMV = mu.est, 
                            Io = -n/mu.est^2 + (2*sum(dados^2))/(mu.est^3), alpha=0.05)
ic.mu.esp <- ic.assintotico(EMV = mu.est, Io = n/mu.est^2, alpha=0.05)

## IC deviance
ic.deviance <- function(dev,intervalo,...){
  ic <- uniroot.all(dev, interval=intervalo, ...)
  return(ic)
}

ic.dev.theta <- ic.deviance(dev=dev.theta, intervalo=c(0,10), 
                            theta.hat=theta.est, y=dados, desloca=qchisq(0.95, df=1))

ic.dev.mu <- ic.deviance(dev=dev.mu, intervalo=c(0,100), 
                         mu.hat=mu.est, y=dados, desloca=qchisq(0.95, df=1))

rbind(ic.theta, ic.dev.theta)
rbind(ic.dev.mu, ic.mu.obs,ic.mu.esp)

## Graficamente
d.theta <- apply(as.matrix(grid.theta),1,dev.theta, theta.hat = theta.est, y=dados)
d.theta.app <- apply(as.matrix(grid.theta),1,dev.app.theta, theta.hat = theta.est, y=dados)

d.mu <- apply(as.matrix(grid.mu),1,dev.mu, mu.hat = mu.est, y=dados)
d.mu.obs <- apply(as.matrix(grid.mu),1,dev.app.mu.obs, mu.hat = mu.est, y=dados)
d.mu.esp <- apply(as.matrix(grid.mu),1,dev.app.mu.esp, mu.hat = mu.est, y=dados)

par(mfrow=c(1,2))
plot(d.theta ~grid.theta,type="l",ylab=expression(D(theta)),xlab=expression(theta))
lines(d.theta.app ~grid.theta, col="red", lty=1)
abline(h=qchisq(0.95,df=1), lty=3)
segments(ic.dev.theta, c(0,0), ic.dev.theta, dev.theta(ic.dev.theta, theta.hat = theta.est, y=dados), lty=2)
segments(ic.theta, c(0,0), ic.theta, dev.app.theta(ic.theta, theta.hat = theta.est, y=dados), col=2, lty=2)
legend("topleft",legend=c("Deviance","Aprox. Quadrática"),col=c("black","red"),lty=1)

plot(d.mu ~grid.mu,type="l",ylab=expression(D(mu)),xlab=expression(mu))
lines(d.mu.obs ~ grid.mu, col=2, lty=2)
lines(d.mu.esp ~ grid.mu, col=4, lty=4)
abline(h=qchisq(0.95,df=1), lty=3)
segments(ic.dev.mu, c(0,0), ic.dev.mu, dev.mu(ic.dev.mu, mu.hat = mu.est, y=dados), lty=2)
segments(ic.mu.obs, c(0,0), ic.mu.obs, dev.app.mu.obs(ic.mu.obs, mu.hat = mu.est, y=dados), col=2, lty=2)
segments(ic.mu.esp, c(0,0), ic.mu.esp, dev.app.mu.esp(ic.mu.esp, mu.hat = mu.est, y=dados), col=4, lty=4)
legend("topleft",legend=c("Deviance","Aprox. Quadrática - Observada","Aprox. Quadrática - Esperada"), col=c(1,2,4),lty=c(1,2,4))


## Reparametrização
V.mu <- 1/(theta.est^2 * n)
ic.mu.theta <- c(1/theta.est - qnorm(0.975)*sqrt(V.mu),
                 1/theta.est + qnorm(0.975)*sqrt(V.mu))

cbind(ic.mu.theta, ic.mu.esp)
