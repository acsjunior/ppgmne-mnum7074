### O modelo:
## Y|theta ~ B(n=16, theta) (verossimilhança)
## theta ~ U[0,1] \equiv B(1,1)      (priori)
##
## MCMC
## Metropolis

# 3 conceitos
#Monte Carlo         : samples/proposal
#Markov Chain        : random walk
#Metropolis Hastings 

## proposta para ir de um valor de \theta_{i} para \theta_{i+1}:
## \theta_{i+1} \sim U(\theta, -qp, qp)

##r(\theta_{i+1}, \theta_{i}) = \frac{f(\theta_{i+1}|y)}{f(\theta_{i}|y}
##     = \frac{f(\theta_{i+1}) f(y|\theta_{i+1})}{f(\theta_{i})f(y|\theta_{i})}

##    if(r \geq 1) aceita novo valor
##    if(r < 1)    aceita novo valor com probabilidade r
##                 na prática, para implementar esta regra: (u \sim U[0,1], aceita se u < r)

##rm(list=ls())

rm(list=ls())

## Os Dados
n <- 16
ydata <- 6

N <- 100000
theta.MCMC <- numeric(N)

## "tunning da proposal"
qp <- 0.1
i <- 1
theta.MCMC[1] <- 0.95

for(i in 2:N){
    newOK <- FALSE
    while(!newOK){
        theta.new <- theta.MCMC[i-1] + runif(1, -qp, qp)
        newOK <- (theta.new  >= 0 & theta.new <= 1)
    }
    ##    acc <- (dunif(theta.new)*dbinom(ydata, n=n, prob=theta.new))/(dunif(theta.new)*dbinom(ydata, n=n, prob=theta.new))
    acc <- (dbinom(ydata, size=n, prob=theta.new))/(dbinom(ydata, size=n, prob=theta.MCMC[i-1]))
    acc <- min(acc, 1)
    aceita <-  runif(1) < acc 
    theta.MCMC[i] <- ifelse(aceita, theta.new, theta.MCMC[i-1])
    ##        print(c(i,acc,theta.MCMC[i]))
    i <- i+1
}

hist(theta.MCMC, freq=FALSE, xlim=c(0,1))
lines(density(theta.MCMC))
curve(dbeta(x,1,1), add=T, col=2)
curve(dbeta(x,7,11), add=T)

plot(theta.MCMC, type="l")
plot(acf(theta.MCMC))

##    Questões
##    - influência do valor inicial (solução: burn-in)
##    - autocorrelação devida a Cadeia de Markov (solução: "thinning" - raleamento)
##    - taxa de aceitação (cada algoritimo tem valores de referẽ nca, para Metrópolis rw por volta de 30%

## burn-in, thining e acc-rate
burnin <- 10000
thinning <- 10
N.mcmc <- 10000
N <- burnin + thinning * N.mcmc
theta.MCMC <- numeric(N.mcmc)

i <- 1 ; j <- 1
acc.conta <- 0
qp <- 0.1
theta.c <- 0.5

for(i in 2:N){
    newOK <- FALSE
    while(!newOK){
        theta.new <- theta.c + runif(1, -qp, qp)
        ##            print(c(theta.c, theta.new))
        newOK <- (theta.new  >= 0 & theta.new <= 1)
    }
    ##    acc <- (dunif(theta.new)*dbinom(ydata, n=n, prob=theta.new))/(dunif(theta.new)*dbinom(ydata, n=n, prob=theta.new))
    acc <- (dbinom(ydata, size=n, prob=theta.new))/(dbinom(ydata, size=n, prob=theta.c))
    acc <- min(acc, 1)
    aceita <-  (runif(1) < acc) 
    theta.c <- ifelse(aceita, theta.new, theta.c)
    acc.conta <- acc.conta + aceita    
    ##               print(c(i,j,acc,aceita, theta.c, theta.new))
    if(i > burnin)
        if(i %% thinning == 0){
            theta.MCMC[j] <- theta.c
            j <- j+1
        }
    i <- i+1
}


hist(theta.MCMC, freq=FALSE, xlim=c(0,1))
lines(density(theta.MCMC))
curve(dbeta(x,1,1), add=T, col=2)
curve(dbeta(x,7,11), add=T)

plot(theta.MCMC, type="l")
plot(acf(theta.MCMC))

## calibrar a aceitação para 30%

##
## Jags
##

datalist <- dump.format(list(x = 6, n = 16))
params <- c("theta")
inicial <- dump.format(list(theta = 0.5))

## Modelo
mod <- "model{
x ~ dbin(theta, n)
theta ~ dbeta(1, 1)
}"

## Ajuste
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = c(inicial, inicial),
    n.chains = 2, burnin = 5000, thin = 5, sample = 10000
)

## Resultados
m.jags
qbeta(c(0.025, 0.5, 0.975), 1 + 6, 1 + 16 - 6)
plot(m.jags)

