## Exemplo: Teste de hipóteses -------------------------------------------------
## Autor: Wagner Hugo Bonat LEG/UFPR -------------------------------------------
## CE085: Estatística Inferencial ----------------------------------------------
rm(list=ls())

# Simulando a amostra
set.seed(123)
x <- rbinom(100, size = 1, prob = 0.3)

# Log-verossimilhança
vero <- function(theta, y){
  ll <- sum(dbinom(y, prob = theta, size = 1, log=TRUE))
  return(ll)
}

# Avaliando a log-verossimilhança no grid
grid.theta <- seq(0.15, 0.45, l = 100)
ll.theta <- apply(as.matrix(grid.theta), 1, vero, y = x)

## Teste da razão de verossimilhança
trv <- function(Est, H0, alpha, ...) {
  critico <- qchisq(1-alpha, df=1)
  est.calc <- Est(H0, ...)
  print(ifelse(est.calc < critico, "Aceita H0", "Rejeita H0"))
  return(c(est.calc,critico))
}

# Teste de Wald
wald <- function(H0, EMV, V.EMV, alpha){
  critico <- qnorm(1-alpha/2)
  Tw <- (EMV - H0)/sqrt(V.EMV)
  print(ifelse(abs(Tw) < critico, "Aceita H0", "Rejeita H0"))
  return(c(Tw,critico))
}

# Teste escore
escore <- function(H0, U, Ie, alpha, ...){
  critico <- qnorm(1-alpha/2)
  Te <- U(H0,...)/sqrt(Ie(H0,...))
  print(ifelse(abs(Te) < critico, "Aceita H0", "Rejeita H0"))
  return(c(Te,critico))
}


# Simulando dados
set.seed(123)
dados <- rbinom(100, prob = 0.3, size = 1)

## Estatistica do TRV caso Bernoulli
Est <- function(H0, y){
  theta.hat <- mean(y)
  n <- length(y)
  lv <- -2*(sum(y)*log(H0/theta.hat) + (n - sum(y))*log( (1-H0)/(1-theta.hat) ) )
  return(lv)
}

## Procedendo com o TRV theta0 = 0.4
trv(Est = Est, H0 = 0.4, alpha = 0.05, y = dados)

## Teste de Wald
V.EMV <- mean(dados)*(1-mean(dados))/length(dados)
wald(H0 = 0.4, EMV = mean(dados), V.EMV = V.EMV, alpha = 0.05)

## Função escore
U <- function(theta, y) {
  n <- length(y)
  esco <- (sum(y)/theta) -  ((n-sum(y))/(1-theta))
  return(esco)
}

# Informação esperada
Ie <- function(theta, y) {
  n <- length(y)
  IE <- n/(theta*(1-theta))
  return(IE)
}

escore(H0 = 0.4, U = U, Ie = Ie, alpha = 0.05, y = dados)
