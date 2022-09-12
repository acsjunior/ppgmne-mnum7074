## Exemplo verossimilhança perfilhada: Distribuição Gama -----------------------
## Autor: Prof. Wagner Hugo Bonat LEG/UFPR -------------------------------------
## CE085: Estatística Inferencial ----------------------------------------------

## Simulando dados da distribuição Gama Yi ~ Ga(a, s) --------------------------


## Parametrização a = shape e s = scale
a = 20
s = 10
Ey <- a*s
Vy <- a*s^2

## Dados simulados -------------------------------------------------------------
set.seed(123)
y <- rgamma(n = 100, shape = a, scale = s)
hist(y)

## Estratégia 1: Usar um algoritmo numérico para maximizar a log-verossimilhança

# Log-verossimilhança
ll_gamma <- function(par, y) {
  out <- sum(dgamma(x = y, shape = par[1], scale = par[2], log = TRUE))
  return(out)
}

# Avaliando a log-verossimilhança em um ponto
ll_gamma(par = c(a, s), y = y)

# Maximizando a log-verossimilhança numéricamente
fit1 <- optim(par = c(a,s), fn = ll_gamma, method = "Nelder-Mead", 
              hessian = TRUE, control = list(fnscale = -1), y = y)

## Estimativas
fit1$par

# Matriz de informação observada aproximada numéricamente
Io <- -fit1$hessian
Io

# Matriz de variância-covariancia assintótica
V_theta <- solve(Io)
V_theta

# Intervalo de confiança assintótico (Wald)
lim_inf <- fit1$par - qnorm(0.975)*sqrt(diag(V_theta))
lim_sup <- fit1$par + qnorm(0.975)*sqrt(diag(V_theta))
cbind(lim_inf, fit1$par, lim_sup)

# Estratégia 2: Usar o pacote bbmle --------------------------------------------
require(bbmle)

# Reescrevendo a log-verossimilhança no formato que o pacote pede
ll_gamma_bbmle <- function(a, s, y) {
  out <- sum(dgamma(x = y, shape = a, scale = s, log = TRUE))
  return(-out)
}

fit_bbmle <- mle2(ll_gamma_bbmle, 
                  start = list("a" = a, "s" = s), 
                  data = list("y" = y))

## Resumo do ajuste
summary(fit_bbmle)

## Intervalo de confiança de Wald (ligeira diferença na aproximação do Hessiano)
confint(fit_bbmle, method = "quad")

## Intervalo de confiança baseado em perfil de verossimilhança
confint(fit_bbmle)

## Gráfico do perfil de verossimilhança
prof <- profile(fit_bbmle)
plot(prof)

## Estratégia 3: Implementar tudo na "mão" ------------------------------------

# Log-verossimilhança
ll_gamma <- function(par, y) {
  out <- sum(dgamma(x = y, shape = par[1], scale = par[2], log = TRUE))
  return(out)
}

# Função escore
fc_escore <- function(par, y) {
  n <- length(y)
  a <- par[1]
  s <- par[2]
  U_a <- -n*log(s) - n*digamma(a) + sum(log(y))
  U_s <- - (n*a)/s + (1/s^2)*sum(y)
  esc <- c(U_a, U_s)
  return(esc)
}

# Resolvendo a escore numéricamente
library("nleqslv")
fit3 <- nleqslv(c(a,s), fc_escore, y = y, method = "Newton")

# Estimativas
fit3$x

## Informação esperada
Ie <- function(par, y) {
  n <- length(y)
  a <- par[1]
  s <- par[2]
  mat <- matrix(NA, 2, 2)
  mat[1,1] <- n*trigamma(a)
  mat[1,2] <- n/s
  mat[2,1] <- n/s
  mat[2,2] <- n*a/(s^2)
  return(mat)
}

# Variancia e covariancia
solve(Ie(par = fit3$x, y = y))

# Comparando com a baseada na informação observada
V_theta

# Intervalo de Wald
V_theta_E <- solve(Ie(par = fit3$x, y = y))
lim_inf_E <- fit1$par - qnorm(0.975)*sqrt(diag(V_theta_E))
lim_sup_E <- fit1$par + qnorm(0.975)*sqrt(diag(V_theta_E))
cbind(lim_inf_E, fit3$x, lim_sup_E)

## Verossimilhança perfilhada para a

perf_a <- function(a, y) {
  optim_s <- function(s, a, y) {ll_gamma(par = c(a, s), y = y)}
  out <- optim(par = 10, fn = optim_s, y = y, a = a, method = "Brent", 
               lower = 0, upper = 100, control = list(fnscale = -1))
  return(out$value)
}

## Vetorizando
perf_a_vec <- Vectorize(FUN = perf_a, vectorize.args = "a")

## Gráfico do perfil de verossimilhança para a
a_grid <- seq(15, 45, l = 100)
ll_perf_a <- perf_a_vec(a = a_grid, y = y)
plot(ll_perf_a ~ a_grid, type = "l")

## Deviance perfilhada para a
dev_perf_a <- function(a, EMV, y) {
  DEV <- -2*(perf_a_vec(a = a, y = y) - ll_gamma(par = EMV, y = y))
  return(DEV)
}

plot(dev_perf_a(a = a_grid, EMV = fit3$x, y = y) ~ a_grid, type = "l")
abline(a = qchisq(0.95, df = 1), b = 0)

## Encontrando o IC numericamente
ic_perf_a <- function(a, EMV, y) {
  out <- dev_perf_a(a = a, EMV = EMV, y = y) - qchisq(0.95, df = 1)
  return(out)
}
uniroot.all(f = ic_perf_a, interval = c(5, 50), y = y, EMV = fit3$x)
