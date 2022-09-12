## Exemplo: Propriedades do EMV ------------------------------------------------
## Autor: Wagner Hugo Bonat LEG/UFPR -------------------------------------------
## CE085 - Estatística Inferencial ---------------------------------------------


# Modelo Poisson

# Simulando dados
set.seed(20)
y <- rpois(20, lambda = 10)


## Diferentes formas de implementar a função de verossimilhança
veroPois <- function(par, dados, tipo, maxlogL) {
  tipo = match.arg(tipo, choices=c("L","LR","logL","dev"))
  ll <- sapply(par, function(p) sum(dpois(dados, lambda=p, log=TRUE)))
  return(switch(tipo, "L" = exp(ll),
                "LR" = exp(ll-maxlogL),
                "logL" = ll,
                "dev" = 2*(maxlogL-ll)))
}


## Gráfico só para testar a função implementada
mll <- sum(dpois(y, lambda = mean(y), log = TRUE))
curve(veroPois(x, dados = y, tipo = "dev", maxlogL = mll), 8, 11, 
      ylab=expression(D(theta)), xlab=expression(theta))

curve(veroPois(x, dados = y, tipo = "LR", maxlogL = mll), 8, 11, 
      ylab=expression(D(theta)), xlab=expression(theta))


## Diferentes representações da função de verossimilhança
par(mfrow = c(1, 4), mar=c(2.6, 2.6, 1.2, 0.5), mgp = c(1.6, 0.6, 0))
mll <- sum(dpois(y, lambda=mean(y), log=TRUE))
curve(veroPois(x, dados=y, tipo="L", maxlogL=mll), 8, 11, 
      ylab=expression(L(theta)), xlab=expression(theta))
curve(veroPois(x, dados=y, tipo="LR", maxlogL=mll), 8, 11,
      ylab=expression(LR(theta)),xlab=expression(theta))
curve(veroPois(x, dados=y, tipo="logL", maxlogL=mll), 8, 11,
      ylab=expression(l(theta)),xlab=expression(theta))
curve(veroPois(x, dados=y, tipo="dev", maxlogL=mll), 8, 11,
      ylab=expression(D(theta)),xlab=expression(theta))



# Estudo de simulação para ilustrar as propriedades assintóticas do EMV
par(mfrow = c(1, 4), mar=c(2.6, 2.6, 1.2, 0.5), mgp = c(1.6, 0.6, 0))
set.seed(123)
emv5 <- replicate(1000, {mean(rpois(5, lam=10))})
emv50 <- replicate(1000, {mean(rpois(50, lam=10))})
emv100 <- replicate(1000, {mean(rpois(100, lam=10))})
emv1000 <- replicate(1000, {mean(rpois(1000, lam=10))})

## Histogramas (distribuição empírica)
hist(emv5,prob=TRUE,ylab="Densidade",xlab=expression(hat(theta)), 
     main=" n = 5", ylim=c(0, 0.3), xlim = c(6, 14.5))
media = 10 ; sd <- sqrt(media/5) ## usando I_E ...
curve(dnorm(x,media,sd), 5, 17, add=TRUE)
abline(v=mean(emv5),col="red")
text(x=mean(emv5) , y= 0.29, label=substitute(avg(hat(theta)) == a, list(a=mean(emv5))), pos=4)
#
hist(emv50,prob=TRUE,ylab="Densidade", xlab=expression(hat(theta)), 
     main=" n = 50",ylim=c(0,0.90), xlim = c(6,14.5))
sd <- sqrt(10^2/(10*50))
curve(dnorm(x,media,sd),5,15, add=TRUE)
abline(v=mean(emv50),col="red")
text(x=mean(emv50), y= 0.895, label=substitute(avg(hat(theta)) == a, list(a=mean(emv50))), pos=4)
#
hist(emv100,prob=TRUE, ylab="Densidade",xlab=expression(hat(lambda)), 
     main=" n = 100",ylim=c(0,1.3), xlim = c(6, 14.5))
sd <- sqrt(10^2/(10*100))
curve(dnorm(x,media,sd),5,15, add=TRUE)
abline(v=mean(emv100),col="red")
text(x= mean(emv100), y= 1.28, label=substitute(avg(hat(theta)) == a, list(a=mean(emv100))), pos=4)
#
hist(emv1000,prob=TRUE, ylab="Densidade",xlab=expression(hat(lambda)), 
     main=" n = 1000",ylim=c(0,1.3), xlim = c(6, 14.5))
sd <- sqrt(10^2/(10*1000))
curve(dnorm(x,media,sd),5,15, add=TRUE)
abline(v=mean(emv1000),col="red")
text(x= mean(emv1000), y= 1.28, label=substitute(avg(hat(theta)) == a, list(a=mean(emv1000))), pos=4)

## FIM --------------------------------------------------------------------------------------